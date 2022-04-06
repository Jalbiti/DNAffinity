from itertools import product
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler


class Dataset:
    def __init__(self, protein, fitxer, model_target='tetramers', randomize_fce=False,
                 chosen_features=('full_fce', 'avg'), score='Median_intensity',
                 selected_tetramers=(0, 1, 2, 3, 4, 5, 6, 7, 8), weighted=False):

        self.protein = protein
        self.pbm35 = fitxer

        if weighted:
            self.pbm35 = self.pbm35[["ID_REF", "VALUE", 'WEIGHT']]
            self.pbm35.columns = ["ID_REF", "VALUE", 'WEIGHT']
        else:
            self.pbm35 = self.pbm35[["ID_REF", "VALUE"]]
            self.pbm35.columns = ["ID_REF", "VALUE"]
        self.pbm35 = self.pbm35.dropna()

        inv_seq_data = self.pbm35.copy()
        inv_seq_data["ID_REF"] = inv_seq_data["ID_REF"].apply(self.inv_seq)
        self.pbm35 = pd.concat([self.pbm35, inv_seq_data])
        self.pbm35.reset_index(inplace=True)
        self.electrostatic, self.tetra_avg, self.tetra_fce, self.tetra_fce_reduced, self.onehot_1mer, self.integer, \
            self.presence_tetra, self.mgw = 8 * [None]

        self.feats = chosen_features
        self.score = score
        self.selected_tetramers = selected_tetramers
        self.randomize = randomize_fce
        self.model_target = model_target           # TODO include inverse sequences as well

        # Get list of all kmers for indices purposes
        self.sequences = list(self.pbm35["ID_REF"])
        self.p = len(self.sequences[0])

        # Get list of all tetramers per each octamer in order to get features
        self.featurize()

    def __str__(self):
        return f"{self.model_target}-based dataset for protein {self.protein} using features {self.feats}"

    def featurize(self):
        """
        Computing all features as specified in the constructor,
        and setting them as the respective attributes
        :return: None
        """
        tetra_fce = {line.split()[0]: np.array([float(x) for x in line.split()[1:]])
                     for line in open('fce_tetramer.dat')
                     if 'SHIFT' not in line}
        if 'avg' in self.feats:
            tetra_avg = {line.split()[0]: np.array([float(x) for x in line.split()[1:]])
                         for line in open('avg_tetramer.dat')
                         if 'SHIFT' not in line}
            self.tetra_avg = np.array([np.concatenate([tetra_avg[otmer[i:i+4]] for i in self.selected_tetramers])
                                       for otmer in self.sequences])
        if 'full_fce' in self.feats or 'diagonal_fce' in self.feats:
            if self.randomize:
                keys = list(tetra_fce.keys())
                permut_keys = np.random.permutation(keys)
                tetra_fce = {key: tetra_fce[val] for key, val in zip(keys, permut_keys)}
            self.tetra_fce = np.array([np.concatenate([tetra_fce[otmer[i:i+4]] for i in self.selected_tetramers])
                                       for otmer in self.sequences])
        # We might want to scramble the matchings to do 'negative control',
        # i.e., remove all physical information from the dataset
        # Let's also keep track of the reduced matrix
        if 'diagonal_fce' in self.feats:
            tetra_fce_reduced = {tt: tetra_fce[tt][list(range(0, 36, 7))] for tt in tetra_fce.keys()}
            self.tetra_fce_reduced = np.array([np.concatenate([tetra_fce_reduced[otmer[i:i+4]]
                                                               for i in self.selected_tetramers])
                                               for otmer in self.sequences])
        if 'onehot_1mer' in self.feats:
            # self.onehot_1mer = np.array([self.onehot_encoding(otmer) for otmer in self.sequences]).astype(np.int8)
            self.onehot_1mer = np.array([self.onehot_encoding(otmer) for otmer in self.sequences])
        if 'integer' in self.feats:
            # define encoding input values
            nucleotides = product('ACGT', repeat=4)
            char_to_int = dict((''.join(c), i) for i, c in enumerate(nucleotides))
            self.integer = np.array([[char_to_int[otmer[i:i+4]] for i in self.selected_tetramers]
                                     for otmer in self.sequences]).astype(int)
        if 'presence_tetramer' in self.feats:
            self.presence_tetra = np.array([self.presence(otmer, 4) for otmer in self.sequences]).astype(np.int8)
        #####
        if 'mgw' in self.feats:
            self.mgw = {line.split()[0]: [float(x) for x in line.split()[1:]]
                        for line in open('mgw_rohs.txt')
                        if 'SHIFT' not in line}
            self.mgw = np.array([np.concatenate([self.mgw[otmer[i:i+4]] for i in self.selected_tetramers])
                                 for otmer in self.sequences])

        if 'electrostatic' in self.feats:
            self.electrostatic = {line.split()[0]: [x for x in line.split()[1:]]
                                  for line in open('electrostatic.txt')
                                  if 'SHIFT' not in line}
            self.electrostatic = np.array([np.concatenate([self.electrostatic[otmer[i:i+4]]
                                                           for i in self.selected_tetramers])
                                           for otmer in self.sequences])

    @staticmethod
    def onehot_encoding(sequence):
        """
        Converts a sequence to a one-hot-encoded binary vector
        :param sequence: str, the DNA sequence to encode
        :return: np.array
        """
        # define encoding input values
        nucleotides = 'ACGT'
        # define mapping of chars to integers and viceversa
        char_to_int = dict((c, i) for i, c in enumerate(nucleotides))
        # integer encode input data
        integer_encoded = [char_to_int[char] for char in sequence]
        # one hot encode
        onehot_encoded = list()
        for value in integer_encoded:
            letter = [0 for _ in range(len(nucleotides))]
            letter[value] = 1
            onehot_encoded.extend(letter)
        return np.array(onehot_encoded)

    @staticmethod
    def presence(sequence, k):
        """
        Converts a sequence to a vector of "presence" features that specifies the count of
        a given k-mer occurrences in the sequence; the output vector has length 4**k
        :param sequence: str, the DNA sequence to encode
        :param k: int, length of the k-mers that will be counted
        :return: np.array
        """
        kmers = np.zeros(4**k)
        positions = {''.join(x): n for n, x in enumerate(list(product('ACTG', repeat=k)))}
        for i in range(len(sequence)-k+1):
            kmers[positions[sequence[i:i+k]]] += 1
        return kmers

    @property
    def scores(self):
        """
        Yields the selected target values/affinity scores, normalized to numbers between 0 and 1
        :return: np.array
        """
        keyw = "VALUE" if self.score == 'Median_intensity' else "Z-score"
        if self.model_target == "octamers":
            vals = self.pbm35[keyw].values
            return MinMaxScaler().fit_transform(vals.reshape(-1, 1))
        elif self.model_target == "tetramers":
            mean_scores = self.mean_score()
            vals = np.array([[mean_scores[i][otmer[i:i+4]] for i in self.selected_tetramers]
                             for otmer in self.sequences])
            return MinMaxScaler().fit_transform(vals)

    @staticmethod
    def inv_seq(seq):
        """
        Yields the complementary DNA sequence in the same 5'->3' direction
        :param seq: str, the DNA sequence
        :return: str
        """
        complementary = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return ''.join([complementary[x] for x in seq[::-1]])

    @property
    def features(self):
        """
        Returns the combined matrix of features, pasted from individual pre-calculated attributes
        :return: np.array
        """
        avail = {'avg': 'tetra_avg', 'full_fce': 'tetra_fce', 'diagonal_fce': 'tetra_fce_reduced',
                 'onehot_1mer': 'onehot_1mer', 'integer': 'integer', 'presence_tetramer': 'presence_tetra',
                 'mgw': 'mgw', 'electrostatic': 'electrostatic'}
        return np.hstack([self.__getattribute__(avail[feat]) for feat in self.feats])

    def mean_score(self):
        """
        Generates a list of dictionaries, each containing the position-wise
        score per tetramer
        :return: list of dicts with {str: float} mappings
        """
        keyw = "VALUE" if self.score == 'Median_intensity' else "Z-score"
        if self.model_target == "tetramers":
            expt_scores = self.pbm35[keyw].values
            from collections import defaultdict
            position_scores = [defaultdict(lambda: 0) for _ in range(self.p - 3)]
            position_counts = [defaultdict(lambda: 0) for _ in range(self.p - 3)]
            for seq, score in zip(self.sequences, expt_scores):
                for i in range(self.p - 3):
                    position_scores[i][seq[i:i+4]] += score
                    position_counts[i][seq[i:i+4]] += 1
        else:
            raise NotImplementedError
        return [{ttmer: position_scores[i][ttmer]/position_counts[i][ttmer]
                 for ttmer in position_scores[i].keys()} for i in range(self.p - 3)]
