import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn import metrics
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Ridge
from sklearn.svm import SVR


class Model:
    regressors = {'random_forest': RandomForestRegressor(n_estimators=10, random_state=None),
                  'linear': LinearRegression(), 'ridge': Ridge(alpha=0.0001), 'svr': SVR()}

    def __init__(self, dataset, len_aln, regressor='random_forest', training_set_size=0.8, weighted=False):
        plt.clf()
        self.data = dataset
        self.X = self.data.features
        self.y = self.data.scores
        if weighted:
            self.weights = self.data.pbm35['WEIGHT']
        else:
            self.weights = None
        self.features = []
        self.l = []
        self.y_pred, self.r2, self.mse = 3 * [None]

        print("features have shape: ", self.X.shape)
        print("target values have shape: ", self.y.shape)
        self.regressor = Model.regressors[regressor]
        # splits w/o randomness for training and testing
        if weighted:
            self.X_train, self.X_test, self.y_train, self.y_test, \
                self.w_train, self.w_test = train_test_split(self.X, self.y, self.weights, train_size=training_set_size,
                                                             random_state=0, shuffle=False)
            self.regressor.fit(self.X_train, self.y_train, self.w_train)
        else:
            self.X_train, self.X_test, self.y_train, \
                self.y_test = train_test_split(self.X, self.y, self.weights, train_size=training_set_size,
                                               random_state=0, shuffle=False)
            self.regressor.fit(self.X_train, self.y_train)

        importances = self.regressor.feature_importances_
        std = np.std([tree.feature_importances_ for tree in self.regressor.estimators_], axis=0)
        indices = np.argsort(importances)[::-1]

        # Print the feature ranking
        print("Feature ranking:")
        for f in range(self.X.shape[1]):
            if indices[f] < 256:
                feat = 'Presence'
            if indices[f] in np.arange(256, 256 + (len_aln - 3) * 6, 1):
                if indices[f] % 6 == 0:
                    feat = 'AVG Shift'
                if indices[f] % 6 == 1:
                    feat = 'AVG Slide'
                if indices[f] % 6 == 2:
                    feat = 'AVG Rise'
                if indices[f] % 6 == 3:
                    feat = 'AVG Tilt'
                if indices[f] % 6 == 4:
                    feat = 'AVG Roll'
                if indices[f] % 6 == 5:
                    feat = 'AVG Twist'
            if indices[f] in np.arange(256 + (len_aln - 3) * 6, 256 + (len_aln - 3) * 12, 1):
                if indices[f] % 6 == 0:
                    feat = 'Diagonal Shift'
                if indices[f] % 6 == 1:
                    feat = 'Diagonal Slide'
                if indices[f] % 6 == 2:
                    feat = 'Diagonal Rise'
                if indices[f] % 6 == 3:
                    feat = 'Diagonal Tilt'
                if indices[f] % 6 == 4:
                    feat = 'Diagonal Roll'
                if indices[f] % 6 == 5:
                    feat = 'Diagonal Twist'
            if indices[f] >= 256 + (len_aln - 3) * 12:
                feat = 'Electro'
            #       for f in range(self.X.shape[1]):
            # print(feat)
            # print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))
            # print(feat,indices[f],importances[indices[f]])
            self.l.append(importances[indices[f]])
            self.features.append([feat, importances[indices[f]]])

    def predict(self, testing_dataset=None):
        if testing_dataset is not None:
            self.X_test, self.y_test = testing_dataset.features, testing_dataset.scores
        self.y_pred = self.regressor.predict(self.X_test).reshape(-1, 1)
        self.y_test = self.y_test.reshape(-1, 1)
        if self.weights is not None:
            self.r2 = np.round(metrics.r2_score(self.y_test, self.y_pred, self.w_test), 3)  # with weights
        else:
            self.r2 = np.round(metrics.r2_score(self.y_test, self.y_pred), 3)           # w/o weights
        self.mse = np.round(metrics.mean_absolute_error(self.y_test, self.y_pred), 6)

    def plot(self):
        plt.rc('xtick', labelsize=15)
        plt.rc('ytick', labelsize=15)
        plt.xlabel('Testing')
        plt.ylabel('Predicted')
        # plt.title(f'MSE = {self.mse}, r$^2$ = {self.r2}')
        plt.xlim = (0, 1)
        plt.ylim = (0, 1)
        plt.scatter(self.y_test, self.y_pred, label=f'r$^2$ = {self.r2}')

    @staticmethod
    def show():
        plt.legend()
        plt.show()
