#options(java.parameters="-Xmx30G")

library(SELEX, quietly = TRUE)

workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)


#selex.defineSample(seqName="seq0", seqFile="./SELEX_workspace/M2B_WT_R0.fastq", sampleName="R0.sample", round=0, varLength=16, leftBarcode="CGC", rightBarcode="CCTGGAA")
#selex.defineSample(seqName="seq1", seqFile="./SELEX_workspace/M2B_WT_R1.fastq", sampleName="R1.sample", round=1, varLength=16, leftBarcode="CGC", rightBarcode="CCTGGAA")
#selex.defineSample(seqName="seq2", seqFile="./SELEX_workspace/M2B_WT_R2.fastq", sampleName="R2.sample", round=2, varLength=16, leftBarcode="CGC", rightBarcode="CCTGGAA")

selex.loadAnnotation("./SELEX_workspace/config.xml", "./SELEX_workspace")

selex.sampleSummary()
print("first summary")
#r0=selex.sample(seqName="R0.lib", sampleName="R0.sample", round=0)
#r1=selex.sample(seqName="R1.lib", sampleName="R1.sample", round=1)
#r2=selex.sample(seqName="R2.lib", sampleName="R2.sample", round=2)

r0 = selex.sample(seqName="R0.lib", sampleName="R0.sample", round=0)
r2 = selex.sample(seqName='R6.lib', sampleName='R6.sample', round=6)

#split the r0 sample into testing and training sets
r0.split=selex.split(sample=r0)

r0.split

#selex.sampleSummary()
print("second summary")

#Find kmax on the test dataset
k=selex.kmax(sample=r0.split$test)

mm = selex.mm(sample=r0.split$train, order=NA, crossValidationSample=r0.split$test, Kmax=k)

# See Markov model R^2 values
selex.mmSummary()
print("third summary")

## Kmer counting with an offset
#t1 =  selex.counts(sample=r2, k=2, offset=20, markovModel=NULL, outputPath = "counts.txt")

## Kmer counting with a Markov model (produces expected counts)
#t2 =  selex.counts(sample=r2, k=4, markovModel=mm, outputPath = "counts.txt")

# Find the optimal motif length
selex.infogain(sample=r2, markovModel = mm)
ig = selex.infogainSummary()
optimalk = ig$K[which(ig$InformationGain == max(ig$InformationGain))]

# Display all available kmer statistics
table = selex.counts(sample = r2, k = optimalk, markovModel = mm)
head(table)
print("fourth summary")
######
selex.run(trainingSample=r0.split$train, crossValidationSample=r0.split$test,infoGainSample=r2, infoRange=10)

# View all stats
selex.summary()
print("fifth summary")


# Perform the default analysis
#selex.run(trainingSample=r0.split$train, crossValidationSample=r0.split$test,infoGainSample=r2, infoRange=10)

# affinities

mm = selex.mm(sample=r0.split$train, order=3, crossValidationSample=r0.split$test)
r2Aff = selex.affinities(sample=r2, minCount=1,  k=10, markovModel=mm)
nrow(r2Aff)
r2Aff_sorted = r2Aff[order(r2Aff$Affinity, decreasing = TRUE), ]
nrow(r2Aff_sorted)
print(length(unique(r2Aff_sorted[,1])))
write.table(r2Aff_sorted, '6.txt', sep = '\t')

#print("head")
#head(r2Aff_sorted, n = 100)
#print("tail")
#tail(r2Aff_sorted, n = 100)



