options(java.parameters="-Xmx30G")

library(SELEX, quietly = TRUE)

workDir = file.path(".", "SELEX_workspace")
selex.config(workingDir=workDir,verbose=FALSE, maxThreadNumber= 4)


#selex.defineSample(seqName="seq0", seqFile="./SELEX_workspace/M2B_WT_R0.fastq", sampleName="R0.sample", round=0, varLength=16, leftBarcode="CGC", rightBarcode="CCTGGAA")
#selex.defineSample(seqName="seq1", seqFile="./SELEX_workspace/M2B_WT_R1.fastq", sampleName="R1.sample", round=1, varLength=16, leftBarcode="CGC", rightBarcode="CCTGGAA")
#selex.defineSample(seqName="seq2", seqFile="./SELEX_workspace/M2B_WT_R2.fastq", sampleName="R2.sample", round=2, varLength=16, leftBarcode="CGC", rightBarcode="CCTGGAA")

selex.loadAnnotation("./SELEX_workspace/config_full.xml", "./SELEX_workspace")

selex.sampleSummary()
print("first summary")
#r0=selex.sample(seqName="R0.lib", sampleName="R0.sample", round=0)
#r1=selex.sample(seqName="R1.lib", sampleName="R1.sample", round=1)
#r2=selex.sample(seqName="R2.lib", sampleName="R2.sample", round=2)

r0 = selex.sample(seqName="R0.lib", sampleName="R0.sample", round=0)
r1 = selex.sample(seqName='R1.lib', sampleName='R1.sample', round=1)
r2 = selex.sample(seqName='R2.lib', sampleName='R2.sample', round=2)
r3 = selex.sample(seqName='R3.lib', sampleName='R3.sample', round=3)
r4 = selex.sample(seqName='R4.lib', sampleName='R4.sample', round=4)
r5 = selex.sample(seqName='R5.lib', sampleName='R5.sample', round=5)
r6 = selex.sample(seqName='R6.lib', sampleName='R6.sample', round=6)


#split the r0 sample into testing and training sets
r0.split=selex.split(sample=r0)

r0.split

mm = selex.mm(sample=r0.split$train, order=3, crossValidationSample=r0.split$test)

## get affs
rAff = selex.affinities(sample=r1, minCount=1, k=10, markovModel=mm)
nrow(rAff)
rAff_sorted = rAff[order(rAff$Affinity, decreasing = TRUE), ]
nrow(rAff_sorted)
print(length(unique(rAff_sorted[,1])))
write.table(rAff_sorted, 'six3/six3_0vs1.txt', sep = '\t')

rAff = selex.affinities(sample=r2, minCount=1, k=10, markovModel=mm)
nrow(rAff)
rAff_sorted = rAff[order(rAff$Affinity, decreasing = TRUE), ]
nrow(rAff_sorted)
print(length(unique(rAff_sorted[,1])))
write.table(rAff_sorted, 'six3/six3_0vs2.txt', sep = '\t')

rAff = selex.affinities(sample=r3, minCount=1, k=10, markovModel=mm)
nrow(rAff)
rAff_sorted = rAff[order(rAff$Affinity, decreasing = TRUE), ]
nrow(rAff_sorted)
print(length(unique(rAff_sorted[,1])))
write.table(rAff_sorted, 'six3/six3_0vs3.txt', sep = '\t')

rAff = selex.affinities(sample=r4, minCount=1, k=10, markovModel=mm)
nrow(rAff)
rAff_sorted = rAff[order(rAff$Affinity, decreasing = TRUE), ]
nrow(rAff_sorted)
print(length(unique(rAff_sorted[,1])))
write.table(rAff_sorted, 'six3/six3_0vs4.txt', sep = '\t')

rAff = selex.affinities(sample=r5, minCount=1, k=10, markovModel=mm)
nrow(rAff)
rAff_sorted = rAff[order(rAff$Affinity, decreasing = TRUE), ]
nrow(rAff_sorted)
print(length(unique(rAff_sorted[,1])))
write.table(rAff_sorted, 'six3/six3_0vs5.txt', sep = '\t')

rAff = selex.affinities(sample=r6, minCount=1, k=10, markovModel=mm)
nrow(rAff)
rAff_sorted = rAff[order(rAff$Affinity, decreasing = TRUE), ]
nrow(rAff_sorted)
print(length(unique(rAff_sorted[,1])))
write.table(rAff_sorted, 'six3/six3_0vs6.txt', sep = '\t')
