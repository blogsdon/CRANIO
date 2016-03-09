#!/bin/sh
#number of cores to reserve for job
nthreads=1

#full s3 path where networks will go
s3="s3://metanetworks/CRANIO/"

#location of data file
dataFile="/shared/CRANIO/cranioRNAseq.csv"

#location of metanetwork synapse scripts
pathv="/shared/metanetworkSynapse/"

#output path for temporary result file prior to pushing to s3/synapse
outputpath="/local/CRANIO/"

#path within s3
s3b="CRANIO"

#id of folder with networks to combine
networkFolderId="syn5650456"

#id of folder on Synapse that file will go to
parentId="syn5650456"

#path to csv file with annotations to add to file on Synapse
annotationFile="/shared/CRANIO/annoFile.txt"

provenanceFile="/shared/CRANIO/provenanceFile.txt"

#path to error output
errorOutput="/shared/CRANIO/Aggregationerror.txt"

#path to out output
outOutput="/shared/CRANIO/Aggregationout.txt"

#job script name
jobname="CRANIOaggregation"

qsub -v s3=$s3,dataFile=$dataFile,pathv=$pathv,outputpath=$outputpath,s3b=$s3b,parentId=$parentId,annotationFile=$annotationFile,provenanceFile=$provenanceFile,networkFolderId=$networkFolderId -pe orte $nthreads -S /bin/bash -V -cwd -N $jobname -e $errorOutput -o $outOutput $pathv/buildConsensus.sh
