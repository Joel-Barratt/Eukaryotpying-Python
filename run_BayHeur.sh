#!/bin/bash -l




euktypDir=/Users/joelbarratt/Documents/NEW_MODULE_3_WORK/PYTHON_DISTANCE_COMP_METHOD/cyclone-master-CYCLONE_MacOS_High_Sierra_BETA_1.001-EUKARYOTYPING/CYCLONE_MacOS_High_Sierra_BETA_1.001/DISTCOMP


#Define input and output dir and file names
pycodeDir=$euktypDir
haplosheetDir="$(dirname "$pycodeDir")"/haplotype_sheets
haplosheetFileName=$(ls $haplosheetDir | sort | tail -1)
haplosheetPathFN=$haplosheetDir/$haplosheetFileName
echo "The newest haplotype file is found:"
echo $haplosheetPathFN
outputDir="$(dirname "$pycodeDir")"/ensemble_matrices
outputFNBase=${haplosheetFileName/.txt/""}
outputhh="${outputDir}/${outputFNBase}_hh_matrixTEST.csv"
outputbh="${outputDir}/${outputFNBase}_bh_matrixTEST.csv"
outputB="${outputDir}/${outputFNBase}_B_matrixTEST.csv"
outputH="${outputDir}/${outputFNBase}_H_matrixTEST.csv"



#Run the command -- first modify "markerList.txt to tell the code which markers are essential" 
python3 $euktypDir/Pycode_distcomp/BayesianHeuristic.py -infile $haplosheetPathFN -outbh $outputbh -outhh $outputhh \
-outB $outputB -outH $outputH -sampmeetlocirequire samplesMeetCutoff.txt \
-expectlocifile $euktypDir/Pycode_distcomp/markerList.txt -expectlocinumber 5

#Immediately the genotypes meeting the minimum data requirements will be printed (i.e., genotypes with any 5 markers (-expectlocinumner 5) within the file "markerList.txt").

echo "Genetic distance computation Complete. "
