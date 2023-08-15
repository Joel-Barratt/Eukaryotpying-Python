# This document provides a set of tips and tricks to ensure the code runs properly on your dataset

1) The name of the haplotypes/markers in the haplotype data sheet is **important**. They need to start with "Mt_" or "Nu_". This prefix will automatically tell the code what the ploidy setting should be.
If the marker name starts with the prefic "Mt_" the ploidy for that marker will automatically be set to "1". For markerss starting with "Nu_", the ploidy will automatically be set to "2".

2) Should work best on Python3
   
3) Use pip3 to install pandas, numpy and orderedset dependencies. These should be the only dependencies required.
