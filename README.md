# qPCR-processing
R script for streamlined and reproducible qPCR data processing and analysis.
Here I utilize R and the tidyverse package library to take data from qPCR machine outputs to data visualization.

## Libraries Used
'tidyverse'
'readxl'

## qPCR data analysis workflow  
1. Average technical replicates  
2. Calculate dCt by subtracting a reference "house-keeping" Target gene's Ct value from each Ct value, matching the Sample and biological replicate Number  
3. Calculate ddCt by subtracting the average dCt of the reference biological Sample (usually wild-type) from the dCt of each Sample, matching the Target gene  
4. Calculate the Fold Change as 2^ddCt  
