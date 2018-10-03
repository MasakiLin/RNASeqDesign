# RNASeqDesign
This is an R package designed for power calculation and study design in RNA-Seq experiments.

To install this package, please copy and paste following codes into your R session:

1. install.packages("devtools")
2. library(devtools)
3. install_github("MasakiLin/RNASeqDesign")

To reproduce the analysis presented in the manuscript, please download "Data.zip" from following dropbox link and "Codes for reproducing results.zip". 

"Codes for reproducing results.zip":
https://www.dropbox.com/s/rpn3lgv8dpew97c/Codes%20for%20reproducing%20results.zip?dl=0

"Data.zip":
https://www.dropbox.com/s/4jb8l6d7o9na1t9/Data.zip?dl=0

Unzip the zip files will result in two directories "/Data" and "/Codes for reproducing results". 

Please follow the R script "Reproduce results.R" under directory "Codes for reproducing results" to reproduce the analysis. Some data required for the analysis are located under directory "Data".

Another R script "functions.R" contains additional sub-functions used by "Reproduce results.R". Please remember to load in the contents within it.
