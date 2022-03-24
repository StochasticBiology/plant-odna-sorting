# plant-odna-sorting

Statistical analysis for the project "Sorting of plant mitochondrial and plastid heteroplasmy is extremely rapid and depends on MSH1 activity", Broz et al.

Requirements: R with `readxl` and `kimura` packages, the latter of which can be found here https://github.com/lbozhilova/kimura

The wrapper code `analysis.R` invokes the other two files `data-tables.R` (to read in raw data from Excel sheets and cast into usable format) and `kimura-functions.R` (statistical analysis helper functions). It then performs maximum likelihood fits of the Kimura distribution to the sets of observed heteroplasmy, conducts likelihood ratio tests to compare distributions, and does various other calculations.
