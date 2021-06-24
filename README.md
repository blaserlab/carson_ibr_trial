# Carson Lab Ibrutinib Clincal Trial

* This is the analysis repository for the Carson Lab Ibrutinib Clinical Trial single cell RNA seq data.

* The source data for this repository is available as a pre-compiled binary file at this url:  https://XXX

* The best plan is to download and install interactively via the RStudio packages tab.  The data repository requries the package "lazyData" as a dependencey. This allows pseudo-lazy loading for packages that exceed the builtin file size limits in R.  Lazy-load the data into your analysis R session with: 

```r
library(lazyData)
requireData("carson.ibr.trial.datapkg")
```


