# Single Cell RNA Seq Data Analysis For: A pilot study of the Bruton's tyrosine kinase inhibitor ibrutinib in combination with nivolumab in patients with metastatic solid tumors

* The source data for this analysis will be available at the time of publication as a pre-compiled R package at this DOI:  doi:10.5061/dryad.1c59zw3vs

* Prior to publication, peer reviewers may access the data via the private URL provided in the manuscript.

* The best plan is to clone the analysis project and then download the data package and install interactively via the RStudio packages tab.  This should install the package "lazyData" as a dependency. This allows pseudo-lazy loading for packages that exceed the built-in file size limits in R.  Within the analysis project run: 

```r
lazyData::requireData("carson.ibr.trial.datapkg")
```

* The data objects should be available within a hidden environment and will occupy memory only when called by the analysis code.

* Data pre-processing scripts are available in the data-raw directory of the data package.

