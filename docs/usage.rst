=========================
A note on CORDEX datasets
=========================

This python package has been evolved from collecting reacurring tasks and codes concerning regional climate model datasets, e.g., in the EURO-CORDEX community. However, it is important to note that, although those cmorized regional model datasets usually have the same grid definitions (mostly a rotated latitude longitude grid mapping), there are some pitfalls when many datasets should be compared and analyzed. Some of the confusion when working with CORDEX datasets might come from the fact, that they do not fullfil CF conventionsa as stricly as they should be, especially when datasets from different institutions are to be compared. This poses a problem to CF aware tools, like, e.g., ``xarray`` which is aware of labels and coordinates when datasets are aligned. The documentation and notebooks contained in this package try to address some of the issues.
