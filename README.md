# cGWASdata
This repository conists all data and code necessery for repoducing all results of cGWAS paper (doi:https://doi.org/10.1101/096982).

## Data desription

The results of the paper was achived using gwas summary data (for more details see paper).
"Data" folder consit gwas summary data, SNP information file and correalation matrices. "Scripts" folder consists functions for applying cGWAS to summary gwas data (desription in file).

- uGWAS
Consists univariae GWAS results filtered by p-value<1e-6. GC correction wasn't applied.

- GGM-cGWAS
Consists GGM-cGWAS results filtered by p-value<1e-6. GC correction /b{was} applied.
