# cGWASdata
This repository conists all data and code necessery for repoducing all results of cGWAS paper (doi:https://doi.org/10.1101/096982).

## Data desription

The results of the paper was achived using gwas summary data (for more details see paper).
"Data" folder consit gwas summary data, SNP information file and correalation matrices. "Scripts" folder consists functions for applying cGWAS to summary gwas data (desription in file).

### GWAS summary data
All gwas stored in text files with headings:

```
SNP	beta	se	Z	P
rs11224105	0.176564410328865	0.0338255949318409	5.21984641170821	1.79071572930154e-07
```

where SNP is the marker name; beta - effect size;  se - standard error of the beta; Z - Z-value (beta/se); P - P-value.
Sample size should be assumed as N=1785.


Information about alleles, frequencies, imputation quality, phisical postion etc could be found in "SNP_information.txt.zip".
Heading:
```
chr	SNP	pos	A1	A2	freq	R2_impute_info	hw	varg_1785
1	rs11804171	713682	A	T	0.948294173319688	0.911636	0.0944081893439549	0.0947969497504853
1	rs2977670	713754	G	C	0.0517239938021306	0.91039	0.0967611574676493	0.0947078794876727
```
where chr is chromosome; SNP - SNP name; pos - position (r37 assembly); A1 - effective allele; A2 - reference allele; freq - effective allele frequency; R2_imputie_info - imputation quality; hw - Hardy-Wainberg eqilibrium P-value; varg_1785 - exact variance of the SNP (needed for calculation of cGWAS).

- uGWAS:
Consists univariae GWAS results filtered by p-value<1e-6. GC correction wasn't applied.

- GGM-cGWAS:
Consists GGM-cGWAS results filtered by p-value<1e-6. GC correction was applied.

- BN-cGWAS:
Consists BN-cGWAS results filtered by p-value<1e-6. GC correction was applied.
