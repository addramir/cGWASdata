# cGWASdata

The repository stores all code and data necessary for reproducing all results of cGWAS paper (doi:https://doi.org/10.1101/096982).

## Data description

The results of the paper were obtained from analysis of GWAS summary data (for more details see paper).
"Data" folder stores GWAS summary data, SNP information file and correlation matrices. "Scripts" folder stores R functions for applying cGWAS to summary gwas data (description in the file). / TODO: provide path to the file

### GWAS summary data

All gwas stored in text files with headings:

```
SNP	beta	se	Z	P
rs11224105	0.176564410328865	0.0338255949318409	5.21984641170821	1.79071572930154e-07
```

where SNP is the marker name; beta - effect size;  se - standard error of the beta; Z - Z-value (beta/se); P - P-value.
Sample size should be assumed as N=1785.


Information about alleles, frequencies, imputation quality, physical position etc could be found in "SNP_information.txt.zip".
Heading:
```
chr	SNP	pos	A1	A2	freq	R2_impute_info	hw	varg_1785
1	rs11804171	713682	A	T	0.948294173319688	0.911636	0.0944081893439549	0.0947969497504853
1	rs2977670	713754	G	C	0.0517239938021306	0.91039	0.0967611574676493	0.0947078794876727
```
where chr is chromosome; SNP - SNP name; pos - position (r37 assembly); A1 - effective allele; A2 - reference allele; freq - effective allele frequency; R2_imputie_info - imputation quality; hw - Hardy-Weinberg equilibrium P-value; varg_1785 - exact variance of the SNP (needed for calculation of cGWAS).

### Folders

- uGWAS:
Stores files with univariate GWAS results filtered by P-value<1e-6. Genomic control (GC) correction wasn't applied.

- GGM-cGWAS:
Stores GGM-cGWAS results filtered by P-value<1e-6. GC correction was applied. GC Lambdas are stored in "GGM_cGWAS_gc_lambda.txt" file.

- BN-cGWAS:
Stores BN-cGWAS results filtered by P-value<1e-6. GC correction was applied. GC Lambdas are stored in "BN_cGWAS_gc_lambda.txt" file.

- uGWAS_snps_from_paper:
Stores uGWAS results for all SNPs mentioned in the paper. GC() correction wasn't applied.

- RData: 
Stores the same results for GGM-cGWAS, BN-cGWAS, uGWAS but in RData format suitable for exact_cGWAS_functions.R cGWAS function.


### Matrices

- 20171207_corr_matrix.txt: Pearson correlation matrix for 151 metabolites.

- 20171207_partial_corr_matrix.txt: partial correlation matrix for 151 metabolites. We used "ppcor" R package for calculations.

- 20171207_partial_corr_pvalues_matrix.txt: partial correlation P-value matrix for 151 metabolites.

- 20171207_biochemical_distances.txt: biochemical distances used for BN-cGWAS. This matrix was produced in work of Krumsiek et.al, 2011 (Krumsiek, J., Suhre, K., Illig, T., Adamski, J., & Theis, F. J. (2011). Gaussian graphical modeling reconstructs pathway reactions from high-throughput metabolomics data. BMC Systems Biology, 5(1), 21. https://doi.org/10.1186/1752-0509-5-21).

### Scripts

Scripts are located in "scripts" folder.

- exact_cGWAS_functions.R: the main function for calculation of cGWAS using uGWAS results and correlation matrices. All descriptions are in file. 

- figure2_corrplot.R: in case you want to produce the same beautiful pictures as Figure 2, feel free to use some code of this script ;)

### Citation

If you use the exact_cGWAS_functions.R and the other procedures, respectively,
please cite:

Tsepilov, Y. A., Sharapov, S. Z., Zaytseva, O. O., Krumsiek, J., Prehn, C., Adamski, J., … Aulchenko, Y. S. (2017). Network based conditional genome wide association analysis of human metabolomics. Submitted

/TODO : Yakov please check the reference