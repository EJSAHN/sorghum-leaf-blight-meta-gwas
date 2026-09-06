# v2.0.2

Adds reproducible OLS transformation/PC diagnostics and matched-estimator
candidate-effect comparisons to the sorghum leaf blight analysis workflow.

The auxiliary command reconstructs the project PC basis, scans raw and
inverse-normal phenotypes under five PC settings, and obtains primary
candidate coefficients using GEMMA mode 3 before comparing them with mode-3
location-omission coefficients. The main location-omission output no longer
compares mode-4 primary coefficients with mode-3 omission coefficients.

Detectable-effect outputs now state residual standard deviations per allele;
columns that incorrectly represented a signal/residual-variance ratio as total
variance explained have been removed. Candidate effect and interval metadata
distinguish descriptive alternative-model estimates from score-test inference.

Primary phenotype construction, SNP QC, kinship construction, GEMMA score/LRT
inference, ACAT/Simes calculations, and exploratory candidate selection are
unchanged. Study data and executable binaries are not included.
