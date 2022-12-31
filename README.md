[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

<p align="center" width="100%">
<img width="40%" src="./cellbridge_logo.png"> 
</p>

# CellBridge

PMCB scRNA-seq pipeline. The pipeline reads-in multiple scRNA-seq samples, performs QC, merges samples, clusters cells into groups, annotates cells, and identifies gene markers.

The pipeline inputs (required) is one folder containing:

A metadata file (a comma separated file in CSV format) containing demographic and experimental information. The metadata file can be either 'sample_based' or 'cell_based'. NOTEs: one column containing the name of samples is required in both metadata types and must be named 'sample'; for 'sample_based' metadata, each row should be associated with one sample; for 'cell_based' metadata, each row should be associated with one cell; the column with cell id information must be named 'cell'; the column name 'sample_id' is reserved for the pipeline.
Sample data: the standard 10X 'cellranger count' outputs; three barcodes, features, and gene-cell count TAR.GZ files or H5 files per each sample. The data corresponding to each sample must be provided in separate folders. The name of each folder must be the same as the sample names provided in the 'sample' column of the metadata file.

---

Single Cell Biology (SCB)
Precision Medicine and Computational Biology (PMCB) 

---

## Contact

For help and questions please contact the [cellbridge's maintenance team](mailto:nima.nouri@sanofi.com).
