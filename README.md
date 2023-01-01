[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

<p align="center" width="100%">
<img width="40%" src="./cellbridge_logo.png"> 
</p>

# CellBridge

PMCB SCB scRNA-seq pipeline. The pipeline reads-in multiple scRNA-seq samples, performs QC per sample, merges samples, clusters cells into groups, dimension reduction, annotates cells, and identifies gene markers.

---

The pipeline inputs (required) is one folder containing:

1) A metadata file (a comma separated file in `CSV` format) containing demographic and experimental information. The metadata file can be either `sample_based` or `cell_based`. NOTEs: 
  * one column containing the name of samples is required in both metadata types and must be named `sample`; 
  * for `sample_based` metadata, each row should be associated with one sample; 
  * for `cell_based` metadata, each row should be associated with one cell; 
  * the column with cell-id information must be named `cell`; 
  * the column name `sample_id` is reserved for the pipeline.
2) Sample data: the standard 10X `cellranger count` outputs; 
  * three barcodes, features, and gene-cell count TAR.GZ files or H5 files per each sample. 
  * data corresponding to each sample must be provided in separate folders. 
  * the name of each folder must be the same as the sample names provided in the `sample` column of the metadata file.

---

To run the docker locally follow the instructions:
* Before you can run the pipeline, you need to get the source code onto your machine. Clone the cellbridge repository using the following command: 
  + `git clone https://<token>@github.com/Sanofi-GitHub/PMCB-CellBridge.git`
  + ask from sargent's maintenance team to get the token.
* In order to build the container image, you’ll need to use the Dockerfile. A Dockerfile is simply a text-based file with no file extension. A Dockerfile contains a script of instructions that Docker uses to create a container image.
  + change directory to the cellbridge directory you just cloned.
  + build the pipeline container image: `docker build -t <image_name> .`. 
  + This step may take 30-45 minutes.
* Now that you have an image, you can run the cellbridge command in a container. To do so, you will use the `docker run` command.
  + to see the pipeline help options use: `docker run --rm -it <image_name> cellbridge --help`.
  + sharing files between the host operating system and the container requires you to bind a directory on the host to one of the container's mount points using the `-v` argument. There is one available mount points defined in the container named `data`.
  + to execute the `cellbridge` command directly inside the container use: `docker run -it --rm -v <local-path-to-data>:/data:z <image_name> cellbridge [options ...]`.

In most of the cases, it is hard to tell the optimal parameter values for best QC in advance. The `only_qc` argument will help users to take a look at the overall metrics of the data in advance. After examination of all the QC metrics, user can run the pipeline with the optimal parameter values.

---

As a result, the pipeline produces one `outputs` folder containing three files and one folder, each of which taged by a 20-characters unique identifier (UI):

* 1) A HTML report (`<project_name>_cellbridge_v<x.y.z>_<UI>_summary.html`), containing quality metric plots, tables, and several other plots providing an overal view of the scRNA-seq data analysis outcomes.
* 2) An RDS object (`<project_name>_cellbridge_v<x.y.z>_<UI>_final-object.rds`) containing the final seurat object with all accosiated metadata and miscellaneous information.
* 3) An RDS object (`<project_name>_cellbridge_v<x.y.z>_<UI>_middle-object.rds`) containing all intermediate files required to repreduce the html summary.
* 4) A folder (`<project_name>_cellbridge_v<x.y.z>_<UI>_spring`) containing required files for SPRING portal visualization.

---

## Contact

For help and questions please contact the [cellbridge's maintenance team](mailto:nima.nouri@sanofi.com).
