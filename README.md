[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

<p align="center" width="100%">
<img width="40%" src="./cellbridge_logo.png"> 
</p>

# CellBridge

**CellBridge** is an automated and versatile workflow meticulously designed to
simplify and expedite the standard procedures entailed in scRNA-seq analysis,
eliminating the need for specialized bioinformatics expertise. CellBridge
harnesses cutting-edge computational methods, integrating an array of advanced
functionalities. It encompasses various crucial steps in scRNA-seq analysis,
starting from the initial conversion of raw unaligned sequencing reads into the
FASTQ format, followed by read alignment, gene expression quantification,
quality control, doublet removal, normalization, batch correction,
dimensionality reduction, clustering, identification of cell markers, and
accurate cell type annotation. CellBridge provides convenient parameterization
of the workflow, while its Docker-based framework ensures reproducibility of
results across diverse computing environments. 

See [ToBridge](https://github.com/Sanofi-GitHub/PMCB-ToBridge) for the
pre-processing of the data.

<p align="center" width="100%">
<img width="95%" src="./pipeline_schematic.png"> 
</p>

---

The pipeline inputs (required) is one folder containing:

1) A metadata file is required, which should be in `CSV` format and contain
demographic and experimental information. The metadata file can be either
`sample_based` or `cell_based`. NOTE: 
  * one column containing the name of samples is required in both metadata types and must be named `sample`; 
  * for `sample_based` metadata, each row should be associated with one sample; 
  * for `cell_based` metadata, each row should be associated with one cell; 
  * the column with cell-id (barcode) information must be named `cell`; 
  * the column name `sample_id` is reserved for the pipeline.
2) CellBridge accepts different types of input data for analysis. The first type
is the widely used output of the `10X` Genomics Cell Ranger pipeline: the trio of
the matrix of UMI counts, the list of cell barcodes, and the list of gene names.
Additionally, Cellbridge accepts Hierarchical Data Format (`HDF5`, `h5`) file
formats that are generated by the Cell Ranger pipeline. CellBridge also accepts
count matrix files in `txt.gz` format, which contain the scRNA-seq gene expression
quantification information. The feature-barcode matrices can be generated by any
microfluidic-, microwell plates-, or droplet-based scRNA-seq technology.
Finally, CellBridge accepts previously processed Seurat RDS (R Data
Serialization) objects as the input. NOTE:
  * data corresponding to each sample must be provided in separate folders. 
  * the name of each folder must be the same as the sample names provided in the `sample` column of the metadata file.

---

To run the docker locally follow the instructions:
* Before you can run the pipeline, you need to get the source code onto your machine. Clone the cellbridge repository using the following command:
+ `git clone https://github.com/Sanofi-GitHub/PMCB-CellBridge.git` + clone a
pecific `x.y.z` tag (release) with `git clone -b <x.y.z>
https://github.com/Sanofi-GitHub/PMCB-CellBridge.git`
* In order to build the container image, you’ll need to use the Dockerfile. A
Dockerfile is simply a text-based file with no file extension. A Dockerfile
contains a script of instructions that Docker uses to create a container image.
+ change directory to the cellbridge directory you just cloned. 
+ build the pipeline container image: `docker build -t <image_name> .`. 
+ This step may take 30-45 minutes.
* Now that you have an image, you can run the `cellbridge` command in a container. 
To do so, you will use the `docker run` command.
+ to see the pipeline help options use: `docker run --rm -it <image_name>
cellbridge --help`. 
+ sharing files between the host operating system and the
container requires you to bind a directory on the host to the container
mount points using the `-v` argument. There is one available mount points
defined in the container named `data`. 
+ to execute the `cellbridge` command directly inside the container use: 
`docker run -it --rm -v <local-path-to-data>:/data:z <image_name> cellbridge [options ...]`.

In most of the cases, it is hard to tell the optimal parameter values for best
QC in advance. The `only_qc` argument will help users to take a look at the
overall metrics of the data in advance. After examination of all the QC metrics,
user can run the pipeline with the optimal parameter values.

---

As a result, the pipeline produces one `outputs` folder containing three files,
each of which is tagged by a 15-character unique identifier (UI).

1) An HTML report (`<project_name>_cellbridge_v<x.y.z>_<UI>_summary.html`),
containing quality metric plots, tables, and several other plots providing an
overal view of the scRNA-seq data analysis outcomes. 
2) An RDS object (`<project_name>_cellbridge_v<x.y.z>_<UI>_final-object.rds`) 
containing the final seurat object with all accosiated metadata and miscellaneous 
information.
3) An RDS object (`<project_name>_cellbridge_v<x.y.z>_<UI>_middle-object.rds`)
containing all intermediate files required to repreduce the html summary.

---

## Contact

For help and questions please contact the [cellbridge's maintenance team](mailto:nima.nouri@sanofi.com).
