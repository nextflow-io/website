title=FA-nf - A Bioinformatics pipeline for functional annotation implemented in Nextflow
date=2019-05-28
type=col8
tags=nextflow,nfcamp,2019,workshop
status=published
~~~~~~

## FA-nf - A Bioinformatics pipeline for functional annotation implemented in Nextflow

### Anna Vlasova
*Bioinformatician, Research Institute of Molecular Pathology (IMP), Austria*

### Toni Hermoso
*Bioinformatician, Centre for Genomics Regulation, Spain*

With the advantages of NGS technologies it became possible to obtain a whole genome sequence and its genome assembly of any novel organism at a relatively low cost and short time. To be able to work with this novel genome assembly, scientists need to know positions of the genic elements, especially protein-coding genes, and their putative function. Therefore, function annotation (FA) is an important step in de-novo genome processing and can provide important information about putative role of concrete gene.  Such annotation usually includes assigning functional domains, i.e. from Pfam or Panther, ontology terms (GO), and specific elements, i.e. cleavage sites. 
 
There are two main outcomes from functional annotation: the first one is an annotation itself, which allows scientists to perform various analysis to understand better genome function. The second one, is an additional quality check for genome assembly and predicted genes. Therefore, it is possible to identify suspicious genes which may belong to another species due to contamination, or non-functional overpredicted genes, erroneously annotated.
 
Here we will present a pipeline for a functional annotation of novel proteins from non-model organisms implemented in Nextflow. The pipeline allows to put together different widely-used tools in the field of functional annotation, including some Java applications and REST API services scripts. This software diversity and complexity is handled thanks to software containers (Docker or Singularity), allowing an easier maintenance and versioning of bundled programs. Data exchange and resulting reports are stored in a database, which can be either sitting on the very filesystem as a single database file using SQLite, or through a preset MySQL DBMS server. For that latter case, we also managed to set up a HPC-compatible MySQL on-demand approach that enabled parallel and subsequent processes of the pipeline to query a single MySQL server instance from different cluster nodes.

The provided output is compatible with commonly-used standard resources for downstream analysis, such as UCSC genome browser or Bioconductor packages. 
This pipeline was used in many genomic projects we collaborated with, among others, those of melon, common bean, wasp, and Iberian lynx.

### Bio 

Anna Vlasova, Bioinformatician at Institute of Molecular Pathology, Vienna, Austria. 
Toni Hermoso Pulido, Bioinformatics Core Facility, Centre for Genomic Regulation, Barcelona. Degree in Biochemistry and PhD in Biotechnology at Autonomous University of Barcelona. Since 2009 he joined CRG as a member of the just established Bioinformatics Core Facility where he has been supporting scientific web services and databases, research training and data analyses at the centre.

### Registration 

To attend Nextflow Camp 2019 register at [this link](https://www.crg.eu/en/event/coursescrg-nextflow-2019).