---
title: DolphinNext: A graphical user interface for distributed data processing of high throughput genomics
date: 2019-05-28
type: col8
tags: nextflow,nfcamp,2019,workshop
status: published
layout: "@layouts/Page.astro"
---

## DolphinNext: A graphical user interface for distributed data processing of high throughput genomics

### Alper Kucukural
*Assistant Professor, UMass Medical School, USA*

Emergence of new biomedical technologies, like next-generation sequencing (NGS) which is producing vast amounts of genomic data every day, is driving a big data revolution in biology. The dramatic increase in the volume, as well as the production rate of genomic data, has now made the data analysis new bottleneck for scientific discovery. Naturally, the need for highly-parallel data processing frameworks is greater than ever. It is also important for these frameworks to have certain design characteristics such as flexibility, portability, and reproducibility. Processing of sequencing data usually involves many different programs, each of which performs a specific step in the overall pipeline. Flexibility ensures that the pipelines can support a variety of use cases or data types without the need to modify existing pipelines or create new ones. Portability gives user the freedom to choose computational resources as he/she deems fit. Reproducibility across computing environments, which warrants credibility of the results, is a particularly important feature in the face of the sheer volume of data and complexity of the pipelines.

There exist several platforms that offer graphical user interfaces for designing and execution of complex pipelines (e.g. Galaxy, GenePattern, GeneProf). Unfortunately, none of these platforms supports parallelism or portability across computing environments. To address these and additional shortcomings discussed in this paper, we have created DolphinNext, an easy-to-use graphical user interface for creating and deploying complex workflows for parallel processing of high throughput genomic data. DolphinNext relies on Nextflow which is a framework enabling scalable and reproducible workflows using software containers. The central idea behind the creation of DolphinNext is to facilitate building and deployment of complex pipelines using a graphically-enabled modular approach. DolphinNext provides:

1. A drag and drop user interface that abstracts Nextflow pipelines and allows users to create pipelines without familiarity with Nextflow.
2. Reproducible pipelines by providing version tracking, and by creating a stand-alone version of any pipeline instance to run independently or to share in the publications.
3. Seamless portability to distributed computational environments such as high performance clusters or cloud based solutions.
4. A graphical user interface to monitor pipeline execution that allows restarting of intermediate processes even after parameter changes.

### Bio

Dr. Kucukural designs and implements reusable, robust and production grade bioinformatics analysis pipelines and pipeline generation tools for processing next-generation sequencing data.
He mainly works on NGS data analysis; RNA-Seq, RIP-Seq, Chip-Seq, CLIP-Seq and derivatives. He implemented algorithms to reduce noise by calling the peaks caused by experimental and alignment biases especially for RIP and CLIP-Seq data.

Dr. Kucukural worked on analyzing deep sequencing data to discover key elements of splicing of pre-mRNAs to have better understanding of post-transcriptional regulations of RNAs. Moreover, he has deep knowledge of finding genome wide mRNA targets of RNA binding proteins (RBPs). Analyzing RNA targets of tdp43 RBP with deep sequencing on Rat and human was another focus of his research to understand the mechanisms of neuro-degenerative diseases such as Alzheimer and ALS.

He also worked for protein structure characterization and prediction. He applied techniques from graph theory on protein structure analysis and implemented the theories from computer sciences to biology to find solutions in drug design and small molecular docking fields. He developed applications using genetic algorithms to discover biomarkers and implemented feature detection methods using various clustering, classification and machine learning algorithms such as hidden markov models and support vector machines.

### Registration

To attend Nextflow Camp 2019 register at [this link](https://www.crg.eu/en/event/coursescrg-nextflow-2019).
