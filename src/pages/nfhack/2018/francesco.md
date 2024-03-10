---
title: Industrial Personalised Immunotherapy Pipeline Development with Nextflow
date: 2018-10-22
tags: nextflow,nfhack18,workshop
layout: "@layouts/Page.astro"
---

## Agile pipelines with Nextflow: how to go from development to production without pain

### Francesco Strozzi
*Head of Bioinformatics at Enterome Bioscience, France*

Nextflow allows creating pipelines in a very agile way, starting from simple blocks and processes to test new software and parameters. Once a single block of the pipeline works as expected and there are no errors, the developer can simply add a new process and start exploring and testing a new layer of the pipeline. The main Nextflow features allowing this agile development are the data flow model, jobs traceability and the possibility to resume workflows, without the need to re-run unchanged and already completed steps.

The advanced logging and reporting system that Nextflow offers are also critical features helping with jobs profiling and monitoring. Finally, the support for multiple executors and the transparent use of software containers allow to prototype pipelines locally and then to port them in production on large HPC systems or on the cloud.

In this talk, agile pipelines development with Nextflow will be presented from the point of view of Enterome Bioscience, a biotech company working on the human gut microbiome. In the company, cutting-edge bioinformatics pipelines are central to support the development of diagnostics and therapeutics on this innovative field.

The ensemble of Nextflow features allows the Bioinformatics Team to transition without pain from a development stage, where pipelines are designed locally or on very small datasets, to a production phase where full datasets are processed with complete pipeline versions on the AWS cloud. A perspective on pipeline unit testing will be also presented, as a topic of possible broader interest for the Nextflow community.

### Deck

<a href='/misc/nfhack18/francesco.pdf'><img src='/img/deck.png' width='45pt' /></a>

### Bio

Bioinformatics Head at Enterome Bioscience, [Francesco Strozzi](https://www.linkedin.com/in/francescostrozzi/) has more than 14 years of experience in bioinformatics and with the management of bioinformatic teams across genomics, biotechnology and biomedicine fields. In his career he focused on the development of data analysis pipelines, data architecture strategies, and the application of computational methods to support and drive research projects and production platforms. Francesco has more than 6 years of experience with metagenomics and computational microbiology to characterise single microorganisms and complex communities, including the human gut, soil and rumen microbiome.


### More information

The event program is available at [this link](https://github.com/nextflow-io/nf-hack18/blob/master/schedule.md). For registration and other information check it out [this page](http://www.crg.eu/en/event/coursescrg-nextflow-reproducible-silico-genomics-0).
