---
title: Shareable and scalable pathogen genomics with Nextflow, containers and Clockwork
date: 2019-05-28
type: col8
tags: nextflow,nfcamp,2019,workshop
status: published
layout: "@layouts/Page.astro"
---

## Shareable and scalable pathogen genomics with Nextflow, containers and Clockwork

### Zamin Iqbal
*Research Group Leader, European Bioinformatics Institute (EMBL-EBI), UK*


Genomics is becoming ubiquitous in microbiology, and there are a plethora of variant analysis tools. I will give a tutorial on our new Clockwork pipeline which we have built for our 100,000 M. tuberculosis genomes project.

There are 4 key components which might be of interest. First: the pipeline is open source, containerised with singularity, and has been successfully been run on Azure, EBI internal cluster, our own laptops, and in training courses in India and Peru. Second: the pipeline includes automated submission to the European Nucleotide Archive (taking metadata in spreadsheet form), removal of human/HIV sequence contamination before submission, variant calling with both samtools and cortex, and state-of-the-art variant adjudication with genome graphs to give high sensitivity and specificity calls, and automated tracking in a database. Third, the pipeline also allows joint genotyping, in order to produce a single "wide" VCF giving consistent genotype calls and likelihoods for all samples at all sites.

I'll talk through how to do this, and also some of the challenges we have seen in our project so far, which has currently reached 50,000 genomes.


### Bio

Zamin Iqbal is a computational biologist at the EBI. He works on fundamental methods for genome analysis (eg cortex, BIGSI) and translation of these methods into clinical and public health practise (e.g drug resistance prediction for TB with Mykrobe).

### Registration

To attend Nextflow Camp 2019 register at [this link](https://www.crg.eu/en/event/coursescrg-nextflow-2019).
