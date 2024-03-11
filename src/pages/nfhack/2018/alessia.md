---
title: Nextflow on the go
date: 2018-10-16
tags: nextflow,nfhack18,workshop
layout: "@layouts/MarkdownPage.astro"
---

## Nextflow on the go

### Alessia Visconti
*Research Fellow, King's College London, UK*

[YAMP](https://github.com/alesssia/YAMP) (Yet Another Metagenomics Pipeline) is a workflow that enables the analysis of whole shotgun metagenomic data while using containerization to ensure computational reproducibility and facilitate collaborative research. Being based on Nextflow, YAMP can be executed on any UNIX-like system and offers seamless support for multiple job schedulers as well as for the Amazon AWS cloud.

Many users have already claimed that Nextflow is the most mature system supporting both local, HPC clusters and cloud computing execution, allowing a smooth and problem-free transition from a local development to an HPC/cloud deployment. Yet, we wondered whether Nextflow could also be executed on portable miniaturised high-performance pieces of hardware, that is, on mobile phones. Indeed, if YAMP could be executed on portable devices, this would empower the analysis of metagenomics data from virtually everywhere, with huge benefits for scientists, patients, and clinicians alike.

The test was carried on during the ARM, Atos and Cavium BioData Hackathon challenge, organised at the Wellcome Genome Campus on July 2-3, 2018, where participants were asked to take advantage of fast-moving mobile technologies to develop innovative solutions for biomedical data processing.

During the two-day Hackathon, Team GoGut successfully and easily ported YAMP onto ARM’s latest mobile architecture, the ARM Cavium ThunderX2, showing that the analysis of microbial sequences, and of biological data at large, can be successfully taken out of centralised data centres, and that Nextflow offers a mature solution also for mobile applications.

### Deck

<a href='https://github.com/alesssia/talks/blob/master/NextflowWorkshop18/ViscontiNextflow18.pdf'><img src='/img/deck.png' width='45pt' /></a>

### Bio

[Alessia Visconti](https://www.researchgate.net/profile/Alessia_Visconti) is a research fellow in Computational Biology at the Department of Twin Research and Genetic Epidemiology, King’s College London. Her research activity deals with the study of the genomics and epigenomics of human diseases, with a focus on melanoma, cognition and neurodevelopmental disorders, and IgA nephropathy. Recently, she became interested in the human microbiome and its connections to human health and diseases.

### More information

The event program is available at [this link](https://github.com/nextflow-io/nf-hack18/blob/master/schedule.md). For registration and other information check it out [this page](http://www.crg.eu/en/event/coursescrg-nextflow-reproducible-silico-genomics-0).
