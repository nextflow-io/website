---
title: FlowCraft: A modular, extensible and flexible tool to build, monitor and report Nextflow pipelines
date: 2018-10-30
tags: nextflow,nfhack18,workshop
layout: "@layouts/Page.astro"
---

## FlowCraft: A modular, extensible and flexible tool to build, monitor and report Nextflow pipelines

### Diogo Nuno Silva

*Postdoctoral researcher, Instituto de Medicina Molecular João Lobo Antunes, Portugal*

Nextflow has emerged as a promising framework for building parallelized, scalable and reproducible workflows using software containers. Amongst its most attractive features are the support for multi-scale containerization and its portability, which allows the same code to be executed on multiple platforms with different batch schedulers and container systems.

Usually, the generation of Nextflow pipelines entails the creation of a main Nextflow file with all the declarations of the pipeline’s processes, its channels and operators, and optional configuration and auxiliary files that are necessary for its execution. However, due to the fast pace and continuous development of software for genomic analyses and the need for a quick response to this change by altering existing pipelines, it would be important to make this process more agile and dynamic. Furthermore, the analyses’ goals also change over time and variations of existing pipelines are needed to answer new questions.

The FlowCraft project aims to address these issues by proposing a new tool that leverages the combination of Nextflow and docker/singularity containers to assemble, monitor and report scientific pipelines. The premise of FlowCraft is simple: the users create a broad set of components written in Nextflow language that can then be freely and easily assembled into ready-to-use pipelines taking full advantage of Nextflow execution. These components are modular pieces of software or scripts (e.g. FastQC or Trimmomatic) that have a set of attributes, including input and output types and other parameters. FlowCraft’s engine then connects them, handling the linking and forking of channels automatically, among several other features.

This modularity has two major advantages: (i) each process on a FlowCraft Nextflow pipeline can be written once and re-used to quickly build fairly complex pipelines for different ends with minimal effort; and (ii) adding/modifying parts of a pipeline does not require modification in the codebase of other components. This means that new additions to any pipeline can be easily introduced, and modifications to existing pipelines only require replacing the corresponding components.

FlowCraft further extends the capabilities of Nextflow by providing an Inspection mode that monitors the execution of a pipeline in real time, either in a terminal using a curses interface, or via the [FlowCraft web application](https://www.youtube.com/watch?v=liG1hlEcs5M). Furthermore, it is currently under development a Report mode that allows the interactive visualization of the results from each component in the web application.

The final aim of this project is to have a flexible, extensible and modular system built on top of Nextflow for building, monitoring and reporting scientific pipelines in any omics field. Flowcraft can be found at [this link](https://github.com/assemblerflow/flowcraft/) and the complete user guide is hosted at [this link](http://flowcraft.readthedocs.io/en/latest/).

### Deck

<a href='https://slides.com/diogosilva-1/nextflow-workshop-2018#/'><img src='/img/deck.png' width='45pt' /></a>


### Bio

Diogo is a freshly minted postdoctoral researcher working on genomics of fungal and bacterial pathogens. He is graduated in Environmental Biology, have a masters in Human and Environmental Biology and a PhD in Bioinformatics and Evolutionary Biology. During his research he have grown a keen fondness for computational biology and all things (bio)informatics. He started hacking on Nextflow as soon as I discovered the framework (about one year ago) and since then he's never looked back on ways of building workflows and scientific pipelines. He is currently active in developing workflows for the analysis of bacterial and fungal genomic data and using Nextflow whenever possible/relevant.


### More information

The event program is available at [this link](https://github.com/nextflow-io/nf-hack18/blob/master/schedule.md). For registration and other information check it out [this page](http://www.crg.eu/en/event/coursescrg-nextflow-reproducible-silico-genomics-0).
