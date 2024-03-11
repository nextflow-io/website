---
title: "U-DAWS: Developing and maintaining reproducible workflows for bioinformatics data"
date: 2019-05-28
type: col8
tags: nextflow,nfcamp,2019,workshop
status: published
layout: "@layouts/MarkdownPage.astro"
---

## U-DAWS: Developing and maintaining reproducible workflows for bioinformatics data

### Venkat Malladi
*Director Biofinformatics Core Facility, UT Southwestern, USA*


The current effort at UT Southwestern is to enable researchers to run their own reproducible analysis and enable data exploration. Therefore a workflow system requires: (1) that computational resource allocation can be defined for each step; (2) a particular step can be run in parallel across all samples; (3) serial steps can be submitted to a queue or be activated on a cloud node when previous steps complete; (4) steps can be restarted in the case of machine failure; (5) workflow versions can run reproducibly; (6) common steps can be shared between different analysis and (7) visualization tools are available to aid researchers in understanding very complex data.

Here we present UTSW Data Analysis Workflows for Sequencing (U-DAWS), a system we developed to develop scalable reproducible sequence analysis workflows using Git, Nextflow, Containers, Shiny and Azure. Git projects with sub-modules allows us to share common steps, maintain version control and test new development with continuous integration. Nextflow provides the features necessary to run workflows on a high-performance compute cluster.  Shiny provides a versatile method for creating easy to use analysis tools using R, in an easy point and click interface for users. Finally, Azure provides us with the ability to widely distribute and scale our pipelines.

In this talk I will describe our architecture and philosophy of building pipelines. Additionally, I will talk about our development of Astrocyte, our GUI interface for running pipelines. Then, I will discuss our efforts with moving our pipelines to Azure.

### Bio

Venkat Malladi is a computational biologist and software engineer with expertise in high-throughput sequencing analysis and reproducible software development. He is currently working on implementing and distribution of genomic analysis pipelines to facilitate reproducibility. He directs the daily operations for the BICF and works with collaborators to identify appropriate data and analysis and provide biological context to results obtained.

### Registration

To attend Nextflow Camp 2019 register at [this link](https://www.crg.eu/en/event/coursescrg-nextflow-2019).
