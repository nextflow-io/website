title=Learning Nextflow in 2022
date=2022-01-21
type=post
description=Resources for learning Nextflow: 2022 Revision
image=img/learn-nextflow-in-2022.jpg
tags=learn,workshop,webinar
status=published
author=Evan Floden
icon=evan.jpg
~~~~~~

A lot has happened since we last wrote about how best to learn Nextflow, over a year ago. Several new resources have been released including a new Nextflow [Software Carpentries](https://carpentries-incubator.github.io/workflows-nextflow/index.html) course and an excellent write-up by [23andMe](https://www.23andme.com).

We have collated some links below from a diverse collection of resources to help you on your journey to learn Nextflow. Nextflow is a community-driven project - if you have any suggestions, please make a pull request to [this page on GitHub](https://github.com/nextflow-io/website/tree/master/content/blog/2022/learn-nextflow-in-2022.md). 

Without further ado, here is the definitive guide for learning Nextflow in 2022. These resources will support anyone in the journey from total beginner to Nextflow expert.

### Prerequisites

Before you start writing Nextflow pipelines, we recommend that you are comfortable with using the command-line and understand the basic concepts of scripting languages such as Python or Perl. Nextflow is widely used for bioinformatics applications, and scientific data analysis. The examples and guides below often focus on applications in these areas. However, Nextflow is now adopted in a number of data-intensive domains such as image analysis, machine learning, astronomy and geoscience.

### Time commitment

We estimate that it will take at least 20 hours to complete the material. How quickly you finish will depend on your background and how deep you want to dive into the content. Most of the content is introductory but there are some more advanced dataflow and configuration concepts outlined in the workshop and pattern sections.

### Contents

* Why learn Nextflow?
* Introduction to Nextflow from 23andMe
* An RNA-Seq hands-on tutorial 
* Nextflow workshop from Seqera Labs
* Software Carpentries Course
* Managing Pipelines in the Cloud
* The nf-core tutorial
* Advanced implementation patterns
* Awesome Nextflow
* Further resources

### 1. Why learn Nextflow?

Nextflow is an open-source workflow framework for writing and scaling data-intensive computational pipelines. It is designed around the Linux philosophy of simple yet powerful command-line and scripting tools that, when chained together, facilitate complex data manipulations. Combined with support for containerization, support for major cloud providers and on-premise architectures, Nextflow simplifies the writing and deployment of complex data pipelines on any infrastructure. 

The following are some high-level motivations on why people choose to adopt Nextflow:

1. Integrating Nextflow in your analysis workflows helps you implement **reproducible** pipelines. Nextflow pipelines follow FAIR guidelines with version-control and containers to manage all software dependencies.
2. Avoid vendor lock-in by ensuring portability. Nextflow is **portable**; the same pipeline written on a laptop can quickly scale to run on an HPC cluster, Amazon and Google cloud services, and Kubernetes. The code stays constant across varying infrastructures allowing collaboration and avoiding lock-in.
3. It is **scalable** allowing the parallelization of tasks using the dataflow paradigm without having to hard-code the pipeline to a specific platform architecture.
4. It is **flexible** and supports scientific workflow requirements like caching processes to prevent re-computation, and workflow reports to better understand the workflows’ executions.
5. It is **growing fast** and has **long-term support** available from Seqera Labs. Developed since 2013 by the same team, the Nextflow ecosystem is expanding rapidly.
6. It is **open source** and licensed under Apache 2.0. You are free to use it, modify it and distribute it.

### 2. Introduction to Nextflow by 23andMe

This informative post begins with the basic concepts of Nextflow and builds towards how Nextflow is used at 23andMe. It includes a detailed use case for how 23andMe run their imputation pipeline in the cloud, processing over 1 million individuals per day with over 10,000 CPUs in a single compute environment.

👉 [Nextflow at 23andMe](https://medium.com/23andme-engineering/introduction-to-nextflow-4d0e3b6768d1)

### 3. A simple RNA-Seq hands-on tutorial

This hands-on tutorial from Seqera Labs will guide you through implementing a proof-of-concept RNA-seq pipeline. The goal is to become familiar with basic concepts, including how to define parameters, using channels to pass data around and writing processes to perform tasks. It includes all scripts, input data and resources and is perfect for getting a taste of Nextflow.

👉 [Tutorial link on GitHub](https://github.com/seqeralabs/nextflow-tutorial)

### 4. Nextflow workshop from Seqera Labs

Here you’ll dive deeper into Nextflow’s most prominent features and learn how to apply them. The full workshop includes an excellent section on containers, how to build them and how to use them with Nextflow. The written materials come with examples and hands-on exercises. Optionally, you can also follow with a series of videos from a live training workshop.

The workshop includes topics on:

* Environment Setup 
* Basic NF Script and Concepts
* Nextflow Processes
* Nextflow Channels
* Nextflow Operators
* Basic RNA-Seq pipeline
* Containers & Conda
* Nextflow Configuration
* On-premise & Cloud Deployment
* DSL 2 & Modules
* [GATK hands-on exercise](https://seqera.io/training/handson/)

👉 [Workshop](https://seqera.io/training) & [YouTube playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmUv4W8ZRlmstkZwhb_fencI).

### 5. Software Carpentry workshop

The [Nextflow Software Carpentry](https://carpentries-incubator.github.io/workflows-nextflow/index.html) workshop (in active development) motivates the use of Nextflow and [nf-core](https://nf-co.re/) as development tools for building and sharing reproducible data science workflows. The intended audience are those with little programming experience, and the course provides a foundation to comfortably write and run Nextflow and nf-core workflows. Adapted from the Seqera training material above, the workshop has been updated by Software Carpentries instructors within the nf-core community to fit [The Carpentries](https://carpentries.org/) style of training. The Carpentries emphasize feedback to improve teaching materials so we would like to hear back from you about what you thought was both well-explained and what needs improvement. Pull requests to the course material are very welcome.

The workshop can be opened on [Gitpod](https://gitpod.io/#https://github.com/carpentries-incubator/workflows-nextflow) where you can try the exercises in an online computing environment at your own pace, with the course material in another window alongside. 

👉 You can find the course in [The Carpentries incubator](https://carpentries-incubator.github.io/workflows-nextflow/index.html).

### 6. Nextflow Tutorial - Variant Calling Edition

The [Nextflow Tutorial - Variant Calling Edition](https://sateeshperi.github.io/nextflow_varcal/nextflow/) has been adapted from the [Nextflow Software Carpentry](https://carpentries-incubator.github.io/workflows-nextflow/index.html) training material above and [Data Carpentry: Wrangling Genomics Lesson](https://datacarpentry.org/wrangling-genomics/). Learners will have the chance to learn along with nextflow basics, nf-core, to convert a variant-calling bash-script into a nextflow workflow and modularize the pipeline using DSL2 modules and sub-workflows. 

The workshop can be opened on [Gitpod](https://gitpod.io/#https://github.com/sateeshperi/nextflow_tutorial.git) where you can try the exercises in an online computing environment at your own pace, with the course material in another window alongside. 

👉 You can find the course in [Nextflow Tutorial - Variant Calling Edition](https://sateeshperi.github.io/nextflow_varcal/nextflow/).

### 7. Managing Pipelines in the Cloud - GenomeWeb Webinar

This on-demand webinar features Phil Ewels from SciLifeLab and nf-core, Brendan Boufler from Amazon Web Services and Evan Floden from Seqera Labs. The wide ranging dicussion covers the significance of scientific workflow, examples of Nextflow in production settings and how Nextflow can be integrated with other processes.

👉 [Watch the webinar](https://seqera.io/webinars-and-podcasts/managing-bioinformatics-pipelines-in-the-cloud-to-do-more-science/)

### 8. Nextflow implementation patterns

This advanced section discusses recurring patterns and solutions to many common implementation requirements. Code examples are available with notes to follow along, as well as a GitHub repository.

👉 [Nextflow Patterns](http://nextflow-io.github.io/patterns/index.html) & [GitHub repository](https://github.com/nextflow-io/patterns).

### 9. nf-core tutorials

A tutorial covering the basics of using and creating nf-core pipelines. It provides an overview of the nf-core framework including:

* How to run nf-core pipelines
* What are the most commonly used nf-core tools
* How to make new pipelines using the nf-core template
* What are nf-core shared modules
* How to add nf-core shared modules to a pipeline
* How to make new nf-core modules using the nf-core module template
* How nf-core pipelines are reviewed and ultimately released

👉 [nf-core usage tutorials](https://nf-co.re/usage/usage_tutorials)
and [nf-core developer tutorials](https://nf-co.re/developers/developer_tutorials)

### 10. Awesome Nextflow 

A collections of awesome Nextflow pipelines.

👉 [Awesome Nextflow](https://github.com/nextflow-io/awesome-nextflow) on GitHub

### 11. Further resources

The following resources will help you dig deeper into Nextflow and other related projects like the nf-core community who maintain curated pipelines and a very active Slack channel. There are plenty of Nextflow tutorials and videos online, and the following list is no way exhaustive. Please let us know if we are missing anything.

#### Nextflow docs

The reference for the Nextflow language and runtime. These docs should be your first point of reference while developing Nextflow pipelines. The newest features are documented in edge documentation pages released every month with the latest stable releases every three months.

👉 Latest [stable](https://www.nextflow.io/docs/latest/index.html) & [edge](https://www.nextflow.io/docs/edge/index.html) documentation.

#### Seqera Labs docs

An index of documentation, deployment guides, training materials and resources for all things Nextflow and Tower.

👉 [Seqera Labs docs](https://seqera.io/docs)

#### nf-core

nf-core is a growing community of Nextflow users and developers. You can find curated sets of biomedical analysis pipelines written in Nextflow and built by domain experts. Each pipeline is stringently reviewed and has been implemented according to best practice guidelines. Be sure to sign up to the Slack channel.

👉 [nf-core website](https://nf-co.re) and [nf-core Slack](https://nf-co.re/join)

#### Nextflow Tower

Nextflow Tower is a platform to easily monitor, launch and scale Nextflow pipelines on cloud providers and on-premise infrastructure. The documentation provides details on setting up compute environments, monitoring pipelines and launching using either the web graphic interface, CLI or API.

👉 [Nextflow Tower](https://tower.nf) and [user documentation](http://help.tower.nf).

#### Nextflow Biotech Blueprint by AWS

A quickstart for deploying a genomics analysis environment on Amazon Web Services (AWS) cloud, using Nextflow to create and orchestrate analysis workflows and AWS Batch to run the workflow processes.

👉 [Biotech Blueprint by AWS](https://aws.amazon.com/quickstart/biotech-blueprint/nextflow/)

#### Nextflow Data Pipelines on Azure Batch

Nextflow on Azure requires at minimum two Azure services, Azure Batch and Azure Storage. Follow the guides below to set up both services on Azure, and to get your storage and batch account names and keys.

👉 [Azure Blog](https://techcommunity.microsoft.com/t5/azure-compute-blog/running-nextflow-data-pipelines-on-azure-batch/ba-p/2150383) and [GitHub repository](https://github.com/microsoft/Genomics-Quickstart/blob/main/03-Nextflow-Azure/README.md).

#### Running Nextflow by Google Cloud

A step-by-step guide to launching Nextflow Pipelines in Google Cloud.

👉 [Nextflow on Google Cloud](https://cloud.google.com/life-sciences/docs/tutorials/nextflow)

###  Community and support 

* Nextflow [Gitter channel](https://gitter.im/nextflow-io/nextflow)
* Nextflow [Forums](https://groups.google.com/forum/#!forum/nextflow)
* Nextflow Twitter [@nextflowio](https://twitter.com/nextflowio?lang=en)
* [nf-core Slack](https://nfcore.slack.com/)
* [Seqera Labs](https://www.seqera.io) and [Nextflow Tower](https://tower.nf)

### Credits

Special thanks to Mahesh Binzer-Panchal for reviewing the latest revision of this post and contributing the Software Carpentry workshop section.
