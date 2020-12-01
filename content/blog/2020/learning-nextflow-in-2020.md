title=Learning Nextflow in 2020
date=2020-12-01
type=post
tags=nextflow,learning,workshop
status=published
author=Evan Floden & Alain Coletta
icon=evan.jpg
~~~~~~

With the year nearly over, we thought it was about time to pull together the best-of-the-best guide for learning Nextflow in 2020. These resources will support anyone in the journey from total noob to Nextflow expert so this holiday season, give yourself or someone you know the gift of learning Nextflow!


### Prerequisites to get started

We recommend that learners are comfortable with using the command line and the basic concepts of a scripting language such as Python or Perl before they start writing pipelines. Nextflow is widely used for bioinformatics applications, and the examples in these guides often focus on applications in these topics. However, Nextflow is now adopted in a number of data-intensive domains such as radio astronomy, satellite imaging and machine learning. No domain expertise is expected.


### Time commitment

We estimate that the speediest of learners can complete the material in around 12 hours. It all depends on your background and how deep you want to dive into the rabbit-hole! Most of the content is introductory with some more advanced dataflow and configuration material in the workshops and patterns sections.


### Overview of the material

* Why learn Nextflow?
* Introduction to Nextflow - AWS HPC Conference 2020 (8m)
* A simple RNA-Seq hands-on tutorial (2h)
* Full-immersion workshop (8h)
* Nextflow advanced implementation Patterns (2h)
* Other resources
* Community and Support


### 1. Why learn Nextflow?
Nextflow is an open-source workflow framework for writing and scaling data-intensive computational pipelines. It is designed around the Linux philosophy of simple yet powerful command-line and scripting tools that, when chained together, facilitate complex data manipulations. Combined with support for containerization, support for major cloud providers and on-premise architectures, Nextflow simplifies the writing and deployment of complex data pipelines on any infrastructure. 

The following are some high-level motivations on why people choose to adopt Nextflow:

1. Integrating Nextflow in your analysis workflows helps you implement **reproducible** pipelines. Nextflow pipelines follow FDA repeatability and reproducibility guidelines with version-control and containers to manage all software dependencies.
2. Avoid vendor lock-in by ensuring portability. Nextflow is **portable**; the same pipeline written on a laptop can quickly scale to run HPC cluster, Amazon and Google cloud services, and Kubernetes. The code stays constant across varying infrastructures allowing collaboration and avoiding lock-in.
3. It is **scalable** allowing the parallelization of tasks using the dataflow paradigm without having to hard-code the pipeline to a specific platform architecture.
4. It is **flexible** and supports scientific workflow requirements like caching processes to prevent re-computation, and workflow reports to better understand the workflows’ executions.
5. It is **growing fast** and has **long-term support**. Developed since 2013 by the same team, the Nextflow ecosystem is expanding rapidly.
6. It is **open source** and licensed under Apache 2.0. You are free to use it, modify it and distribute it.


### 2. Introduction to Nextflow from the HPC on AWS Conference 2020

This short YouTube video provides a general overview of Nextflow, the motivations behind its development and a demonstration of some of the latest features.

<iframe width="560" height="315" src="https://www.youtube.com/embed/SYhDkUgcOXo" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>


### 3. A simple RNA-Seq hands-on tutorial

This hands-on tutorial will guide you through implementing a proof-of-concept RNA-seq pipeline. The goal is to become familiar with basic concepts, including how to define parameters, use channels for data and write processes to perform tasks. It includes all scripts, data and resources and is perfect for getting a flavor for Nextflow.

[Project repository on GitHub](https://github.com/nextflow-io/rnaseq-nf)


### 4. Full-immersion workshop
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

[Workshop](https://seqera.io/training) & [YouTube playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmUv4W8ZRlmstkZwhb_fencI).

### 5. Nextflow implementation Patterns
This advanced section discusses recurring patterns and solutions to many common implementation requirements. Code examples are available with notes to follow along with as well as a GitHub repository.

[Nextflow Patterns](http://nextflow-io.github.io/patterns/index.html) & [GitHub repository](https://github.com/nextflow-io/patterns).


### Other resources
The following resources will help you dig deeper into Nextflow and other related projects like the nf-core community who maintain curated pipelines and a very active Slack channel. There are plenty of Nextflow tutorials and videos online, and the following list is no way exhaustive. Please let us know if we are missing something.

#### Nextflow docs
The reference for the Nextflow language and runtime. The docs should be your first point of reference when something is not clear. Newest features are documented in edge documentation pages released every month with the latest stable releases every three months.

Latest [stable](https://www.nextflow.io/docs/latest/index.html) & [edge](https://www.nextflow.io/docs/edge/index.html) documentation.

#### nf-core
nf-core is a growing community of Nextflow users and developers. You can find curated sets of biomedical analysis pipelines built by domain experts with Nextflow, that have passed tests and have been implemented according to best practice guidelines. Be sure to sign up to the Slack channel.

[nf-core website](https://nf-co.re)


#### Tower Docs
Nextflow Tower is a platform to easily monitor, launch and scale Nextflow pipelines on cloud providers and on-premise infrastructure. The documentation provides details on setting up compute environments, monitoring pipelines and launching using either the web graphic interface or API.

[Nextflow Tower documentation](http://help.tower.nf)


#### Nextflow Biotech Blueprint by AWS

A quickstart for deploying a genomics analysis environment on Amazon Web Services (AWS) cloud, using Nextflow to create and orchestrate analysis workflows and AWS Batch to run the workflow processes.

[Biotech Blueprint by AWS](https://aws.amazon.com/quickstart/biotech-blueprint/nextflow/)


#### Running Nextflow by Google Cloud
Google Cloud Nextflow step-by-step guide to launching Nextflow Pipelines in Google Cloud.

[Nextflow on Google Cloud ](https://cloud.google.com/life-sciences/docs/tutorials/nextflow)

#### Awesome Nextflow 

A collections of Nextflow based pipelines and other resources. 

[Awesome Nextflow](https://github.com/nextflow-io/awesome-nextflow)


###  Community and support 

* Nextflow [Gitter channel](https://gitter.im/nextflow-io/nextflow)
* Nextflow [Forums](https://groups.google.com/forum/#!forum/nextflow)
* [nf-core Slack](https://nfcore.slack.com/)
* Twitter [@nextflowio](https://twitter.com/nextflowio?lang=en)
* [Seqera Labs](https://www.seqera.io) technical support & consulting


Nextflow is a community-driven project. The list of links below has been collated from a diverse collection of resources and experts to guide you in learning Nextflow. If you have any suggestions, please make a pull request to this page on GitHub.

Also stay tuned for our upcoming post, where we will discuss the ultimate Nextflow development environment. 

