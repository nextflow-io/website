title= Learn Nextflow in 2023
date=2023-02-24
type=post
description=An updated compilation of the best Nextflow learning resources for 2023.
image=img/learning_nextflow_2023.png
tags=nextflow, tower
status=published
author=Evan Floden
icon=evan.jpg
~~~~~~

In 2023, the world of Nextflow is more exciting than ever! With new resources constantly being released, there is no better time to dive into this powerful tool. From a new [Software Carpentries’](https://carpentries-incubator.github.io/workflows-nextflow/index.html) course to an [in-depth write-up by 23andMe](https://medium.com/23andme-engineering/introduction-to-nextflow-4d0e3b6768d1) to [new tutorials on Wave and Fusion](https://github.com/seqeralabs/wave-showcase), the options for learning Nextflow are endless.

We've compiled a list of the best resources in 2023 to make your journey to Nextflow mastery as seamless as possible. And remember, Nextflow is a community-driven project. If you have suggestions or want to contribute to this list, head to the [GitHub page](https://github.com/nextflow-io/) and make a pull request.

## Before you start

Before learning Nextflow, you should be comfortable with the Linux command line and be familiar with some basic scripting languages, such as Perl or Python. The beauty of Nextflow is that task logic can be written in your language of choice. You will just need to learn Nextflow’s domain-specific language (DSL) to control overall flow.

Nextflow is widely used in bioinformatics, so many tutorials focus on life sciences. However, Nextflow can be used for almost any data-intensive workflow, including image analysis, ML model training, astronomy, and geoscience applications.

So, let's get started! These resources will guide you from beginner to expert and make you unstoppable in the field of scientific workflows.

## Contents

- [Why Learn Nextflow](#why-learn-nextflow)
- [Meet the Tutorials!](#meet-the-tutorials)
    1. [Basic Nextflow Community Training](#introduction-to-nextflow-by-community)
    2. [Hands-on Nextflow Community Training](#nextflow-hands-on-by-community)
    3. [Advanced Nextflow Community Training](#advanced-nextflow-by-community)
    4. [Software Carpentry workshop](#software-carpentry-workshop)
    5. [An introduction to Nextflow course by Uppsala University](#intro-nexflow-by-uppsala)
    6. [Introduction to Nextflow workshop by VIB](#intro-nextflow-by-vib)
    7. [Nextflow Training by Curtin Institute of Radio Astronomy (CIRA)](#nextflow-training-cira)
    8. [Managing Pipelines in the Cloud - GenomeWeb Webinar](#managing-pipelines-in-the-cloud-genomeweb-webinar)
    9. [Nextflow implementation patterns](#nextflow-implementation-patterns)
    10. [nf-core tutorials](#nf-core-tutorials)
    11. [Awesome Nextflow](#awesome-nextflow)
    12. [Wave showcase: Wave and Fusion tutorials](#wave-showcase-wave-and-fusion-tutorials)
    13. [Building Containers for Scientific Workflows](#building-containers-for-scientific-workflows)
    14. [Best Practices for Deploying Pipelines with Nextflow Tower](#best-practices-for-deploying-pipelines-with-nextflow-tower)
- [Cloud integration tutorials](#cloud-integration-tutorials)
    1. [Nextflow and AWS Batch  Inside the Integration](#nextflow-and-aws-batch-inside-the-integration)
    2. [Nextflow and Azure Batch  Inside the Integration](#nextflow-and-azure-batch-inside-the-integration)
    3. [Get started with Nextflow on Google Cloud Batch](#get-started-with-nextflow-on-google-cloud-batch)
    4. [Nextflow and K8s Rebooted: Running Nextflow on Amazon EKS](#nextflow-and-k8s-rebooted-running-nextflow-on-amazon-eks)
- [Additional resources](#additional-resources)
    1. [Nextflow docs](#nextflow-docs)
    2. [Seqera Labs docs](#seqera-labs-docs)
    3. [nf-core](#nf-core)
    4. [Nextflow Tower](#nextflow-tower)
    5. [Nextflow on AWS](#nextflow-on-aws)
    6. [Nextflow Data pipelines on Azure Batch](#nextflow-data-pipelines-on-azure-batch)
    7. [Running Nextflow with Google Life Sciences](#running-nextflow-with-google-life-sciences)
    8. [Bonus: Nextflow Tutorial - Variant Calling Edition](#bonus-nextflow-tutorial-variant-calling-edition)
- [Community and support](#community-and-support)

<h2 id="why-learn-nextflow">Why Learn Nextflow</h2>

There are hundreds of workflow managers to choose from. In fact, Meir Wahnon and several of his colleagues have gone to the trouble of compiling an awesome-workflow-engines list. The workflows community initiative is another excellent source of information about workflow engines.

- Using Nextflow in your analysis workflows helps you implement reproducible pipelines. Nextflow pipelines follow [FAIR guidelines](https://www.go-fair.org/fair-principles/) (findability, accessibility, interoperability, and reuse). Nextflow also supports version control and containers to manage all software dependencies.
- Nextflow is portable; the same pipeline written on a laptop can quickly scale to run on an HPC cluster, Amazon AWS, Microsoft Azure, Google Cloud Platform, or Kubernetes. With features like [configuration profiles](https://nextflow.io/docs/latest/config.html?#config-profiles), code can be written so that it is 100% portable across different on-prem and cloud infrastructures enabling collaboration and avoiding lock-in.
- It is massively **scalable**, allowing the parallelization of tasks using the dataflow paradigm without hard-coding pipelines to specific platforms, workload managers, or batch services.
- Nextflow is **flexible**, supporting scientific workflow requirements like caching processes to avoid redundant computation and workflow reporting to help understand and diagnose workflow execution patterns.
- It is **growing fast**, and **support is available** from [Seqera Labs](https://seqera.io). The project has been active since 2013 with a vibrant developer community, and the Nextflow ecosystem continues to expand rapidly.
- Finally, Nextflow is open source and licensed under Apache 2.0. You are free to use it, modify it, and distribute it.

<h2 id="meet-the-tutorials">Meet the Tutorials!</h2>

Some of the best publicly available tutorials are listed below:

<h3 id="introduction-to-nextflow-by-community">1. Basic Nextflow Community Training </h3>

Basic training for all things Nextflow. Perfect for anyone looking to get to grips with using Nextflow to run analyses and build workflows. This is the primary Nextflow training material used in most Nextflow and nf-core training events. It covers a large number of topics, with both theoretical and hands-on chapters.

[Basic Nextflow Community Training](https://training.nextflow.io/basic_training/)

We run a free online training event for this course approximately every six months. Videos are streamed to YouTube and questions are handled in the nf-core Slack community. You can watch the recording of the most recent training ([September, 2023](https://nf-co.re/events/2023/training-basic-2023)) in the [YouTube playlist](https://youtu.be/ERbTqLtAkps?si=6xDoDXsb6kGQ_Qa8) below:

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/watch?v=ERbTqLtAkps&list=PL3xpfTVZLcNiLFLiDqk_H5b3TBwvgO_-W" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

<h3 id="nextflow-hands-on-by-community">2. Hands-on Nextflow Community Training </h3>

A "learn by doing" tutorial with less focus on theory, instead leading through exercises of slowly increasing complexity. This course is quite short and hands-on, great if you want to practice your Nextflow skills.

[Hands-on Nextflow Community Training](https://training.nextflow.io/hands_on/)

You can watch the recording of the most recent training ([September, 2023](https://nf-co.re/events/2023/training-hands-on-2023/)) below:

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/x5klpxczAXA?si=moNZUFGd4veMdtC8" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

<h3 id="advanced-nextflow-by-community">3. Advanced Nextflow Community Training </h3>

An advanced material exploring the advanced features of the Nextflow language and runtime, and how to use them to write efficient and scalable data-intensive workflows. This is the Nextflow training material used in advanced training events.

[Advanced Nextflow Community Training](https://training.nextflow.io/advanced/)

You can watch the recording of the most recent training ([September, 2023](https://nf-co.re/events/2023/training-sept-2023/)) below:

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/nPAH9owvKvI?si=-1If8F5DcLqa-cd2" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

<h3 id="software-carpentry-workshop">4. Software Carpentry workshop</h3>

The [Nextflow Software Carpentry](https://carpentries-incubator.github.io/workflows-nextflow/index.html) workshop (still being developed) explains the use of Nextflow and [nf-core](https://nf-co.re/) as development tools for building and sharing reproducible data science workflows. The intended audience is those with little programming experience. The course provides a foundation to write and run Nextflow and nf-core workflows comfortably. Adapted from the Seqera training material above, the workshop has been updated by Software Carpentries instructors within the nf-core community to fit The Carpentries training style. [The Carpentries](https://carpentries.org/) emphasize feedback to improve teaching materials, so we would like to hear back from you about what you thought was well-explained and what needs improvement. Pull requests to the course material are very welcome.
The workshop can be opened on Gitpod where you can try the exercises in an online computing environment at your own pace while referencing the course material in another window alongside the tutorials.

The workshop can be opened on [Gitpod](https://gitpod.io/#https://github.com/carpentries-incubator/workflows-nextflow) where you can try the exercises in an online computing environment at your own pace while referencing the course material in another window alongside the tutorials.

You can find the course in [The Carpentries incubator](https://carpentries-incubator.github.io/workflows-nextflow/index.html).


<h3 id="intro-nexflow-by-uppsala">5. An introduction to Nextflow course from Uppsala University</h3>

This 5-module course by Uppsala University covers the basics of Nextflow, from running Nextflow pipelines, writing your own pipelines and even using containers and conda.

The course can be viewed [here](https://uppsala.instructure.com/courses/51980/pages/nextflow-1-introduction?module_item_id=328997).

<h3 id="intro-nextflow-by-vib">6. Introduction to Nextflow workshop by VIB</h3>

Workshop materials by VIB (mainly) in DSL2 aiming to get familiar with the Nextflow syntax by explaining basic concepts and building a simple RNAseq pipeline. Highlights also reproducibility aspects with adding containers (docker & singularity).

The course can be viewed [here](https://vibbits-nextflow-workshop.readthedocs.io/en/latest/).

<h3 id="nextflow-training-cira">7. Nextflow Training by Curtin Institute of Radio Astronomy (CIRA)</h3>

This training was prepared for physicists and has examples applied to astronomy which may be interesting for Nextflow users coming from this background!

The course can be viewed [here](https://carpentries-incubator.github.io/Pipeline_Training_with_Nextflow/).

<h3 id="managing-pipelines-in-the-cloud-genomeweb-webinar">8. Managing Pipelines in the Cloud - GenomeWeb Webinar</h3>

This on-demand webinar features Phil Ewels from SciLifeLab, nf-core (now also Seqera Labs), Brendan Boufler from Amazon Web Services, and Evan Floden from Seqera Labs. The wide-ranging discussion covers the significance of scientific workflows, examples of Nextflow in production settings, and how Nextflow can be integrated with other processes.

[Watch the webinar](https://seqera.io/events/managing-bioinformatics-pipelines-in-the-cloud-to-do-more-science/)

<h3 id="nextflow-implementation-patterns">9. Nextflow implementation patterns</h3>

This advanced documentation discusses recurring patterns in Nextflow and solutions to many common implementation requirements. Code examples are available with notes to follow along and a GitHub repository.

[Nextflow Patterns](http://nextflow-io.github.io/patterns/index.html) & [GitHub repository](https://github.com/nextflow-io/patterns).

<h3 id="nf-core-tutorials">10. nf-core tutorials</h3>

A set of tutorials covering the basics of using and creating nf-core pipelines developed by the team at [nf-core](https://nf-co.re/). These tutorials provide an overview of the nf-core framework, including:

- How to run nf-core pipelines
- What are the most commonly used nf-core tools
- How to make new pipelines using the nf-core template
- What are nf-core shared modules
- How to add nf-core shared modules to a pipeline
- How to make new nf-core modules using the nf-core module template
- How nf-core pipelines are reviewed and ultimately released

[nf-core usage tutorials](https://nf-co.re/docs/usage/tutorials) and [nf-core developer tutorials](https://nf-co.re/docs/contributing/tutorials).

<h3 id="awesome-nextflow">11. Awesome Nextflow</h3>

A collection of awesome Nextflow pipelines compiled by various contributors to the open-source Nextflow project.

[Awesome Nextflow](https://github.com/nextflow-io/awesome-nextflow) and GitHub

<h3 id="wave-showcase-wave-and-fusion-tutorials">12. Wave showcase: Wave and Fusion tutorials</h3>

Wave and the Fusion file system are new Nextflow capabilities introduced in November 2022. Wave is a container provisioning and augmentation service fully integrated with the Nextflow ecosystem. Instead of viewing containers as separate artifacts that need to be integrated into a pipeline, Wave allows developers to manage containers as part of the pipeline itself.

Tightly coupled with Wave is the new Fusion 2.0 file system. Fusion implements a virtual distributed file system and presents a thin client, allowing data hosted in AWS S3 buckets (and other object stores in the future) to be accessed via the standard POSIX filesystem interface expected by most applications.

Wave can help simplify development, improve reliability, and make pipelines easier to maintain. It can even improve pipeline performance. The optional Fusion 2.0 file system offers further advantages, delivering performance on par with FSx for Lustre while enabling organizations to reduce their cloud computing bill and improve pipeline efficiency throughput. See the [blog article](https://seqera.io/blog/breakthrough-performance-and-cost-efficiency-with-the-new-fusion-file-system/) released in February 2023 explaining the Fusion file system and providing benchmarks comparing Fusion to other data handling approaches in the cloud.

[Wave showcase](https://github.com/seqeralabs/wave-showcase) on GitHub

<h3 id="building-containers-for-scientific-workflows">13. Building Containers for Scientific Workflows</h3>

While not strictly a guide about Nextflow, this article provides an overview of scientific containers and provides a tutorial involved in creating your own container and integrating it into a Nextflow pipeline. It also provides some useful tips on troubleshooting containers and publishing them to registries.

[Building Containers for Scientific Workflows](https://seqera.io/blog/building-containers-for-scientific-workflows/)

<h3 id="best-practices-for-deploying-pipelines-with-nextflow-tower">14. Best Practices for Deploying Pipelines with Nextflow Tower</h3>

When building Nextflow pipelines, a best practice is to supply a nextflow_schema.json file describing pipeline input parameters. The benefit of adding this file to your code repository, is that if the pipeline is launched using Nextflow, the schema enables an easy-to-use web interface that users through the process of parameter selection. While it is possible to craft this file by hand, the nf-core community provides a handy schema build tool. This step-by-step guide explains how to adapt your pipeline for use with Nextflow Tower by using the schema build tool to automatically generate the nextflow_schema.json file.

[Best Practices for Deploying Pipelines with Nextflow Tower](https://seqera.io/blog/best-practices-for-deploying-pipelines-with-nextflow-tower/)

<h2 id="cloud-integration-tutorials">Cloud integration tutorials</h2>

In addition to the learning resources above, several step-by-step integration guides explain how to run Nextflow pipelines on your cloud platform of choice. Some of these tutorials extend to the use of [Nextflow Tower](https://cloud.tower.nf/). Organizations can use the Tower Cloud Free edition to launch pipelines quickly in the cloud. Organizations can optionally use Tower Cloud Professional or run self-hosted or on-premises Tower Enterprise environments as requirements grow. This year, we added Google Cloud Batch to the cloud services supported by Nextflow.

<h3 id="nextflow-and-aws-batch-inside-the-integration">1. Nextflow and AWS Batch — Inside the Integration</h3>

This three-part series of articles provides a step-by-step guide explaining how to use Nextflow with AWS Batch. The [first of three articles](https://seqera.io/blog/nextflow-and-aws-batch-inside-the-integration-part-1-of-3/) covers AWS Batch concepts, the Nextflow execution model, and explains how the integration works under the covers. The [second article](https://seqera.io/blog/nextflow-and-aws-batch-inside-the-integration-part-2-of-3/) in the series provides a step-by-step guide explaining how to set up the AWS batch environment and how to run and troubleshoot open-source Nextflow pipelines. The [third article](https://seqera.io/blog/nextflow-and-aws-batch-using-tower-part-3-of-3/) builds on what you've learned, explaining how to integrate workflows with Nextflow Tower and share the AWS Batch environment with other users by "publishing" your workflows to the cloud.

Nextflow and AWS Batch — Inside the Integration ([part 1 of 3](https://seqera.io/blog/nextflow-and-aws-batch-inside-the-integration-part-1-of-3/), [part 2 of 3](https://seqera.io/blog/nextflow-and-aws-batch-inside-the-integration-part-2-of-3/), [part 3 of 3](https://seqera.io/blog/nextflow-and-aws-batch-using-tower-part-3-of-3/))

<h3 id="nextflow-and-azure-batch-inside-the-integration">2. Nextflow and Azure Batch — Inside the Integration</h3>

Similar to the tutorial above, this set of articles does a deep dive into the Nextflow Azure Batch integration. [Part 1](https://seqera.io/blog/nextflow-and-azure-batch-part-1-of-2/) covers Azure Batch and essential concepts, provides an overview of the integration, and explains how to set up Azure Batch and Storage accounts. It also covers deploying a machine instance in the Azure cloud and configuring it to run Nextflow pipelines against the Azure Batch service.

[Part 2](https://seqera.io/blog/nextflow-and-azure-batch-working-with-tower-part-2-of-2/) builds on what you learned in part 1 and shows how to use Azure Batch from within Nextflow Tower Cloud. It provides a walkthrough of how to make the environment set up in part 1 accessible to users through Tower's intuitive web interface.

Nextflow and Azure Batch — Inside the Integration ([part 1 of 2](https://seqera.io/blog/nextflow-and-azure-batch-part-1-of-2/), [part 2 of 2](https://seqera.io/blog/nextflow-and-azure-batch-working-with-tower-part-2-of-2/))

<h3 id="get-started-with-nextflow-on-google-cloud-batch">3. Get started with Nextflow on Google Cloud Batch</h3>

This excellent article by Marcel Ribeiro-Dantas provides a step-by-step tutorial on using Nextflow with Google’s new Google Cloud Batch service. Google Cloud Batch is expected to replace the Google Life Sciences integration over time. The article explains how to deploy the Google Cloud Batch and Storage environments in GCP using the gcloud CLI. It then goes on to explain how to configure Nextflow to launch pipelines into the newly created Google Cloud Batch environment.

[Get started with Nextflow on Google Cloud Batch](https://nextflow.io/blog/2023/nextflow-with-gbatch.html)

<h3 id="nextflow-and-k8s-rebooted-running-nextflow-on-amazon-eks">4. Nextflow and K8s Rebooted: Running Nextflow on Amazon EKS</h3>

While not commonly used for HPC workloads, Kubernetes has clear momentum. In this educational article, Ben Sherman provides an overview of how the Nextflow / Kubernetes integration has been simplified by avoiding the requirement for Persistent Volumes (PVs) and Persistent Volume Claims (PVCs). This detailed guide provides step-by-step instructions for using Amazon EKS as a compute environment complete with how to configure IAM Roles for Kubernetes Services Accounts (IRSA), now an Amazon EKS best practice.

[Nextflow and K8s Rebooted: Running Nextflow on Amazon EKS](https://seqera.io/blog/deploying-nextflow-on-amazon-eks/)

<h2 id="additional-resources">Additional resources</h2>

The following resources will help you dig deeper into Nextflow and other related projects like the nf-core community which maintain curated pipelines and a very active Slack channel. There are plenty of Nextflow tutorials and videos online, and the following list is no way exhaustive. Please let us know if we are missing anything.

<h3 id="nextflow-docs">1. Nextflow docs</h3>

The reference for the Nextflow language and runtime. These docs should be your first point of reference while developing Nextflow pipelines. The newest features are documented in edge documentation pages released every month, with the latest stable releases every three months.

Latest [stable](https://www.nextflow.io/docs/latest/index.html) & [edge](https://www.nextflow.io/docs/edge/index.html) documentation.

<h3 id="seqera-labs-docs">2. Seqera Labs docs</h3>

An index of documentation, deployment guides, training materials, and resources for all things Nextflow and Tower.

[Seqera Labs docs](https://seqera.io/docs/)

<h3 id="nf-core">3. nf-core</h3>

nf-core is a growing community of Nextflow users and developers. You can find curated sets of biomedical analysis pipelines written in Nextflow and built by domain experts. Each pipeline is stringently reviewed and has been implemented according to best practice guidelines. Be sure to sign up for the Slack channel.

[nf-core website](https://nf-co.re/) and [nf-core Slack](https://nf-co.re/join)

<h3 id="nextflow-tower">4. Nextflow Tower</h3>

Nextflow Tower is a platform to easily monitor, launch, and scale Nextflow pipelines on cloud providers and on-premise infrastructure. The documentation provides details on setting up compute environments, monitoring pipelines, and launching using either the web graphic interface, CLI, or API.

[Nextflow Tower](https://tower.nf/) and [user documentation](http://help.tower.nf/).

<h3 id="nextflow-on-aws">5. Nextflow on AWS</h3>

Part of the Genomics Workflows on AWS, Amazon provides a quickstart for deploying a genomics analysis environment on Amazon Web Services (AWS) cloud, using Nextflow to create and orchestrate analysis workflows and AWS Batch to run the workflow processes. While this article is packed with good information, the procedure outlined in the more recent [Nextflow and AWS Batch – Inside the integration](https://seqera.io/blog/nextflow-and-aws-batch-inside-the-integration-part-1-of-3/) series, may be an easier place to start. Some of the steps that previously needed to be performed manually have been updated in the latest integration.

[Nextflow on AWS Batch](https://docs.opendata.aws/genomics-workflows/orchestration/nextflow/nextflow-overview.html)

<h3 id="nextflow-data-pipelines-on-azure-batch">6. Nextflow Data Pipelines on Azure Batch</h3>

Nextflow on Azure requires at minimum two Azure services, Azure Batch and Azure Storage. Follow the guide below developed by the team at Microsoft to set up both services on Azure, and to get your storage and batch account names and keys.

[Azure Blog](https://techcommunity.microsoft.com/t5/azure-compute-blog/running-nextflow-data-pipelines-on-azure-batch/ba-p/2150383) and [GitHub repository](https://github.com/microsoft/Genomics-Quickstart/blob/main/03-Nextflow-Azure/README.md).

<h3 id="running-nextflow-with-google-life-sciences">7. Running Nextflow with Google Life Sciences</h3>

A step-by-step guide to launching Nextflow Pipelines in Google Cloud. Note that this integration process is specific to Google Life Sciences – an offering that pre-dates Google Cloud Batch. If you want to use the newer integration approach, you can also check out the Nextflow blog article [Get started with Nextflow on Google Cloud Batch](https://nextflow.io/blog/2023/nextflow-with-gbatch.html).

[Nextflow on Google Cloud](https://cloud.google.com/life-sciences/docs/tutorials/nextflow]

<h3 id="bonus-nextflow-tutorial-variant-calling-edition">8. Bonus: Nextflow Tutorial - Variant Calling Edition</h3>

This [Nextflow Tutorial - Variant Calling Edition](https://sateeshperi.github.io/nextflow_varcal/nextflow/) has been adapted from the [Nextflow Software Carpentry training material](https://carpentries-incubator.github.io/workflows-nextflow/index.html) and [Data Carpentry: Wrangling Genomics Lesson](https://datacarpentry.org/wrangling-genomics/). Learners will have the chance to learn Nextflow and nf-core basics, to convert a variant-calling bash script into a Nextflow workflow, and modularize the pipeline using DSL2 modules and sub-workflows.

The workshop can be opened on [Gitpod](https://gitpod.io/#https://github.com/sateeshperi/nextflow_tutorial.git), where you can try the exercises in an online computing environment at your own pace, with the course material in another window alongside.

You can find the course in [Nextflow Tutorial - Variant Calling Edition](https://sateeshperi.github.io/nextflow_varcal/nextflow/).

<h2 id="community-and-support">Community and support</h2>

- [Seqera Community Forum](https://community.seqera.io)
- [Nextflow Slack](https://nextflow.slack.com/)
- Nextflow [Forums](https://groups.google.com/forum/#!forum/nextflow) (active, but the Nextflow Slack gets more traffic)
- Nextflow Twitter [@nextflowio](https://twitter.com/nextflowio?lang=en)
- [nf-core Slack](https://nfcore.slack.com/)
- [Seqera Labs](https://www.seqera.io/) and [Nextflow Tower](https://tower.nf/)
- [Nextflow patterns](https://github.com/nextflow-io/patterns)
- [Nextflow Snippets](https://github.com/mribeirodantas/NextflowSnippets)
