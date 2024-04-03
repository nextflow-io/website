---
title: Leveraging nf-test for enhanced quality control in nf-core
date: 2024-04-03
type: post
description: Reproducibility is an important attribute of all good science. This is specially true in the realm of bioinformatics, where software is hopefully being constantly updated, and pipelines are ideally being maintained. This blog post covers nf-test in the nf-core context.
image: /img/blog-2024-03-28--share.jpg
tags: nextflow,nf-core,nf-test,ambassador_post
status: published
author: Carson Miller
icon: Carson.png
author2: Sateesh Peri
icon2: SateeshPeri.png
---

<div class="alert alert-info">
This post has been written by our valued community members. It’s important to note that the opinions shared here are solely those of the authors.
</div>

# The ever-changing landscape of bioinformatics

Reproducibility is an important attribute of all good science. This is especially true in the realm of bioinformatics, where software is **hopefully** being updated, and pipelines are **ideally** being maintained. Improvements and maintenance are great, but they also bring about an important question: Do bioinformatics tools and pipelines continue to run successfully and produce consistent results despite these changes? Fortunately for us, there is an existing approach to ensure software reproducibility: testing.

<!-- end-archive-description -->

# The Wonderful World of Testing

> "Software testing is the process of evaluating and verifying that a software product does what it is supposed to do,"
> Lukas Forer, co-creator of nf-test.

Software testing has two primary purposes: determining whether an operation continues to run successfully after changes are made, and comparing outputs across runs to see if they are consistent. Testing can alert the developer that an output has changed so that an appropriate fix can be made. Admittedly, there are some instances when altered outputs are intentional (i.e., improving a tool might lead to better, and therefore different, results). However, even in these scenarios, it is important to know what has changed, so that no unintentional changes are introduced during an update.

# Writing effective tests

Although having any test is certainly better than having no tests at all, there are several considerations to keep in mind when adding tests to pipelines and/or tools to maximize their effectiveness. These considerations can be broadly categorized into two groups:

1. Which inputs/functionalities should be tested?
2. What contents should be tested?

## Consideration 1: Testing inputs/functionality

Generally, software will have a default or most common use case. For instance, the nf-core [FastQC](https://nf-co.re/modules/fastqc) module is commonly used to assess the quality of paired-end reads in FastQ format. However, this is not the only way to use the FastQC module. Inputs can also be single-end/interleaved FastQ files, BAM files, or can contain reads from multiple samples. Each input type is analyzed differently by FastQC, and therefore, to increase your test coverage (["the degree to which a test or set of tests exercises a particular program or system"](https://www.geeksforgeeks.org/test-design-coverage-in-software-testing/)), a test should be written for each possible input. Additionally, different settings can change how a process is executed. For example, in the [bowtie2/align](https://nf-co.re/modules/bowtie2_align) module, aside from input files, the `save_unaligned` and `sort_bam` parameters can alter how this module functions and the outputs it generates. Thus, tests should be written for each possible scenario. When writing tests, aim to consider as many variations as possible. If some are missed, don't worry! Additional tests can be added later. Discovering these different use cases and how to address/test them is part of the development process.

## Consideration 2: Testing outputs

Once test cases are established, the next step is determining what specifically should be evaluated in each test. Generally, these evaluations are referred to as assertions. Assertions can range from verifying whether a job has been completed successfully to comparing the output channel/file contents between runs. Ideally, tests should incorporate all outputs, although there are scenarios where this is not feasible (for example, outputs containing timestamps or paths). In such cases, it's often best to include at least a portion of the contents from the problematic file or, at the minimum, the name of the file to ensure that it is consistently produced.

# Testing in nf-core

nf-core is a community-driven initiative that aims to provide high-quality, Nextflow-based bioinformatics pipelines. The community's emphasis on reproducibility makes testing an essential aspect of the nf-core ecosystem. Until recently, tests were implemented using pytest for modules/subworkflows and test profiles for pipelines. These tests ensured that nf-core components could run successfully following updates. However, at the pipeline level, they did not check file contents to evaluate output consistency. Additionally, using two different testing approaches lacked the standardization nf-core strives for. An ideal test framework would integrate tests at all Nextflow development levels (functions, modules, subworkflows, and pipelines) and comprehensively test outputs.

# New and Improved Nextflow Testing with nf-test

Created by [Lukas Forer](https://github.com/lukfor) and [Sebastian Schönherr](https://github.com/seppinho), nf-test has emerged as the leading solution for testing Nextflow pipelines. Their goal was to enhance the evaluation of reproducibility in complex Nextflow pipelines. To this end, they have implemented several notable features, creating a robust testing platform:

1. **Comprehensive Output Testing**: nf-test employs [snapshots](https://www.nf-test.com/docs/assertions/snapshots/) for handling complex data structures. This feature evaluates the contents of any specified output channel/file, enabling comprehensive and reliable tests that ensure data integrity following changes.
2. **A Consistent Testing Framework for All Nextflow Components**: nf-test provides a unified framework for testing everything from individual functions to entire pipelines, ensuring consistency across all components.
3. **A DSL for Tests**: Designed in the likeness of Nextflow, nf-test's intuitive domain-specific language (DSL) uses 'when' and 'then' blocks to describe expected behaviors in pipelines, facilitating easier test script writing.
4. **Readable Assertions**: nf-test offers a wide range of functions for writing clear and understandable [assertions](https://www.nf-test.com/docs/assertions/assertions/), improving the clarity and maintainability of tests.
5. **Boilerplate Code Generation**: To accelerate the testing process, nf-test and nf-core tools feature commands that generate boilerplate code, streamlining the development of new tests.

# But wait… there's more!

The merits of having a consistent and comprehensive testing platform are significantly amplified with nf-test's integration into nf-core. This integration provides an abundance of resources for incorporating nf-test into your Nextflow development. Thanks to this collaboration, you can utilize common nf-test commands via nf-core tools and easily install nf-core modules/subworkflows that already have nf-test implemented. Moreover, an [expanding collection of examples](https://nf-co.re/docs/contributing/tutorials/nf-test_assertions) is available to guide you through adopting nf-test for your projects.

# Adding nf-test to pipelines

Several nf-core pipelines have begun to adopt nf-test as their testing framework. Among these, [nf-core/methylseq](https://nf-co.re/methylseq/) was the first to implement pipeline-level nf-tests as a proof-of-concept. However, since this initial implementation, nf-core maintainers have identified that the existing nf-core pipeline template needs modifications to better support nf-test. These adjustments aim to enhance compatibility with nf-test across components (modules, subworkflows, workflows) and ensure that tests are included and shipped with each component. A more detailed blog post about these changes will be published in the future.
Following these insights, [nf-core/fetchngs](https://nf-co.re/fetchngs) has been at the forefront of incorporating nf-test for testing modules, subworkflows, and at the pipeline level. Currently, fetchngs serves as the best-practice example for nf-test implementation within the nf-core community. Other nf-core pipelines actively integrating nf-test include [mag](https://nf-co.re/mag), [sarek](https://nf-co.re/sarek), [readsimulator](https://nf-co.re/readsimulator), and [rnaseq](https://nf-co.re/rnaseq).

# Pipeline development with nf-test

**For newer nf-core pipelines, integrating nf-test as early as possible in the development process is highly recommended**. An example of a pipeline that has benefitted from the incorporation of nf-tests throughout its development is [phageannotator](https://github.com/nf-core/phageannotator). Although integrating nf-test during pipeline development has presented challenges, it has offered a unique opportunity to evaluate different testing methodologies and has been instrumental in identifying numerous development errors that might have been overlooked using the previous test profiles approach. Additionally, investing time early on has significantly simplified modifying different aspects of the pipeline, ensuring that functionality and output remain unaffected.
For those embarking on creating new Nextflow pipelines, here are a few key takeaways from our experience:

1. **Leverage nf-core modules/subworkflows extensively**. Devoting time early to contribute modules/subworkflows to nf-core not only streamlines future development for you and your PR reviewers but also simplifies maintaining, linting, and updating pipeline components through nf-core tools. Furthermore, these modules will likely benefit others in the community with similar research interests.
2. **Prioritize incremental changes over large overhauls**. Incremental changes are almost always preferable to large, unwieldy modifications. This approach is particularly beneficial when monitoring and updating nf-tests at the module, subworkflow, and pipeline levels. Introducing too many changes simultaneously can overwhelm both developers and reviewers, making it challenging to track what has been modified and what requires testing. Aim to keep changes straightforward and manageable.
3. **Facilitate parallel execution of nf-test to generate and test snapshots**. By default, nf-test runs each test sequentially, which can make the process of running multiple tests to generate or updating snapshots time-consuming. Implementing scripts that allow tests to run in parallel—whether via a workload manager or in the cloud—can significantly save time and simplify the process of monitoring tests for pass or fail outcomes.

# Community and contribution

nf-core is a community that relies on consistent contributions, evaluation, and feedback from its members to improve and stay up-to-date. This holds true as we transition to a new testing framework as well. Currently, there are two primary ways that people have been contributing in this transition:

1. **Adding nf-tests to new and existing nf-core modules/subworkflows**. There has been a recent emphasis on migrating modules/subworkflows from pytest to nf-test because of the advantages mentioned previously. Fortunately, the nf-core team has added very helpful [instructions](https://nf-co.re/docs/contributing/modules#migrating-from-pytest-to-nf-test) to the website, which has made this process much more streamlined.
2. **Adding nf-tests to nf-core pipelines**. Another area of focus is the addition of nf-tests to nf-core pipelines. This process can be quite difficult for large, complex pipelines, but there are now several examples of pipelines with nf-tests that can be used as a blueprint for getting started ([fetchngs](https://github.com/nf-core/fetchngs/tree/master), [sarek](https://github.com/nf-core/sarek/tree/master), [rnaseq](https://github.com/nf-core/rnaseq/tree/master), [readsimulator](https://github.com/nf-core/readsimulator/tree/master), [phageannotator](https://github.com/nf-core/phageannotator)).

> These are great areas to work on & contribute in nf-core hackathons

The nf-core community added a significant number of nf-tests during the recent [hackathon in March 2024](https://nf-co.re/events/2024/hackathon-march-2024). Yet the role of the community is not limited to adding test code. A robust testing infrastructure requires nf-core users to identify testing errors, additional test cases, and provide feedback so that the system can continually be improved. Each of us brings a different perspective, and the development-feedback loop that results from collaboration brings about a much more effective, transparent, and inclusive system than if we worked in isolation.

# Future directions

Looking ahead, nf-core and nf-test are poised for tighter integration and significant advancements. Anticipated developments include enhanced testing capabilities, more intuitive interfaces for writing and managing tests, and deeper integration with cloud-based resources. These improvements will further solidify the position of nf-core and nf-test at the forefront of bioinformatics workflow management.

# Conclusion

The integration of nf-test within the nf-core ecosystem marks a significant leap forward in ensuring the reproducibility and reliability of bioinformatics pipelines. By adopting nf-test, developers and researchers alike can contribute to a culture of excellence and collaboration, driving forward the quality and accuracy of bioinformatics research.

Special thanks to everyone in the #nf-test channel in the nf-core slack workspace for their invaluable contributions, feedback, and support throughout this adoption. We are immensely grateful for your commitment and look forward to continuing our productive collaboration.

<br>
<div class="footer-wrapper">
<img src="/img/nextflow_ambassador_logo.svg" height="70" class="pull-right" style="padding: 1rem 1rem;">

> This post was contributed by a Nextflow Ambassador, passionate individuals who support the Nextflow community. Interested in becoming an ambassador? Read more about it <a href="https://www.nextflow.io/ambassadors.html">here</a>.
</div>
