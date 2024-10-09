---
title: "Daisy-chaining workflows: the nf-cascade concept"
date: 2024-10-09
type: post
description: ""
image: /img/blog-2024-10-09--share.png
tags: nextflow,ambassador_post
status: published
author: Mahesh Binzer-Panchal
icon: mahesh_binzer-panchal.png
community_post: true
ambassador_post: true
---

Nextflow workflows are widely used across various fields, performing tasks ranging from data retrieval and cleaning to processing and downstream analyses. Generally, workflows which feed into other workflows are run independently. They are also often not written in a way they can be included directly as subworkflows into another workflow, making it challenging to connect them seamlessly. [nf-cascade](https://github.com/mahesh-panchal/nf-cascade) is a proof-of-concept pipeline that demonstrates how multiple existing workflows can be integrated into a single workflow. Specifically, the proof-of-concept connects the nf-core pipelines: nf-core/fetchngs, nf-core/taxprofiler, nf-core/mag, and nf-core/funcscan. The code can be generalised further to run and connect any existing workflow, without any modifications necessary to the workflows you want to run.

The inspiration for nf-cascade came from a question on the Nextflow Slack workspace. A member of the community wanted to benchmark a set of tools and compare them with their own workflow. They needed a way to run their Nextflow workflow within another Nextflow workflow while still utilising container technologies and cluster submission systems.

I proposed running Nextflow natively using a process with [exec:](https://www.nextflow.io/docs/latest/process.html#native-execution). Processes that use exec: execute Groovy code natively, meaning they run locally in the same environment as the Nextflow runner. This allows the child Nextflow workflow to access the same environment as the parent workflow, enabling the use of the same tools for both (e.g., the same container platform and same job scheduler).

One challenge was resuming a child workflow, as the parent workflow would create a new working directory, causing the child workflow to restart from the beginning and lose all progress. To address this, a folder was created in the work directory to execute the child workflow. Upon successful completion, the results are published to the task work directory, allowing the child workflow to resume from where it left off. This approach could also be implemented as a standard Nextflow process.

Realising that any Nextflow workflow could be run this way and connected together, I decided to implement nf-cascade to demonstrate how workflows could be daisy-chained. This feature has been requested many times in the Nextflow and nf-core communities. I had previously [attempted to chain workflows](https://github.com/mahesh-panchal/test_nfcore_workflow_chain) when nf-core used a different pipeline template, but this ended up being a lot of work, and put a large maintenance burden on the developer. At the time, I was also less familiar with native processes so the idea had not occurred to me.

To get started, copy one of the implementations of the NEXTFLOW_RUN process from the repository. The main branch focuses on nf-core pipelines with specific requirements, while the generalise branch is designed to run any Nextflow workflow. Input files should be treated as in any workflow, using appropriate channel factories or file to create file objects. These can then be passed through channels to the NEXTFLOW_RUN process. The output is the results folder. The map and filter channel operators can be used to extract files from the results folder and transform them into input for the next Nextflow workflow. A guided example can be found [on the wiki](https://github.com/mahesh-panchal/nf-cascade/wiki/Guided-Example-%E2%80%90-Nf%E2%80%90core-style-workflow).

The benefits of this method include the straightforward daisy-chaining of workflows. A child workflow runs as it would from the command line, requiring no modifications for integration. However, there are some downsides: you are limited to what the child workflow allows, meaning you must wait for the entire pipeline to finish even if you only need a portion of it. Configuration can also be tricky, particularly in creating profiles to ensure all workflows run under the same execution conditions.

nf-cascade demonstrates how multiple existing workflows can be connected into a single workflow with minimal effort. While there are some caveats, as mentioned above, it is generally applicable and effective. It provides a suitable interim solution for connecting workflows together as a single workflow.

