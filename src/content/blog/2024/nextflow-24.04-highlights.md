---
title: "Nextflow 24.04 - Release highlights"
date: 2024-05-27
type: post
description: >
  A highlight of some of the goodies that have just been
  released in the 24.04 stable release of Nextflow 24.04 stable.
image: /img/blog-nextflow-24.04-highlights.png
tags: nextflow
status: published
author: Paolo Di Tommaso
icon: paolo.jpg
author2: Ben Sherman
icon2: ben.jpg
author3: Phil Ewels
icon3: phil.jpg
---

We release an "edge" version of Nextflow every month and a "stable" version every six months. The stable releases are recommended for production usage and represent a significant milestone. The [release changelogs](https://github.com/nextflow-io/nextflow/releases) contain a lot of detail, so we thought we'd highlight some of the goodies that have just been released in Nextflow 24.04 stable. Let's get into it!

:::tip
We also did a podcast episode about some of these changes!
Check it out here: [Channels Episode 41](/podcast/2024/ep41_nextflow_2404.html).
:::

## Table of contents

- [New features](#new-features)
  - [Seqera Containers](#seqera-containers)
  - [Workflow output definition](#workflow-output-definition)
  - [Topic channels](#topic-channels)
  - [Process eval outputs](#process-eval-outputs)
  - [Resource limits](#resource-limits)
  - [Job arrays](#job-arrays)
- [Enhancements](#enhancements)
  - [Colored logs](#colored-logs)
  - [AWS Fargate support](#aws-fargate-support)
  - [OCI auto pull mode for Singularity and Apptainer](#oci-auto-pull-mode-for-singularity-and-apptainer)
  - [Support for GA4GH TES](#support-for-ga4gh-tes)
- [Fusion](#fusion)
  - [Enhanced Garbage Collection](#enhanced-garbage-collection)
  - [Increased File Handling Capacity](#increased-file-handling-capacity)
  - [Correct Publishing of Symbolic Links](#correct-publishing-of-symbolic-links)
- [Other notable changes](#other-notable-changes)

## New features

### Seqera Containers

A new flagship community offering was revealed at the Nextflow Summit 2024 Boston - **Seqera Containers**. This is a free-to-use container cache powered by [Wave](https://seqera.io/wave/), allowing anyone to request an image with a combination of packages from Conda and PyPI. The image will be built on demand and cached (for at least 5 years after creation). There is a [dedicated blog post](https://seqera.io/blog/introducing-seqera-pipelines-containers/) about this, but it's worth noting that the service can be used directly from Nextflow and not only through [https://seqera.io/containers/](https://seqera.io/containers/)

In order to use Seqera Containers in Nextflow, simply set `wave.freeze` _without_ setting `wave.build.repository` - for example, by using the following config for your pipeline:

```groovy
wave.enabled = true
wave.freeze = true
wave.strategy = 'conda'
```

Any processes in your pipeline specifying Conda packages will have Docker or Singularity images created on the fly (depending on whether `singularity.enabled` is set or not) and cached for immediate access in subsequent runs. These images will be publicly available. You can view all container image names with the `nextflow inspect` command.

### Workflow output definition

The workflow output definition is a new syntax for defining workflow outputs:

```groovy
nextflow.preview.output = true // [!code ++]

workflow {
  main:
  ch_foo = foo(data)
  bar(ch_foo)

  publish:
  ch_foo >> 'foo' // [!code ++]
}

output { // [!code ++]
  directory 'results' // [!code ++]
  mode 'copy' // [!code ++]
} // [!code ++]
```

It essentially provides a DSL2-style approach for publishing, and will replace `publishDir` once it is finalized. It also provides extra flexibility as it allows you to publish _any_ channel, not just process outputs. See the [Nextflow docs](https://nextflow.io/docs/latest/workflow.html#publishing-outputs) for more information.

:::info
This feature is still in preview and may change in a future release.
We hope to finalize it in version 24.10, so don't hesitate to share any feedback with us!
:::

### Topic channels

Topic channels are a new channel type introduced in 23.11.0-edge. A topic channel is essentially a queue channel that can receive values from multiple sources, using a matching name or "topic":

```groovy
process foo {
  output:
  val('foo'), topic: 'my-topic' // [!code ++]
}

process bar {
  output:
  val('bar'), topic: 'my-topic' // [!code ++]
}

workflow {
  foo()
  bar()

  Channel.topic('my-topic').view() // [!code ++]
}
```

Topic channels are particularly useful for collecting metadata from various places in the pipeline, without needing to write all of the channel logic that is normally required (e.g. using the `mix` operator). See the [Nextflow docs](https://nextflow.io/docs/latest/channel.html#topic) for more information.

### Process `eval` outputs

Process `eval` outputs are a new type of process output which allows you to capture the standard output of an arbitrary shell command:

```groovy
process sayHello {
  output:
  eval('bash --version') // [!code ++]

  """
  echo Hello world!
  """
}

workflow {
  sayHello | view
}
```

The shell command is executed alongside the task script. Until now, you would typically execute these supplementary commands in the main process script, save the output to a file or environment variable, and then capture it using a `path` or `env` output. The new `eval` output is a much more convenient way to capture this kind of command output directly. See the [Nextflow docs](https://nextflow.io/docs/latest/process.html#output-type-eval) for more information.

#### Collecting software versions

Together, topic channels and eval outputs can be used to simplify the collection of software tool versions. For example, for FastQC:

```groovy
process FASTQC {
  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path('*.html'), emit: html
  tuple val("${task.process}"), val('fastqc'), eval('fastqc --version'), topic: versions // [!code ++]

  """
  fastqc $reads
  """
}

workflow {
  Channel.topic('versions') // [!code ++]
    | unique()
    | map { process, name, version ->
      """\
      ${process.tokenize(':').last()}:
        ${name}: ${version}
      """.stripIndent()
    }
    | collectFile(name: 'collated_versions.yml')
    | CUSTOM_DUMPSOFTWAREVERSIONS
}
```

This approach will be implemented across all nf-core pipelines, and will cut down on a lot of boilerplate code. Check out the full prototypes for nf-core/rnaseq [here](https://github.com/nf-core/rnaseq/pull/1109) and [here](https://github.com/nf-core/rnaseq/pull/1115) to see them in action!

### Resource limits

The **resourceLimits** directive is a new process directive which allows you to define global limits on the resources requested by individual tasks. For example, if you know that the largest node in your compute environment has 24 CPUs, 768 GB or memory, and a maximum walltime of 72 hours, you might specify the following:

```groovy
process.resourceLimits = [ cpus: 24, memory: 768.GB, time: 72.h ]
```

If a task requests more than the specified limit (e.g. due to [retry with dynamic resources](https://nextflow.io/docs/latest/process.html#dynamic-computing-resources)), Nextflow will automatically reduce the task resources to satisfy the limit, whereas normally the task would be rejected by the scheduler or would simply wait in the queue forever! The nf-core community has maintained a custom workaround for this problem, the `check_max()` function, which can now be replaced with `resourceLimits`. See the [Nextflow docs](https://nextflow.io/docs/latest/process.html#resourcelimits) for more information.

### Job arrays

**Job arrays** are now supported in Nextflow using the `array` directive. Most HPC schedulers, and even some cloud batch services including AWS Batch and Google Batch, support a "job array" which allows you to submit many independent jobs with a single job script. While the individual jobs are still executed separately as normal, submitting jobs as arrays where possible puts considerably less stress on the scheduler.

With Nextflow, using job arrays is a one-liner:

```groovy
process.array = 100
```

You can also enable job arrays for individual processes like any other directive. See the [Nextflow docs](https://nextflow.io/docs/latest/process.html#array) for more information.

:::tip
On Google Batch, using job arrays also allows you to pack multiple tasks onto the same VM by using the `machineType` directive in conjunction with the `cpus` and `memory` directives.
:::

## Enhancements

### Colored logs

<div class="row"><div class="col-lg-6">

**Colored logs** have come to Nextflow! Specifically, the process log which is continuously printed to the terminal while the pipeline is running. Not only is it more colorful, but it also makes better use of the available space to show you what's most important. But we already wrote an entire [blog post](https://nextflow.io/blog/2024/nextflow-colored-logs.html) about it, so go check that out for more details!

</div><div class="col-lg-6">

![New coloured output from Nextflow](/img/blog-nextflow-colored-logs/nextflow_coloured_logs.png)

</div></div>

### AWS Fargate support

Nextflow now supports **AWS Fargate** for AWS Batch jobs. See the [Nextflow docs](https://nextflow.io/docs/latest/aws.html#aws-fargate) for details.

### OCI auto pull mode for Singularity and Apptainer

Nextflow now supports OCI auto pull mode both Singularity and Apptainer. Historically, Singularity could run a Docker container image converting to the Singularity image file format via the Singularity pull command and using the resulting image file in the exec command. This adds extra overhead to the head node running Nextflow for converting all container images to the Singularity format.

Now Nextflow allows specifying the option `ociAutoPull` both for Singularity and Apptainer. When enabling this setting Nextflow delegates the pull and conversion of the Docker image directly to the `exec` command.

```groovy
singularity.ociAutoPull = true
```

This results in the running of the pull and caching of the Singularity images to the compute jobs instead of the head job and removing the need to maintain a separate image files cache.

See the [Nextflow docs](https://nextflow.io/docs/latest/config.html#scope-singularity) for more information.

### Support for GA4GH TES

The [Task Execution Service (TES)](https://ga4gh.github.io/task-execution-schemas/docs/) is an API specification, developed by [GA4GH](https://www.ga4gh.org/), which attempts to provide a standard way for workflow managers like Nextflow to interface with execution backends. Two noteworthy TES implementations are [Funnel](https://github.com/ohsu-comp-bio/funnel) and [TES Azure](https://github.com/microsoft/ga4gh-tes).

Nextflow has long supported TES as an executor, but only in a limited sense, as TES did not support some important capabilities in Nextflow such as glob and directory outputs and the `bin` directory. However, with TES 1.1 and its adoption into Nextflow, these gaps have been closed. You can use the TES executor with the following configuration:

```groovy
plugins {
  id 'nf-ga4gh'
}

process.executor = 'tes'
tes.endpoint = '...'
```

See the [Nextflow docs](https://nextflow.io/docs/latest/executor.html#ga4gh-tes) for more information.

:::note
To better facilitate community contributions, the nf-ga4gh plugin will soon be moved from the Nextflow repository into its own repository, `nextflow-io/nf-ga4gh`. To ensure a smooth transition with your pipelines, make sure to explicitly include the plugin in your configuration as shown above.
:::

## Fusion

[Fusion](https://seqera.io/fusion/) is a distributed virtual file system for cloud-native data pipeline and optimized for Nextflow workloads. Nextflow 24.04 now works with a new release, Fusion 2.3. This brings a few notable quality-of-life improvements:

### Enhanced Garbage Collection

Fusion 2.3 features an improved garbage collection system, enabling it to operate effectively with reduced scratch storage. This enhancement ensures that your pipelines run more efficiently, even with limited temporary storage.

### Increased File Handling Capacity

Support for more concurrently open files is another significant improvement in Fusion 2.3. This means that larger directories, such as those used by Alphafold2, can now be utilized without issues, facilitating the handling of extensive datasets.

### Correct Publishing of Symbolic Links

In previous versions, output files that were symbolic links were not published correctly — instead of the actual file, a text file containing the file path was published. Fusion 2.3 addresses this issue, ensuring that symbolic links are published correctly.

These enhancements in Fusion 2.3 contribute to a more robust and efficient filesystem for Nextflow users.

## Other notable changes

- Add native retry on spot termination for Google Batch ([`ea1c1b`](https://github.com/nextflow-io/nextflow/commit/ea1c1b70da7a9b8c90de445b8aee1ee7a7148c9b))
- Add support for instance templates in Google Batch ([`df7ed2`](https://github.com/nextflow-io/nextflow/commit/df7ed294520ad2bfc9ad091114ae347c1e26ae96))
- Allow secrets to be used with `includeConfig` ([`00c9f2`](https://github.com/nextflow-io/nextflow/commit/00c9f226b201c964f67d520d0404342bc33cf61d))
- Allow secrets to be used in the pipeline script ([`df866a`](https://github.com/nextflow-io/nextflow/commit/df866a243256d5018e23b6c3237fb06d1c5a4b27))
- Add retry strategy for publishing ([`c9c703`](https://github.com/nextflow-io/nextflow/commit/c9c7032c2e34132cf721ffabfea09d893adf3761))
- Add `k8s.cpuLimits` config option ([`3c6e96`](https://github.com/nextflow-io/nextflow/commit/3c6e96d07c9a4fa947cf788a927699314d5e5ec7))
- Removed `seqera` and `defaults` from the standard channels used by the nf-wave plugin. ([`ec5ebd`](https://github.com/nextflow-io/nextflow/commit/ec5ebd0bc96e986415e7bac195928b90062ed062))

You can view the full [Nextflow release notes on GitHub](https://github.com/nextflow-io/nextflow/releases/tag/v24.04.0).
