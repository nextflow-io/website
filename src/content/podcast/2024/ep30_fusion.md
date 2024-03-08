---
title: Fusion file system
episode: 30
description: Super-fast data access for Nextflow pipelines
date: 2024-02-06
type: podcast
subtype: Technical discussion
youtubeid: mMKnCS2Oy2A
image: /img/podcast_ep30.png
tags: nextflow,opensource,fusion
author: Developer advocates
icon: logo_podcast_channels.jpg
---

In this episode of Channels, [Phil Ewels](https://twitter.com/tallphil) talks to [Paolo Di Tommaso](https://twitter.com/PaoloDiTommaso) (creator of Nextflow, Seqera CTO & cofounder) and [Jordi Deu Pons](https://github.com/jordeu) (software engineer @ Seqera) about Fusion - a file system written specifically for Nextflow.

<!-- end-archive-description -->

We talk about how _"Fusion is not yet another FUSE driver"_ and how it's heavily optimised for Nextflow data pipelines.

Specifically, Fusion is:

- Designed for single job execution, runs in the job container
- Able to do pre-fetching parallel download, with async parallel upload
- Has support for file links over object storage
- Eases data transfer pressure on Nextflow driver app
- Is (almost) zero-config to use
- No need anymore for custom AMI to run Nextflow on AWS Batch

We chat about how it's different to other comparable products, such
as AWS Mountpoint, Goofys, AWS FSx and others and pick over some
benchmark results in detail.
We also clarify two super important points about Fusion:

- The difference between Fusion v1 / v2 (they're totally different tools)
- What AWS NVMe disks are, and why they matter

We touch on some super-powers which are unique to Fusion:
it's multi-cloud and multi-region abilities, abililty to work on HPC
and wrap up by looking to the future to see what's on the horizon for Fusion in 2024.

If you'd like to read more about Fusion, please see the following links:

- [Fusion homepage](https://seqera.io/fusion/)
- [Nextflow docs for Fusion](https://www.nextflow.io/docs/latest/fusion.html)
- [Community forum](https://community.seqera.io/c/fusion/9)


Finally, Phil mentioned some recent and upcoming community content:

- New blog post: [Optimizing Nextflow for HPC and Cloud at Scale](https://nextflow.io/blog/2024/optimizing-nextflow-for-hpc-and-cloud-at-scale.html)
- The upcoming online / distributed [nf-core hackathon](https://nf-co.re/events/2024/hackathon-march-2024), March 18-20, 2024
