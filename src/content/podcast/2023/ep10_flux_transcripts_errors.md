---
title: The Flux executor, bytesize transcripts and improved error messaging
episode: 10
description: The Flux executor, bytesize transcripts and improved error messaging.
date: 2023-01-18
type: podcast
subtype: News and Views
youtubeid: 5UDIGQVo9hA
image: img/podcast_ep10.jpg
tags: news and views,opensource,community
author:  Developer advocates
icon: logo_podcast_channels.jpg
---

In this News and Views episode, [Chris Hakkaart](https://twitter.com/chris_hakk) and [Marcel Ribeiro-Dantas](https://twitter.com/mribeirodantas) discuss the hottest topics in the Nextflow world.

<!-- end-archive-description -->

#### The Flux executor

- [Flux](https://flux-framework.org/) is a job scheduler akin to Slurm, Sun Grid Engine, or similar but with some extra features and has been added to Nextflow as an executor.
- Flux was originally designed as an HPC job scheduler but is being extended it to more cloud environments and will be one to watch in the future.
- A great example of a feature that was added to Nextflow by the community [nextflow-io/nextflow#3412](https://github.com/nextflow-io/nextflow/pull/3412).

#### Transcripts for bytesize seminars

- Video transcripts create a better user experience by increasing accessibility for users who are deaf or hard of hearing, improving audience comprehension, and enabling flexible viewing in sound-sensitive environments.
- [Franziska Bonath](https://github.com/FranBonath) has been working towards making the [nf-core bytesize](https://nf-co.re/events) seminars more accessible by adding transcripts that can be used for subtitles to all videos.
- All transcripts for all bytesize seminars have now been made and are now under review.
- A massive thank you to Fran, to all those who have presented, and also those who have spent time reviewing the transcripts.

#### Improved error messaging

- Understanding error messages is a common challenge faced by the community.
- There has been a lot of work on improving Nextflow error messages.
    - Examples include improved messaging for errors created by passing a channel to a channel and operators that also work on folders and gave error messages that would refer to a file.
- An example of community feedback guiding development.

#### Upcoming events

- There's a [nf-core/bytesize talk on January 25th](https://nf-co.re/events/2023/bytesize_funcscan) by Jasmin Frangenberg about the nf-core/funcscan pipeline.
    - A bioinformatics best-practice analysis pipeline for the screening of functional components of nucleotide sequences.
- The [Nextflow / nf-core training](https://nf-co.re/events/2023/training-march-2023) is being held March 13-16.
- The [nf-core hackathon](https://nf-co.re/events/2023/hackathon-march-2023) is being held March 27-29.
