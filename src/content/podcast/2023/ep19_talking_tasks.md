---
title: "Nextflow runtime updates: Talking tasks"
episode: 19
description: We discuss task provenance, cloud cache databases and task grouping.
date: 2023-07-31
type: podcast
subtype: News and Views
youtubeid: HlOb_10KV-o
image: /img/podcast_ep19.jpg
tags: news and views,opensource,community
author: Developer advocates
icon: logo_podcast_channels.jpg
---

In this News and Views episode, [Phil Ewels](https://twitter.com/tallphil) and [Chris Hakkaart](https://twitter.com/chris_hakk) talk to Seqera software engineer [Ben Sherman](https://github.com/bentsherman) about what's going on "under the hood" with Nextflow. We touch on task provenance, cloud cache databases and task / sub-workflow grouping. Gotta cache 'em all!

<!-- end-archive-description -->

#### Nextflow runtime

We covered the following topics in this episode:

- Task provenance
  - [Discussion topic](https://github.com/nextflow-io/nextflow/discussions/3447) on the Nextflow GitHub repo
  - [Pull-request #3802](https://github.com/nextflow-io/nextflow/pull/3802) to add a task-graph
  - Related:
    - [Bruno Grande's](https://github.com/BrunoGrandePhD) Nextflow plugin, [nf-prov](https://github.com/Sage-Bionetworks-Workflows/nf-prov)
    - [Quilt Data](https://quiltdata.com/) plugin: [nf-quilt](https://github.com/quiltdata/nf-quilt)
    - Summit 2022 [recap blog post](https://www.nextflow.io/blog/2022/nextflow-summit-2022-recap.html) and [hackathon recap talk](https://summit.nextflow.io/2022/program/oct-14-nf-core-hackathon-report/).
- Cloud cache database
  - [Pull-request #4097](https://github.com/nextflow-io/nextflow/pull/4097) to add support for cloud-based cache directories
- Task and sub-workflow grouping
  - [Issue #2527](https://github.com/nextflow-io/nextflow/issues/2527) on the Nextflow repo
  - [Draft PR #3909](https://github.com/nextflow-io/nextflow/pull/3909) to add Task Grouping.

#### Upcoming events

- [Nextflow Summit 2023](https://summit.nextflow.io/) in Barcelona (October 16-20) and in Boston (November 28-30).
  - Call for abstracts is currently open, just! Get them in quick!
- Two community training events are open for registration:
  - [Foundational training](https://nf-co.re/events/2023/training-basic-2023) - September 6-8, 2023
  - [Advanced training](https://nf-co.re/events/2023/training-sept-2023) - September 27-28, 2023
