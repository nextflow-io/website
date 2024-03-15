---
title: "Coming soon to Nextflow: From command-line to cleanup"
episode: 15
description: "We talk to Nextflow developer Ben Sherman about some exciting new features coming soon to Nextflow: a new CLI v2 interface; submission of job arrays and intermediate file cleanup during runs."
date: 2023-06-05
type: podcast
subtype: News and Views
youtubeid: ni5A1LlCVYI
image: /img/podcast_ep15.jpg
tags: news and views,opensource,community
author: Developer advocates
icon: logo_podcast_channels.jpg
---

In this News and Views episode, [Phil Ewels](https://twitter.com/tallphil) talks to [Ben Sherman](https://github.com/bentsherman) about upcoming features in core Nextflow development: a new CLI v2 interface; submission of job arrays and intermediate file cleanup during runs.

<!-- end-archive-description -->

#### CLI v2

- The new command-line interface will be invoked with the `nf` command
- It will co-exist with the current `nextflow` CLI, at least for a while
- Expect argument / option parsing more in line with linux norms
  - No longer a single-dash for core flags and double-dash for pipeline parameters
  - Single dash for short flags, double for long-form flags
  - Core and pipeline options separated with `--`

See the pull-request here: [Add CLI v2 `#3600`](https://github.com/nextflow-io/nextflow/pull/3600)

#### Job Arrays

- Soon you'll be able to submit job arrays to HPC and cloud
- Group the submission of any number (10s-1000s) of jobs
- Still independent jobs, but submitted in one go

See the pull-request here: [Implement job array using TaskRun as carrier `#3905`](https://github.com/nextflow-io/nextflow/pull/3905)

#### Intermediate file cleanup

- New ability to track which intermediate files are no longer needed
- Allows files to be deleted from the work directory _during a run_
- Reduces the peak storage required for large pipeline runs

See the pull-request here: [Add support for temporary output paths `#3818`](https://github.com/nextflow-io/nextflow/pull/3818)
