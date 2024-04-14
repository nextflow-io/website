---
title: "Nextflow plugins: development walkthrough"
episode: 35
description: A walk through of how to create and publish a plugin from scratch.
date: 2024-04-16
type: podcast
subtype: Tutorial
youtubeid: 2NSRA78C-jg
image: /img/podcast_ep35.png
tags: nextflow,opensource,plugins
author: Developer advocates
icon: logo_podcast_channels.jpg
---

Ever wondered about how Nextflow plugins work?
Join us for this walk through tutorial, as
[Ben Sherman](https://github.com/bentsherman/) guides
[Phil Ewels](https://github.com/ewels/) to create and publish a plugin from scratch.
Right from the first line of code to creating a release and publishing the plugin.

<!-- end-archive-description -->

Nextflow plugins have been getting more attention lately and we've been getting a lot of questions in the community.
There is better documentation and resources planned,
but in the mean time we hope that this can be a useful guide for anyone curious in how they work, and interested in getting started.

:::info
**This podcast episode is quite different to normal.**
Rather than an interview or a discussion, it's more of a tutorial.
It'll work best if you follow the video, though hopefully it's still useful with just sound.

We'd love to hear what you think of it!
After this we will be back to our usual routine in the podcast,
but we'd love to hear what you thought of this taster!
:::

Links to resources used / mentioned in this episode:

- [Nextflow Documentation: Plugins](https://nextflow.io/docs/latest/plugins.html)
- Example plugin used as a base: [nextflow-io/nf-hello](https://github.com/nextflow-io/nf-hello/)
- Phil's plugin created during recording: [ewels/channels-demo](https://github.com/ewels/channels-demo/)
- Provenance reporting plugin: [nextflow-io/nf-prov](https://github.com/nextflow-io/nf-prov)
- Schema validation plugin: [nextflow-io/nf-schema](https://github.com/nextflow-io/nf-schema/)

Time stamps:

<ul class="list-unstyled">
  <li><code>00:00:00</code> Welcome</li>
  <li><code>00:03:27</code> Nextflow docs</li>
  <li><code>00:08:34</code> Starting by forking nf-hello</li>
  <li><code>00:10:32</code> Overview of project files</li>
  <li><code>00:17:13</code> Trying a first compile</li>
  <li><code>00:18:12</code> Different publishing methods</li>
  <li><code>00:19:54</code> nf-boost local publish method</li>
  <li><code>00:23:05</code> Trying the new compile</li>
  <li><code>00:24:00</code> Running locally with the plugin</li>
  <li><code>00:27:32</code> Looking at the nf-hello plugin code</li>
  <li><code>00:39:29</code> Deleting files we don't need</li>
  <li><code>00:42:34</code> Finding event names in the Nextflow source</li>
  <li><code>00:47:46</code> Writing some custom functions</li>
  <li><code>00:49:41</code> Testing our new code</li>
  <li><code>00:51:51</code> Modifying the test pipeline to create files</li>
  <li><code>00:54:14</code> Coding up the JSON output</li>
  <li><code>00:59:20</code> Looking at nf-prov code for BCO files</li>
  <li><code>01:05:30</code> Testing JSON output</li>
  <li><code>01:09:39</code> Automatic work dir cleanup / nf-boost</li>
  <li><code>01:12:17</code> Publishing a release on GitHub</li>
  <li><code>01:19:07</code> Custom plugin repositories</li>
  <li><code>01:26:47</code> Publishing for all Nextflow users</li>
  <li><code>01:29:52</code> Conclusion and end</li>
</ul>

<style>
.dl-horizontal dt {
  max-width: 50px;
}
</style>
