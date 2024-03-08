---
title: Nextflow meets Dockstore
date: 2018-09-18
type: post
tags: nextflow,ga4gh,nf-core,dockstore
author: Paolo Di Tommaso
icon: paolo.jpg
---

<div class='text-muted' style='margin-bottom:2em'>
<i>This post is co-authored with Denis Yuen, lead of the Dockstore project at the Ontario Institute for Cancer Research</i>
</div>

One key feature of Nextflow is the ability to automatically pull and execute a workflow application directly from a sharing platform such as GitHub. We realised this was critical to allow users to properly track code changes and releases and, above all, to enable the [seamless sharing of workflow projects](/blog/2016/best-practice-for-reproducibility.html).

Nextflow never wanted to implement its own centralised workflow registry because we thought that in order for a registry to be viable and therefore useful, it should be technology agnostic and it should be driven by a consensus among the wider user community.

This is exactly what the [Dockstore](https://dockstore.org/) project is designed for and for this reason we are thrilled to announce that Dockstore has just released the support for Nextflow workflows in its latest release!

### Dockstore in a nutshell

Dockstore is an open platform that collects and catalogs scientific data analysis tools and workflows, starting from the genomics community. It’s developed by the [OICR](https://oicr.on.ca/) in collaboration with [UCSC](https://ucscgenomics.soe.ucsc.edu/) and it is based on the [GA4GH](https://www.ga4gh.org/) open standards and the FAIR principles i.e. the idea to make research data and applications findable, accessible, interoperable and reusable ([FAIR](https://www.nature.com/articles/sdata201618)).

<img src='/img/dockstore.png' alt="Dockstore logo" style='float:right; width: 150pt; padding: .5em;' />


In Dockstore’s initial release of support for Nextflow, users will be able to register and display Nextflow workflows. Many of Dockstore’s cross-language features will be available such as [searching](https://dockstore.org/search?descriptorType=nfl&searchMode=files), displaying metadata information on authorship from Nextflow’s config ([author and description](https://www.nextflow.io/docs/latest/config.html?highlight=author#scope-manifest)), displaying the [Docker images](https://dockstore.org/workflows/github.com/nf-core/hlatyping:1.1.1?tab=tools) used by a workflow, and limited support for displaying a visualization of the [workflow structure](https://dockstore.org/workflows/github.com/nf-core/hlatyping:1.1.1?tab=dag).

The Dockstore team will initially work to on-board the high-quality [nf-core](https://github.com/nf-core) workflows curated by the Nextflow community. However, all developers that develop Nextflow workflows will be able to login, contribute, and maintain workflows starting with our standard [workflow tutorials](https://docs.dockstore.org/docs/publisher-tutorials/workflows/).

Moving forward, the Dockstore team hopes to engage more with the Nextflow community and integrate Nextflow code in order to streamline the process of publishing Nextflow workflows and draw better visualizations of Nextflow workflows. Dockstore also hopes to work with a cloud vendor to add browser based launch-with support for Nextflow workflows.

Finally, support for Nextflow workflows in Dockstore will also enable the possibility of cloud platforms that implement [GA4GH WES](https://github.com/ga4gh/workflow-execution-service-schemas) to run Nextflow workflows.


### Conclusion

We welcome the support for Nextflow workflows in the Dockstore platform. This is a valuable contribution and presents great opportunities for workflow developers and the wider scientific community.

We invite all Nextflow developers to register their data analysis applications in the Dockstore platform to make them accessible and reusable to a wider community of researchers.
