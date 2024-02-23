title=Interview: Josh Chorlton
episode=32
description=Infra and backend for Nextflow and MultiQC in clinical practice.
date=2024-03-05
type=podcast
subtype=Interview
youtubeid=Gjw7QC8cP_M
image=img/podcast_ep32.png
tags=nextflow,multiqc,opensource
status=published
author=Developer advocates
icon=logo_podcast_channels.jpg
~~~~~~

In this episode of Channels, [Phil Ewels](https://twitter.com/tallphil) talks to [Josh Chorlton](https://joshchorlton.com/) - CTO & cofounder of [BugSeq Bioinformatics Inc](https://bugseq.com/).

We talk about how BugSeq got into using Nextflow and MultiQC and the tips and tricks that they've employed to push scale and performace of their tools to the limit.

<!-- end-archive-description -->

Josh leads engineering at BugSeq, serving hundreds of public health and clinical labs globally. He oversees infrastructure that processes thousands of samples every month, in multiple regulatory jurisdictions, for labs of all sizes. He leads deployment for a suite of automated pipelines, with an emphasis on correctness, speed and security. Josh brings to BugSeq a decade of experience managing critical infrastructure at Silicon Valley companies including Stripe, Snap, Uber and Docker. He has a passion for solving challenging, impactful problems with scalable technology to ultimately leave the world in a better state.

...you can see why we wanted to interview him!

In this deeply technical chat, we cover topics like testing and CI/CD, passing around structured data objects between Groovy and Python using Protobufs, building Docker images at scale and the ins and outs of using MultiQC for custom clinical data reports.

* Testing / CI strategies
    * Types of errors and various error strategies
    * BugSeq uses [pytest-workflow](https://pytest-workflow.readthedocs.io/en/stable/) for CI testing
    * Blog post mentioned by Josh: [Providing Reliable, Reproducible and Valid Results with Bioinformatic Versioning](https://docs.bugseq.com/blog/2023/11/15/providing-reliable-reproducible-and-valid-results-with-bioinformatic-versioning/)
* Using protobufs to share constructs between languages
    * [Protobufs](https://protobuf.dev/) (Protocol buffers) are language-neutral, platform-neutral extensible mechanisms for serializing structured data
    * Jake talks us through how protobufs allow BugSeq us to share type-checked objects between our Nextflow code and the python scripts within their pipeline processes
* Docker images
    * Strategy for building hundreds of custom images, with thorough testing and CI
    * BugSeq now uses [Bazel](https://bazel.build/) to help disentangle the work of solving conda envs, building base libraries, and so on during Docker image creation
* MultiQC
    * BugSeq uses MultiQC as the bridge between bioinformaticians and customers
    * Extensive use of [Custom Content](https://multiqc.info/docs/custom_content/) allows custom pipeline outputs to be presented alongside QC metrics parsed by MultiQC from standard bioinformatics tools
* Running bioinformatics analyses as a service
    * BugSeq handles troubleshooting and customer requests with an engineer pager rotation system
    * Customers expect turnaround time within 1-2 hours, so production needs to be rock-solid. But they also expect rapid new features!
    * The CI/CD testing infrastructure is critical to allow BugSeq to ship rapidly with confidence

Huge thanks to Josh for joining us on the podcast. Let's hope we can have him back in early 2025 and see where we're getting to on his wishlist for Nextflow and MultiQC!

Apologies to viewers for the funky camera focus issues. Phil will try to get that sorted before next time.
