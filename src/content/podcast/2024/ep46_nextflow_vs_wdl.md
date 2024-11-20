---
title: "Nextflow meets WDL: Comparing modern workflow languages"
episode: 46
description: Exploring what Nextflow can learn from WDL
date: 2024-09-17
type: podcast
subtype: Technical discussion
youtubeid: qeJqdV_d8vI
image: /img/podcast_ep46.png
tags: Nextflow
author: Developer advocates
icon: logo_podcast_channels.jpg
---

In this episode of Channels, the Nextflow podcast, host Geraldine Van der Auwera from Seqera and her colleague Ben Sherman dive into the differences and similarities between Nextflow and the [Workflow Description Language](https://openwdl.org/) (WDL). They discuss the origins and development of WDL at the Broad Institute, the challenges of working with different workflow languages, and how Nextflow's channel-based data flow model compares to WDL's approach. The conversation covers key features, type systems, and runtime specifications, highlighting the strengths and limitations of both languages. Ben hints at upcoming improvements and enhancements for Nextflow, aimed at making the language more robust and user-friendly. The episode offers valuable insights for bioinformaticians and developers navigating the complexities of workflow management.

<!-- end-archive-description -->

**Some key points from the episode**

- Tasks vs. processes: basically the same thing
- Channels: why they're awesome once you 'get it' and how to make it easier to leverage their power
- WDL's scatter functions vs. Nextflow's implicit parallelization
- How many operators does Nextflow really need?
- The benefits of a DSL and the appeal of static types
- Don't make me care about what order to pass inputs in
- The evolution of WDL's runtime block vs Nextflow's config files
- Portability across infrastructure backends: who wore it best?
- Cloud-native vs. HPC-to-cloud

## Podcast overview

### Introduction

The comparison between Nextflow and WDL (Workflow Description Language) is a fascinating one for anyone involved in bioinformatics and data science workflows. Both languages have their unique strengths and offer different approaches to solving similar problems. This blog post delves into the intricacies of both languages, highlighting their differences, similarities, and areas of potential improvement.

### Meet the Contributors

Ben Sherman, a software engineer at Seqera and one of the main developers of Nextflow, brings a wealth of experience to the discussion. Geraldine, Lead Developer Advocate at Seqera, offers insights from her extensive past use of WDL. Together, they provide a comparative overview of these two powerful workflow languages.

### Introduction to WDL

WDL, or Workflow Description Language, was developed at the Broad Institute starting around 2014. This language emerged during a period when existing workflow solutions were struggling to scale with increasing data volumes and the need for greater reproducibility. The Broad Institute, a nonprofit research institution in Cambridge, MA, needed a more effective software stack for its genome analysis pipeline, which was moving to the cloud.
This led to the development of WDL as a user-friendly, cloud-native DSL (domain-specific language) designed for high reproducibility. Today, WDL is an open-source project supported by two main engines: Cromwell (focused on high-scale execution) and MiniWDL (geared towards prototyping).

### WDL Development History

From its inception, WDL was intended to simplify the process of sharing genomic analysis pipelines. With WDL, researchers can share readable, user-friendly scripts, with some limitations. As it evolved, WDL introduced more flexibility while keeping complexity in check. This balance is crucial for maintaining an approachable yet powerful workflow language.

### Key Differences Between WDL and Nextflow

One of the major differences between Nextflow and WDL lies in their approach to parallelization and data handling. WDL uses explicit task invocation and a scatter function for batch processing. This approach is straightforward but can be limiting in complex scenarios.
On the other hand, Nextflow relies on channels to handle data flow between processes. This abstraction simplifies the user’s code and allows for more compact, readable scripts. Channels can hold single values, lists, or result sets, enabling processes to consume and produce data efficiently. While this implicit handling of data flow can be challenging to learn, it offers significant rewards in process efficiency and code simplicity.

### Advanced Features of Nextflow

Nextflow processes and channels set it apart from most workflow managers. Inspired by Unix’s pipeline philosophy, Nextflow allows for powerful data manipulation through operators like map, flatMap, and filter. These operators enable users to create flexible, scalable workflows with fewer lines of code.
However, Nextflow's power can also lead to over-engineering in inexperienced hands. The development team continually strives to simplify Nextflow while maintaining its flexibility. For example, integrating static typing and record types could ease the learning curve by clarifying data structures and inputs.

### Challenges and Improvements in Nextflow

Many of Nextflow's operators might be rethought or merged with standard library functions to simplify the language. Adding features like named parameter passing, similar to WDL, could enhance user-friendliness. Balancing power and simplicity is key to preventing users from creating unnecessarily complex workflows.

### The Future of Nextflow and WDL

Both Nextflow and WDL have come a long way and continue to evolve. Lessons from WDL's strong typing and static typing system offer valuable insights for Nextflow's development. Borrowing best practices from each other can drive future improvements.

### Conclusion

The comparison between Nextflow and WDL highlights the strengths of both languages and areas for potential enhancement. While WDL excels in its simplicity and user-friendliness, Nextflow shines in its flexibility and powerful data manipulation capabilities. As both languages continue to develop, they promise to offer even more robust solutions for bioinformatics and data science workflows. Stay tuned for more updates and innovations in both Nextflow and WDL!

## Resources for newcomers to Nextflow

If you're new to Nextflow and you'd like to give it a try, we recommend starting with the [Hello Nextflow tutorial series](https://training.nextflow.io/hello_nextflow/), which will guide you through the process of writing and running your first Nextflow pipelines. This web-based training uses Gitpod and does not require installing anything on your own computer.
