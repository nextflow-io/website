title=Nextflow Summit 2022 Recap
date=2022-11-03
type=post
description=Catch up on exciting developments from this year’s Nextflow Summit held in Barcelona
image=img/nextflow-summit-2022-recap.jpg
tags=nextflow,tower,cloud
status=published
author=Noel Ortiz
icon=noel.jpg
~~~~~~

## Three days of Nextflow goodness in Barcelona

After a three-year COVID-related hiatus from in-person events, Nextflow developers and users found their way to Barcelona this October for the 2022 Nextflow Summit. Held at Barcelona’s iconic Agbar tower, this was easily the most successful Nextflow community event yet!

The week-long event kicked off with 50 people participating in a hackathon organized by nf-core beginning on October 10th. The [hackathon](https://nf-co.re/events/2022/hackathon-october-2022) tackled several cutting-edge projects with developer teams focused on various aspects of nf-core including documentation, subworkflows, pipelines, DSL2 conversions, modules, and infrastructure. The Nextflow Summit began mid-week attracting nearly 600 people, including 165 attending in person and another 433 remotely. The [YouTube live streams](https://summit.nextflow.io/stream/) have now collected over two and half thousand views. Just prior to the summit, three virtual Nextflow training events were also run with separate sessions for the Americas, EMEA, and APAC in which 835 people participated.

## An action-packed agenda

The three-day Nextflow Summit featured 33 talks delivered by speakers from academia, research, healthcare providers, biotechs, and cloud providers. This year’s speakers came from the following organizations:

- Amazon Web Services
- Center for Genomic Regulation
- Centre for Molecular Medicine and Therapeutics, University of British Columbia
- Chan Zukerberg Biohub
- Curative
- DNA Nexus
- Enterome
- Google
- Janelia Research Campus
- Microsoft
- Oxford Nanopore
- Quadram Institute BioScience
- Seqera Labs
- Quantitative Biology Center, University of Tübingen
- Quilt Data
- UNC Lineberger Comprehensive Cancer Center
- Università degli Studi di Macerata
- University of Maryland
- Wellcome Sanger Institute
- Wyoming Public Health Laboratory

## Some recurring themes

While there were too many excellent talks to cover individually, a few themes surfaced throughout the summit. Not surprisingly, SARS-Cov-2 was a thread that wound through several talks. Tony Zeljkovic from Curative led a discussion about [unlocking automated bioinformatics for large-scale healthcare](https://www.youtube.com/watch?v=JZMaRYzZxGU&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=8), and Thanh Le Viet of Quadram Institute Bioscience discussed [large-scale SARS-Cov-2 genomic surveillance at QIB](https://www.youtube.com/watch?v=6jQr9dDaais&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=30). Several speakers discussed best practices for building portable, modular pipelines. Other common themes were data provenance & traceability, data management, and techniques to use compute and storage more efficiently. There were also a few talks about the importance of dataflows in new application areas outside of genomics and bioinformatics.

## Data provenance tracking

In the Thursday morning keynote, Rob Patro﹘Associate Professor at the University of Maryland Dept. of Computer Science and CTO and co-founder of Ocean Genomics﹘described in his talk “[What could be next(flow)](https://www.youtube.com/watch?v=vNrKFT5eT8U&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=6),” how far the Nextflow community had come in solving problems such as reproducibility, scalability, modularity, and ease of use. He then challenged the community with some complex issues still waiting in the wings. He focused on data provenance as a particularly vexing challenge explaining how tremendous effort currently goes into manual metadata curation.

Rob offered suggestions about how Nextflow might evolve, and coined the term “augmented execution contexts” (AECs) drawing from his work on provenance tracking – answering questions such as “what are these files, and where did they come from.” This thinking is reflected in [tximeta](https://github.com/mikelove/tximeta), a  project co-developed with Mike Love of UNC. Rob also proposed ideas around automating data format conversions analogous to type casting in programming languages explaining how such conversions might be built into Nextflow channels to make pipelines more interoperable.

In his talk with the clever title “[one link to rule them all](https://www.youtube.com/watch?v=dttkcuP3OBc&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=13),” Aneesh Karve of Quilt explained how every pipeline run is a function of the code, environment, and data, and went on to show how Quilt could help dramatically simplify data management with dataset versioning, accessibility, and verifiability. Data provenance and traceability were also front and center when Yih-Chii Hwang of DNAnexus described her team’s work around [bringing GxP compliance to Nextflow workflows](https://www.youtube.com/watch?v=RIwpJTDlLiE&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=21).

## Data management and storage

Other speakers also talked about challenges related to data management and performance. Angel Pizarro of AWS gave an interesting talk comparing the [price/performance of different AWS cloud storage options](https://www.youtube.com/watch?v=VXtYCAqGEQQ&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=12). [Hatem Nawar](https://www.youtube.com/watch?v=jB91uqUqsRM&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=9) (Google) and [Venkat Malladi](https://www.youtube.com/watch?v=GAIL8ZAMJPQ&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=20) (Microsoft) also talked about cloud economics and various approaches to data handling in their respective clouds. Data management was also a key part of Evan Floden’s discussion about Nextflow Tower where he discussed Tower Datasets, as well as the various cloud storage options accessible through Nextflow Tower. Finally, Nextflow creator Paolo Di Tommaso unveiled new work being done in Nextflow to simplify access to data residing in object stores in his talk “[Nextflow and the future of containers](https://www.youtube.com/watch?v=PTbiCVq0-sE&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=14)”.

## Compute optimization

Another recurring theme was improving compute efficiency. Several talks discussed using containers more effectively, leveraging GPUs & FPGAs for added performance, improving virtual machine instance type selection, and automating resource requirements. Mike Smoot of Illumina talked about Nextflow, Kubernetes, and DRAGENs and how Illumina’s FPGA-based Bio-IT Platform can dramatically accelerate analysis. Venkat Malladi discussed efforts to suggest optimal VM types based on different standardized nf-core labels in the Azure cloud (process_low, process_medium, process_high, etc.) Finally, Evan Floden discussed [Nextflow Tower](https://www.youtube.com/watch?v=yJpN3fRSClA&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=22) and unveiled an exciting new [resource optimization feature](https://seqera.io/blog/optimizing-resource-usage-with-nextflow-tower/) that can intelligently tune pipeline resource requests to radically reduce cloud costs and improve run speed. Overall, the Nextflow community continues to make giant strides in improving efficiency and managing costs in the cloud.

## Beyond genomics

While most summit speakers focused on genomics, a few discussed data pipelines in other areas, including statistical modeling, analysis, and machine learning. Nicola Visonà from Università degli Studi di Macerata gave a fascinating talk about [using agent-based models to simulate the first industrial revolution](https://www.youtube.com/watch?v=PlKJ0IDV_ds&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=27). Similarly, Konrad Rokicki from the Janelia Research Campus explained how Janelia are using [Nextflow for petascale bioimaging data](https://www.youtube.com/watch?v=ZjSzx1I76z0&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=18) and why bioimage processing remains a large domain area with an unmet need for reproducible workflows.

## Summit Announcements

This year’s summit also saw several exciting announcements from Nextflow developers. Paolo Di Tommaso, during his talk on [the future of containers](https://www.youtube.com/watch?v=PTbiCVq0-sE&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=14), announced the availability of [Nextflow 22.10.0](https://github.com/nextflow-io/nextflow/releases/tag/v22.10.0). In addition to various bug fixes, the latest Nextflow release introduces an exciting new technology called Wave that allows containers to be built on the fly from Dockerfiles or Conda recipes saved within a Nextflow pipeline. Wave also helps to simplify containerized pipeline deployment with features such as “container augmentation”; enabling developers to inject new container scripts and functionality on the fly without needing to rebuild the base containers such as a cloud-native [Fusion file system](https://www.nextflow.io/docs/latest/fusion.html). When used with Nextflow Tower, Wave also simplifies authentication to various public and private container registries. The latest Nextflow release also brings improved support for Kubernetes and enhancements to documentation, along with many other features.

Several other announcements were made during [Evan Floden’s talk](https://www.youtube.com/watch?v=yJpN3fRSClA&list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32&index=22&t=127s), such as:

- MultiQC is joining the Seqera Labs family of products
- Fusion – a distributed virtual file system for cloud-native data pipelines
- Nextflow Tower support for Google Cloud Batch
- Nextflow Tower resource optimization
- Improved Resource Labels support in Tower with integrations for cost accounting with all major cloud providers
- A new Nextflow Tower dashboard coming soon, providing visibility across workspaces

## Thank you to our sponsors

The summit organizers wish to extend a sincere thank you to the event sponsors: AWS, Google Cloud, Seqera Labs, Quilt Data, Oxford Nanopore Technologies, and Element BioSciences. In addition, the [Chan Zuckerberg Initiative](https://chanzuckerberg.com/eoss/) continues to play a key role with their EOSS grants funding important work related to Nextflow and the nf-core community. The success of this year’s summit reminds us of the tremendous value of community and the critical impact of open science software in improving the quality, accessibility, and efficiency of scientific research.

## Learning more

For anyone who missed the summit, you can still watch the sessions or view the training sessions at your convenience:

- Watch post-event recordings of the [Nextflow Summit on YouTube](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmUdAJlHowo7zL2pN3x97d32)
- View replays of the recent online [Nextflow and nf-core training](https://nf-co.re/events/2022/training-october-2022)

For additional detail on the summit and the preceding nf-core events, also check out an excellent [summary of the event](https://mribeirodantas.xyz/blog/index.php/2022/10/27/nextflow-and-nf-core-hot-news/) written by Marcel Ribeiro-Dantas in his blog, the [Dataist Storyteller](https://mribeirodantas.xyz/blog/index.php/2022/10/27/nextflow-and-nf-core-hot-news/)!
