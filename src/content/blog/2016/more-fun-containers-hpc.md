---
title: More fun with containers in HPC
date: 2016-12-20
type: post
tags: aws,pipelines,nextflow,genomic,docker,singularity
author: Paolo Di Tommaso
icon: paolo.jpg
---

Nextflow was one of the [first workflow framework](https://www.nextflow.io/blog/2014/nextflow-meets-docker.html)
to provide built-in support for Docker containers. A couple of years ago we also started
to experiment with the deployment of containerised bioinformatic pipelines at CRG,
using Docker technology (see [here]((https://www.nextflow.io/blog/2014/using-docker-in-hpc-cluster.html)) and [here](https://www.nextplatform.com/2016/01/28/crg-goes-with-the-genomics-flow/)).

We found that by isolating and packaging the complete computational workflow environment
with the use of Docker images, radically simplifies the burden of maintaining complex
dependency graphs of real workload data analysis pipelines.

Even more importantly, the use of containers enables replicable results with minimal effort
for the system configuration. The entire computational environment can be archived in a
self-contained executable format, allowing the replication of the associated analysis at
any point in time.

This ability is the main reason that drove the rapid adoption of Docker in the bioinformatic
community and its support in many projects, like for example [Galaxy](https://galaxyproject.org),
[CWL](http://commonwl.org), [Bioboxes](http://bioboxes.org), [Dockstore](https://dockstore.org) and many others.

However, while the popularity of Docker spread between the developers, its adaption in
research computing infrastructures continues to remain very low and it's very unlikely
that this trend will change in the future.

The reason for this resides in the Docker architecture, which requires a daemon running
with root permissions on each node of a computing cluster. Such a requirement raises many
security concerns, thus good practices would prevent its use in shared HPC cluster or
supercomputer environments.

### Introducing Singularity

Alternative implementations, such as [Singularity](http://singularity.lbl.gov), have
fortunately been promoted by the interested in containers technology.

Singularity is a containers engine developed at the Berkeley Lab and designed for the
needs of scientific workloads. The main differences with Docker are: containers are file
based, no root escalation is allowed nor root permission is needed to run a container
(although a privileged user is needed to create a container image), and there is no
separate running daemon.

These, along with other features, such as support for autofs mounts, makes Singularity a
container engine better suited to the requirements of HPC clusters and supercomputers.

Moreover, although Singularity uses a container image format different to that of Docker,
they provide a conversion tool that allows Docker images to be converted to the
Singularity format.

### Singularity in the wild

We integrated Singularity support in Nextflow framework and tested it in the CRG
computing cluster and the BSC [MareNostrum](https://www.bsc.es/discover-bsc/the-centre/marenostrum) supercomputer.

The absence of a separate running daemon or image gateway made the installation
straightforward when compared to Docker or other solutions.

To evaluate the performance of Singularity we carried out the [same benchmarks](https://peerj.com/articles/1273/)
we performed for Docker and compared the results of the two engines.

The benchmarks consisted in the execution of three Nextflow based genomic pipelines:

1. [Rna-toy](https://github.com/nextflow-io/rnatoy/tree/peerj5515): a simple pipeline for RNA-Seq data analysis.
2. [Nmdp-Flow](https://github.com/nextflow-io/nmdp-flow/tree/peerj5515/): an assembly-based variant calling pipeline.
3. [Piper-NF](https://github.com/cbcrg/piper-nf/tree/peerj5515): a pipeline for the detection and mapping of long non-coding RNAs.

In order to repeat the analyses, we converted the container images we used to perform
the Docker benchmarks to Singularity image files by using the [docker2singularity](https://github.com/singularityware/docker2singularity) tool
*(this is not required anymore, see the update below)*.

The only change needed to run these pipelines with Singularity was to replace the Docker
specific settings with the following ones in the configuration file:

    singularity.enabled = true
    process.container = '<the image file path>'

Each pipeline was executed 10 times, alternately by using Docker and Singularity as
container engine. The results are shown in the following table (time in minutes):

<style>
table#benchmark { width: 100%; margin-bottom: 1em; margin-top: 1em; border-top: 1px solid #999; border-bottom: 1px solid #999 }
table#benchmark th { text-align: center; background-color: #eee; padding: 2px 5px }
table#benchmark td { text-align: right; padding: 2px 5px; padding-right: 15px }
table#benchmark .r { text-align: right }
table#benchmark .l { text-align: left }
</style>
<table id='benchmark'>
<tr>
<th class='l'>Pipeline</th>
<th >Tasks</th>
<th colspan=2 >Mean task time</th>
<th colspan=2 >Mean execution time</th>
<th colspan=2 >Execution time std dev</th>
<th colspan=2 >Ratio</th>
</tr>

<tr>
<th>&nbsp;</th>
<th>&nbsp;</th>
<th>Singularity</th>
<th>Docker</th>
<th>Singularity</th>
<th>Docker</th>
<th>Singularity</th>
<th>Docker</th>
<th>&nbsp;</th>
</tr>

<tr>
<td class='l'>RNA-Seq</td>
<td>9</td>
<td>73.7</td>
<td>73.6</td>
<td>663.6</td>
<td>662.3</td>
<td>2.0</td>
<td>3.1</td>
<td>0.998</td>
</tr>

<tr>
<td class='l'>Variant call</td>
<td>48</td>
<td>22.1</td>
<td>22.4</td>
<td>1061.2</td>
<td>1074.4</td>
<td>43.1</td>
<td>38.5</td>
<td>1.012</td>
</tr>

<tr>
<td class='l'>Piper-NF</td>
<td>98</td>
<td>1.2</td>
<td>1.3</td>
<td>120.0</td>
<td>124.5</td>
<td>6.9 </td>
<td>2.8</td>
<td>1.038</td>
</tr>

</table>


The benchmark results show that there isn't any significative difference in the
execution times of containerised workflows between Docker and Singularity. In two
cases Singularity was slightly faster and a third one it was almost identical although
a little slower than Docker.


### Conclusion

In our evaluation Singularity proved to be an easy to install,
stable and performant container engine.

The only minor drawback, we found when compared to Docker, was the need to define the
host path mount points statically when the Singularity images were created. In fact,
even if Singularity supports user mount points to be defined dynamically when the
container is launched, this feature requires the overlay file system which was not
supported by the kernel available in our system.

Docker surely will remain the *de facto* standard engine and image format for containers
due to its popularity and [impressive growth](http://www.coscale.com/blog/docker-usage-statistics-increased-adoption-by-enterprises-and-for-production-use).

However, in our opinion, Singularity is the tool of choice for the execution of
containerised workloads in the context of HPC, thanks to its focus on system security
and its simpler architectural design.

The transparent support provided by Nextflow for both Docker and Singularity technology
guarantees the ability to deploy your workflows in a range of different platforms (cloud,
cluster, supercomputer, etc). Nextflow transparently manages the deployment of the
containerised workload according to the runtime available in the target system.


#### Credits

Thanks to Gabriel Gonzalez (CRG), Luis Exposito (CRG) and Carlos Tripiana Montes (BSC)
for the support installing Singularity.


**Update** Singularity, since version 2.3.x, is able to pull and run Docker images from the Docker Hub.
This greatly simplifies the interoperability with existing Docker containers. You only need
to prefix the image name with the `docker://` pseudo-protocol to download it as a Singularity image,
for example:

    singularity pull --size 1200 docker://nextflow/rnatoy
