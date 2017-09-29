title=Nexflow Hackathon 2017
date=2017-09-30
type=post
tags=nextflow,docker,hackathon
status=published
author=Evan Floden 
icon=evan.jpg
~~~~~~

Last week saw the inaugural Nextflow meeting organised at the Centre for Genomic Regulation 
(CRG) in Barcelona. The event combined talks, demos, a tutorial/workshop for beginners as 
well as two hackathons sessions for more advanced users. 

Nearly 50 participants attended over the two days which included an entertaining tapas course 
during the first evening!

One of the main objectives of the event was to bring together Nextflow users to work 
together on common interest projects. Several propositions were proposed for the hackathon 
sessions and in the end five diverse ideas were chosen for communal development ranging from 
developing new pipelines through the addition of new features into Nextflow. 

The outcomes of each the projects, which can be found in the issues section 
of [this GitHub repository](https://github.com/nextflow-io/hack17), have been summarised below.

### Nextflow HTML tracing reports

The HTML tracing project aims to generate a rendered version of the Nextflow trace file to 
enable fast sorting and visualisation of task/process execution statistics. 

Currently the data in the trace includes information such as CPU duration, memory usage and 
completion status of each task, however wading through the file is often not convenient 
when a large number of tasks have been executed. 

[Phil Ewels](https://github.com/ewels) proposed the idea and led the coordination effort 
with the outcome being a very impressive working prototype which can be found in the Nextflow 
branch `html-trace`. 

An  image of the example report is shown below with the interactive HTML available 
[here](/misc/nf-trace-report.html). It is expected to merged into the main branch of Nextflow 
with documentation into a near future release.

<img alt='Nextflow HTML execution report' width='760' src='/img/nf-trace-report-min.png' style='margin:1em auto'/>

### Nextflow pipeline for 16S microbial data

The H3Africa Bioinformatics Network have been developing several pipelines which are used 
across the participating centers. The diverse computing resources available has led to 
members wanting workflow solutions with a particular focus on portability. 

With this is mind, Scott Hazelhurst proposed a project for a 16S Microbial data analysis 
pipeline which had [previously been developed using CWL](https://github.com/h3abionet/h3abionet16S/tree/master). 

The participants made a new [branch](https://github.com/h3abionet/h3abionet16S/tree/nextflow) 
of the original pipeline and ported it into Nextflow. 

The pipeline will continue to be developed with the goal of acting as a comparison between 
CWL and Nextflow. It is thought this can then be extended to other pipelines by both those 
who are already familiar with Nextflow as well as used as a tool in training newer users.

### Nextflow modules prototyping

*Toolboxing* allows users to incorporate software into their pipelines in an efficient and 
reproducible manner. Various software repositories are becoming increasing popular, 
highlighted by the over 5,000 tools available in the [Galaxy toolshed](https://toolshed.g2.bx.psu.edu/). 

Projects such as [Biocontainers](http://biocontainers.pro/) aim to wrap up the execution 
environment using containers. [Myself](https://github.com/skptic) and [Johan Viklund](https://github.com/viklund) 
wished to piggyback  off existing repositories and settled on [Dockstore](https://dockstore.org) 
which is an open platform compliant with the [GA4GH](http://genomicsandhealth.org) initiative. 

The majority of tools in Dockstore are written in the CWL and therefore required a parser 
between the CWL CoomandLineTool class and Nextflow processes. Johan was able to develop 
a parser which generates Nextflow processes for several Dockstore tools. 

As these resources such as Dockstore become mature and standardised, it will be 
possible to automatically generate a *Nextflow Store* and enable efficient incorporation 
of tools into workflows.

[IMAGE]

*Example showing a Nextflow process generated from the Dockstore CWL repository for the tool BAMStats.*

### Nextflow pipeline for de novo assembly of nanopore reads

[Nanopore sequencing](https://en.wikipedia.org/wiki/Nanopore_sequencing) is an exciting 
and emerging technology which promises to change the landscape of nucleotide sequencing. 

With keen interest in Nanopore specific pipelines, [Hadrien Gourlé](https://github.com/HadrienG) 
lead the hackathon project for Nanoflow. 

[Nanoflow](https://github.com/HadrienG/nanoflow) is a de novo assembler of bacterials genomes 
from nanopore reads using Nextflow. 

During the two days the participants developed the pipeline for adapter trimming as well 
as assembly and consensus sequence generation using either 
[Canu](https://github.com/marbl/canu) and [Miniasm](https://github.com/lh3/miniasm). 

The future plans are to finalise the pipeline to include a polishing step and a genome 
annotation step. 

### Nextflow AWS Batch integration

Nextflow already has experimental support for [AWS Batch](https://aws.amazon.com/batch/) 
and the goal of this project proposed by [Francesco Strozzi](https://github.com/fstrozzi) 
was to improve this support, add features and test the implementation on real world pipelines. 

Earlier work from [Paolo Di Tommaso](https://github.com/pditommaso) of the Nextflow 
repository, highlighted several challenges to using AWS Batch with Nextflow. 

The major obstacle described by [Tim Dudgeon](https://github.com/tdudgeon) was the requirement 
for each Docker container to have a version of the Amazon Web Services Command Line tools 
(aws-cli) installed. 

Solution worked on was to install the AWS CLI tools on a custom AWS image that is used as the 
Docker host machine, and then mount the directory that contains the necessary items into 
each of the Docker containers as a volume. Early testing suggests this approach works and 
with hopes of providing a more elegant solution in future iterations. 

The code and documentation for AWS Batch has been prepared and will be tested further 
before being rolled into a official Nextflow release in the near future.

### Conclusion

The event was seen as an overwhelming success and special thanks must be made to all the 
participants. As the Nextflow community continues to grow, it is hoped that these types 
of meetings can we can become regular occasions. 

In the meantime we have put together a short video containing some of the highlights 
of the two days. 

We hope to see you all again in Barcelona or at new events around the world soon!

<iframe width="760" height="428" src="https://www.youtube.com/embed/s7SqYMRiY8w?rel=0" frameborder="0" allowfullscreen></iframe>
