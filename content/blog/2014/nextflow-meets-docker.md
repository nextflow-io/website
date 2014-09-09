title=Reproducibility in Science - Nextflow meets Docker 
date=2014-09-09
type=post
tags=docker,github,reproducibility,data-analysis
status=published
author=Maria Chatzou
icon=maria.png
~~~~~~
The scientific world nowadays operates on the basis of published articles. 
These are used to report novel discoveries to the rest of the scientific community. 

But have you ever wondered what a scientific article is? It is a:

1. defeasible argument for claims, supported by 	
2. exhibited, reproducible data and methods, and	
3. explicit references to other work in that domain;	
4. described using domain-agreed technical terminology,	
5. which exists within a complex ecosystem of technologies, people and activities.

Hence the very essence of Science relies on the ability of scientists to reproduce and 
build upon each other’s published results.

So how much can we rely on published data? In a recent report in Nature, researchers at the 
Amgen corporation found that only 11% of the academic research in the literature was 
reproducible by their groups [[1](http://www.nature.com/nature/journal/v483/n7391/full/483531a.html)]. 

While many factors are likely at play here, perhaps the most basic requirement for 
reproducibility holds that the materials reported in a study can be uniquely identified
and obtained, such that experiments can be reproduced as faithfully as possible. 
This information is meant to be documented in the "materials and methods" of journal articles, 
but as many can attest, the information provided there is often not adequate for this task. 

### Promoting Computational Research Reproducibility

Encouragingly scientific reproducibility has been at the forefront of many news stories 
and there exist numerous initiatives to help address this problem. Particularly, when it 
comes to producing reproducible computational analyses, some publications are starting 
to publish the code and data used for analysing and generating figures. 

For example, many articles in Nature and in the new Elife journal (and others) provide a 
"source data" download link next to figures. Sometimes Elife might even have an option 
to download the source code for figures.

This is a good start, but there are still lots of problems. For example, if one wants 
to re-execute a data analyses from these papers, he/she will have to download the 
scripts and the data, to only realize that he/she has not all the required libraries, 
or that it only runs on, for example, an Ubuntu version he/she doesn't have, or some 
paths are hard-coded to match the authors' machines. 

If it's not easy to run and doesn't run out of the box the chances that a researcher 
will actually ever run most of these scripts is close to zero, especially if they lack 
the time or expertise to manage the required installation of third-party libraries, 
tools or implement from scratch state-of-the-art data processing algorithms.

### Here it comes Docker

[Docker](http://www.docker.com) containers technology is a solution to many of the computational 
research reproducibility problems.  Basically, it is kind of a lightweight virtual machine 
where you can set up a computing environment including all the libraries, code and data that you need, 
into a single *image*. 

This image can be distributed publicly and can seamlessly run on any major Linux operating system. 
No need for the user to mess with installation, paths, etc. 

They just run the Docker image you provided, and everything is set up to work out of the box.
Researchers have already started discussing this (e.g. [here](http://www.bioinformaticszen.com/post/reproducible-assembler-benchmarks/), 
[here](https://bcbio.wordpress.com/2014/03/06/improving-reproducibility-and-installation-of-genomic-analysis-pipelines-with-docker/) 
and [here](http://melissagymrek.com/science/2014/08/29/docker-reproducible-research.html)).
 
### Docker and Nextflow: a perfect match 

A big advantage of Docker containers compared to *traditional* machine virtualisation technology
is that it doesn't need a complete copy of the operating system, thus it has a minimal 
startup time. This makes it possible to virtualise single applications or to launch the execution 
of plenty of containers, that can run in parallel, in order to speedup a large computation. 

Nextflow is a data-driven toolkit for computational pipelines, which aims to simplify the deployment of 
distributed and highly parallelised pipelines for scientific applications.  
  
The latest version integrates the support for Docker containers that enables the deployment 
of self-contained and truly reproducible pipelines. 
  
### How they work together 

A Nextflow pipeline is made up by putting together several processes. Each process 
can be written in any scripting language that can be executed by the Linux platform 
(BASH, Perl, Ruby, Python, etc). Parallelisation is automatically managed 
by the framework and it is implicitly defined by the processes input and 
output declarations. 

By integrating Docker with Nextflow, every pipeline process can be executed independently
in its own container, this guarantee that each of them run in a predictable 
manner without worrying about the configuration of the target execution platform. Moreover the 
minimal overhead added by Docker allow us to spawn multiple containers executions in a parallel 
manner with a negligible performance loss when compared to a platform *native* execution.   


### An example

As a proof of concept of the Docker integration with Nextflow you can try out the 
pipeline example at this [link](https://github.com/nextflow-io/examples/blob/master/blast-parallel.nf). 

It splits a protein sequences' FASTA file in chunks of *n* entries, executes a BLAST query 
for each of them, then extracts the top 10 matching sequences and 
finally aligns the results with the T-Coffee multiple sequence aligner. 

In a common scenario you would need to install and configure the tools required by this 
script: BLAST and T-Coffee. Moreover you should provide a formatted protein database in order
to execute the BLAST search.

By using  Docker with Nextflow you only need to have the Docker engine installed in your 
computer and a Java VM. In order to try this example out follow these steps: 

Install the latest version of Nextflow by entering the following command in your shell terminal:
   
     curl -fsSL get.nextflow.io | bash
     
Then download the required Docker image with this command: 

     docker pull nextflow/examples
     
You can check the content of the image looking at the [Dockerfile](https://github.com/nextflow-io/examples/blob/master/Dockerfile)
used to create it.

Now you are ready to run the demo by launching the pipeline execution as shown below: 

    nextflow run examples/blast-parallel.nf -with-docker
    
    
This will run the pipeline printing the final alignment out in the terminal screen. 
You can also provide your own protein sequences' FASTA file by adding, in the above command line, 
the option ``--query <file>`` and change the splitting chunk size with ``--chunk n`` option. 

Note: the result doesn't have a real biological meaning since it uses a very small protein database. 


### Conclusion 

The mix of Docker, GitHub and Nextflow technologies makes it possible to deploy 
self-contained and truly replicable pipelines. It requires zero configuration and 
enables the reproducibility of data analysis pipelines in any system in which a Java VM and 
the Docker engine are available.


### Learn how to do it!

Follow our documentation for a quick start using Docker with Nextflow at 
the following link http://www.nextflow.io/docs/latest/docker.html




       
