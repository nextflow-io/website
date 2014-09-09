title=Reproducibility in Science - Nextflow meets Docker 
date=2014-09-09
type=post
tags=docker,github,reproducibility
status=draft
author=Maria Chatzou
icon=maria.png
~~~~~~
The scientific world nowadays operates on the basis of published articles. 
These are used to report novel discoveries to the rest of the scientific community. 

But have you ever wondered what a scientific article is? It is a:

1. defeasible argument for claims, supported by 	
2. exhibited, reproducible data and methods, and	
3. explicit references to other work in the domain;	
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
In fact, for only 54% of the papers investigated by Vasilevsky and his colleagues the resources 
could be obtained or re-created based on information provided in the publication [[2](https://peerj.com/articles/148/)].

Encouragingly scientific reproducibility has been at the forefront of many news stories 
and there exist numerous initiatives to help address this problem. Particularly, when it 
comes to producing reproducible computational analyses, some publications are starting 
to actually publish the code and data used for analysis and to generate figures. 

For example, many articles in Nature and in the new Elife journal (and others) provide a 
"source data" download link next to figures. Sometimes Elife might even have an option 
to download the source code for figures.

This is a good start, but there are still lots of problems. For example, if one wants 
to re-execute a data analyses from these papers, he/she will have to download the 
scripts and data, to only realize that he/she has not all the required libraries, 
or that it only runs on, for example, an Ubuntu version he/she doesn't have, or some 
paths are hard-coded to match the authors' machines. 

If it's not easy to run and doesn't run "out of the box" the chances that a researcher 
will actually ever run most of these scripts is close to zero, especially if they lack 
the time and expertise to manage the required installations of third-party libraries, 
tools or implement from scratch state-of-the-art data processing algorithms.

### Here it comes Docker

This is why we think that Docker containers technology represent a big opportunity for 
published code and data to be useful. It not only has to be reproducible, 
but it also has to be easy to reproduce and interact with the data. 

[Docker](http://www.docker.com) containers technology, is a solution to many of the coumputational research reproducibility problems. 
Basically, it is kind of like a lightweight virtual machine where you can set up a compute 
environment, including all the libraries, code, data, you need, in a single "image". 

That image can be distributed publicly and can seamlessly run on theoretically 
any major Linux operating system. No need for the user to mess with installation, paths, etc. 

They just run the Docker image you provided, and everything is set up to run out of the box.
Many others (e.g. [here](http://www.bioinformaticszen.com/post/reproducible-assembler-benchmarks/) 
and [here](https://bcbio.wordpress.com/2014/03/06/improving-reproducibility-and-installation-of-genomic-analysis-pipelines-with-docker/)) 
have already started discussing this.
 
### Docker and Nextflow: a perfect match 

A big advantage of Docker containers compared to *traditional* machine virtualisation technology, 
is that it has a minimal startup time. This make it possible to virtualise single application 
executions which can be launched in parallel in their own isolated process container, in 
order to speedup a large computation. 

Nextflow is data-driven toolkit for computational pipelines which aims to simplify the deployment of 
distributed and highly parallelised pipelines for scientific application.  
  
The latest version integrates the support for Docker containers that enables the deployment 
of self-contained and truly reproducible pipelines. 
  
### How they work together 

A Nextflow pipeline is made up by putting together several processes. Each process 
can be written in any scripting language that can be executed by the Linux platform 
(BASH, Perl, Ruby, Python, etc). Parallelisation is automatically managed 
by the framework and it is implicitly defined by the processes input and 
output declarations. 

By integrating Docker with Nextflow, a pipeline processes can be executed independently
in its own isolated container, this guarantee that each of them run in a predictable 
manner without worry about the configuration of the target execution platform. Moreover the minimal overhead added by Docker allow us to spawn multiple containers execution in parallel with a negligible performance loss if compared to a platform *native* execution.   


### An example

As a proof of concept of the Nextflow integration with Docker you can try out the 
following example pipeline at this link 

https://github.com/nextflow-io/examples/blob/master/small_parallel.nf

It splits a FASTA file chunks on sequences, execute a BLAST query for each of this 
chunks, then is a separate process extract the top 10 matching sequences and 
finally aligns all the results with the T-Coffee multiple sequence aligner. 

In order to run it you will need to have a Java VM (7+) and Docker (1.0+) installed in 
your computer. 

if you already don't have it, install the latest version of Nextflow by entering the 
following command in your shell terminal:
   
     curl -fsSL get.nextflow.io | bash
     
     
Then download the Docker image used by this example with the following command: 

     docker pull nextflow/examples
     
Note: can checkout the content of this image at this [link](https://github.com/nextflow-io/examples/blob/master/Dockerfile)

Now you are ready to the demo bu entering the following command: 

    nextflow run examples/small_parallel.nf -with-docker
    
    
This will launch the pipeline execution printing the final alignment out in the terminal screen. 
You can also provide your own protein sequence FASTA file by adding in the above command line 
the option ``--query <file>`` and/or change the splitting chunk size with ``--chunk n``	 


### Conclusion 

This mix of technologies provide Docker, GitHub and Nextflow makes it possible to write 
self-contained and truly reproducible pipeline which requires zero configuration and 
be reproduced in any system having a Java VM and Docker engine installed.


### Learn exactly how to do it!

Follow our documentation guidelines for a quick start with Docker and Nextflow:




       
