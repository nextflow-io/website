title=Using Docker for scientific data analysis in an HPC cluster   
date=2014-11-06
type=post
tags=docker,reproducibility,data-analysis,hpc
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~

Scientific data analysis pipelines are rarely composed by a single piece of software. 
In a real world scenario, computational pipelines are made up of multiple stages, each of which 
can execute many different scripts, system commands and external tools deployed in a hosting computing 
environment, usually an HPC cluster. 

As I work as a research engineer in a bioinformatics lab I experience on a daily basis the 
difficulties related on keeping such a piece of software consistent. 

Computing environments can change frequently in order to test new pieces of software or 
maybe because system libraries need to be updated. For this reason replicating the results 
of a data analysis over time can be a challenging task. 
    
[Docker](http://www.docker.com) has emerged recently as a new type of virtualisation technology that allows one
to create a self-contained runtime environment. There are plenty of examples 
showing the benefits of using it to run application services, like web servers 
or databases. 

However it seems that few people have considered using Docker for the deployment of scientific 
data analysis pipelines on distributed cluster of computer, in order to simplify the development, 
the deployment and the replicability of this kind of applications.  

For this reason I wanted to test the capabilities of Docker to solve these problems in the 
cluster available in our [institute](http://www.crg.eu).

## Method

The Docker engine has been installed in each node of our cluster, that runs a [Univa grid engine](http://www.univa.com/products/grid-engine.php) resource manager. 
A Docker private registry instance has also been installed in our internal network, so that images 
can be pulled from the local repository in a much faster way when compared to the public 
[Docker registry](http://registry.hub.docker.com). 

Moreover the Univa grid engine has been configured with a custom [complex](http://www.gridengine.eu/mangridengine/htmlman5/complex.html) 
resource type. This allows us to request a specific Docker image as a resource type while 
submitting a job execution to the cluster. 

The Docker image is requested as a *soft* resource, by doing that the UGE scheduler 
tries to run a job to a node where that image has already been pulled, 
otherwise a lower priority is given to it and it is executed, eventually, by a node where 
the specified Docker image is not available. This will force the node to pull the required 
image from the local registry at the time of the job execution. 

This environment has been tested with [Piper-NF](https://github.com/cbcrg/piper-nf), a genomic pipeline for the 
detection and mapping of long non-coding RNAs. 

The pipeline runs on top of Nextflow, which takes care of the tasks parallelisation and submits 
the jobs for execution to the Univa grid engine. 

The Piper-NF code wasn't modified in order to run it using Docker.
Nextflow is able to handle it automatically. The Docker containers are run in such a way that 
the tasks result files are created in the hosting file system, in other 
words it behaves in a completely transparent manner without requiring extra steps or affecting 
the flow of the pipeline execution. 

It was only necessary to specify the Docker image (or images) to be used in the Nextflow 
configuration file for the pipeline. You can read more about this at [this link](https://www.nextflow.io/docs/latest/docker.html).

## Results 

To benchmark the impact of Docker on the pipeline performance a comparison was made running 
it with and without Docker.

For this experiment 10 cluster nodes were used. The pipeline execution launches around 100 jobs, 
and it was run 5 times by using the same dataset with and without Docker. 

The average execution time without Docker was 28.6 minutes, while the average 
pipeline execution time, running each job in a Docker container, was 32.2 minutes.
Thus, by using Docker the overall execution time increased by something around 12.5%. 

It is important to note that this time includes both the Docker bootstrap time, 
and the time overhead that is added to the task execution by the virtualisation layer. 

For this reason the actual task run time was measured as well i.e. without including the 
Docker bootstrap time overhead. In this case, the aggregate average task execution time was 57.3 minutes
and 59.5 minutes when running the same tasks using Docker. Thus, the time overhead 
added by the Docker virtualisation layer to the effective task run time can be estimated 
to around 4% in our test.

Keeping the complete toolset required by the pipeline execution within a Docker image dramatically
reduced configuration and deployment problems. Also storing these images into the private and 
[public](https://registry.hub.docker.com/repos/cbcrg/) repositories with a unique tag allowed us 
to replicate the results without the usual burden required to set-up an identical computing environment.    


## Conclusion 

The fast start-up time for Docker containers technology allows one to virtualise a single process or 
the execution of a bunch of applications, instead of a complete operating system. This opens up new possibilities, 
for example the possibility to "virtualise" distributed job executions in an HPC cluster of computers. 

The minimal performance loss introduced by the Docker engine is offset by the advantages of running 
your analysis in a self-contained and dead easy to reproduce runtime environment, which guarantees 
the consistency of the results over time and across different computing platforms. 
  
  
#### Credits 

Thanks to Arnau Bria and the all scientific systems admins team to manage the Docker installation 
in the CRG computing cluster. 

