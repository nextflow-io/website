title=Running scientific data analysis with Docker  
date=2014-09-09
type=post
tags=docker,github,reproducibility,data-analysis
status=draft
author=Paolo Di Tommaso
icon=pablo.jpg
~~~~~~

Scientific data analysis pipelines are rarely composed by a single piece of software. 
In a real world scenario, computational pipelines are made up of multiple stages, each of which 
can execute many different scripts, system commands and external tools deployed in the hosting computing environment, usually an HPC cluster. 

Working as a research engineer in a bioinformatics lab I've direct experience on how 
difficult is keeping such a piece of software consistent. 

Computing environments are updated frequently, to test new piece of software or because 
a new version of the same is required. For this reason it's almost impossible to replicate 
data analysis results over the time. 
    

Docker has emerged recently as a new type of virtualisation technology that allows one
to create a self-contained, executable computing environment. There are plenty of examples 
showing the benefits of using Docker to run application services, like web servers 
or databases, however it looks that still few haven taken it into consideration to "standardise" 
the deployment of scientific data analysis in a cluster of computers, i.e. a deployment scenario
commonly available in a research institute in order to simply the development, the deployment 
and the replicability of computational pipelines. 

For this reason I decided to evaluate Docker in the cluster available in our lab institute.  

## Method

The Docker engine has been installed in each node in our cluster. A Docker private registry 
instance has also been installed in our internal network so that images can be pulled from the 
local repository in a much faster way when compared to the public [Docker registry](http://registry.hub.docker.com). 

Moreover the [Univa grid engine](http://www.univa.com/products/grid-engine.php) has been configured with custom [complex](http://www.gridengine.eu/mangridengine/htmlman5/complex.html) 
resource type. This allows us to request a specific Docker image as a resource type while 
submitting a job execution to the cluster. 

The Docker image is requested as a `soft` resource, by doing that the UGE scheduler 
tries to run the job [in] to a node where that image has already been already pulled, 
otherwise a lower priority is given to it and it is executed, eventually, by a node where 
the specified Docker image is not available. This will force the node to pull the image 
from the local registry at the time of the job execution. 

We tested it with [Piper-NF](https://github.com/cbcrg/piper-nf), a genomic pipeline for the 
detection and mapping of long non-coding RNAs. 

The pipeline runs on top of Nextflow, which takes care of the tasks parallelisation and submits 
the jobs for execution to the Univa grid engine. 

The Piper-NF code wasn't modified in order to run it by using Docker.
Nextflow is able to handle it automatically. The container is run in such a way that 
the tasks result files are created in the hosting file system, in other 
words it behaves in a completely transparent manner without requiring extra steps or affecting 
the flow of the pipeline execution. 

It was only necessary to specify the Docker image (or images) to be used in the Nextflow 
configuration file for the pipeline. You can read more about this at [this link](http://www.nextflow.io/docs/latest/docker.html).

## Results 

(?) To benchmark how Docker affects the performance of the pipeline execution I compared it to 
the same pipeline executed with the same dataset, without using Docker "virtualisation".

The pipeline execution launches around 100 jobs, 10 cluster nodes were used. I ran it 5 times 
by using the same dataset with and without Docker. 

The average execution time without Docker was 28.6 minutes, while the average 
pipeline executing time, running each job in a Docker container, was 32.2 minutes.
Thus, by using Docker the overall execution time increased by something around 12%. 

It is important to note that this time includes both the Docker bootstrap time, 
and the time overhead that is added to the task execution by its virtualisation layer. 

For this reason also the actual task run time was measured i.e. without including the Docker bootstrap 
time overhead. In this case, the aggregate average task execution time was 57.3 minutes, while 
the time was 59.5 minutes when running the same tasks by using Docker. Thus, the time overhead 
added by the Docker virtualisation layer to the actual tasks run time can be estimated to around 4%
in our pipeline.

Keeping the complete toolset required by the pipeline execution into a Docker image dramatically
reduced configuration and deployment problems. Also storing these images into the private and 
[public repository](https://registry.hub.docker.com/repos/cbcrg/) with a unique tag allowed us 
to replicate the result with the identical computing environment without the usual burden 
of setting-up the required computing environment.    


## Conclusion 

The fast start-up time for Docker containers technology allows one to virtualise single a process or 
a bunch of applications, instead of a complete operating system. This opens up new possibilities, 
for example the possibility to "virtualise" distributed job executions in an HPC computers cluster. 

The minimal performance loss introduced by the Docker engine virtualisation is offset by 
the advantages of running your analysis in a self-contained and dead easy to reproduce computing environment, 
which guarantees the consistency of the results over time and across different computing platforms. 
  

