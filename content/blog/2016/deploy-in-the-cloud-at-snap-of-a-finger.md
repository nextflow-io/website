title=Deploy your computational pipelines in the cloud at the snap-of-a-finger 
date=2016-09-01
type=post
tags=aws,cloud,pipelines,nextflow,genomic,docker
status=published
author=Paolo Di Tommaso 
icon=paolo.jpg
~~~~~~

<p class="text-muted" style='font-size: 1.2em; padding-bottom: 10px'>
<i>Learn how to deploy and run a computational pipeline in the Amazon AWS cloud with ease 
thanks to Nextflow and Docker containers</i>
</p>

Nextflow is a framework that simplifies the writing of parallel and distributed computational
pipelines in a portable and reproducible manner across different computing platforms, from 
a laptop to a cluster of computers. 

Indeed, the original idea, when this project started three years ago, was to 
implement a tool that would allow researchers in 
[our lab](http://www.crg.eu/es/programmes-groups/comparative-bioinformatics) to smoothly migrate 
their data analysis applications in the cloud when needed - without having 
to change or adapt their code. 

However to date Nextflow has been used mostly to deploy computational workflows within on-premise 
computing clusters or HPC data-centers, because these infrastructures are easier to use 
and provide, on average, cheaper cost and better performance when compared to a cloud environment. 

A major obstacle to efficient deployment of scientific workflows in the cloud is the lack
of a performant POSIX compatible shared file system. These kinds of applications
are usually made-up by putting together a collection of tools, scripts and 
system commands that need a reliable file system to share with each other the input and 
output files as they are produced, above all in a distributed cluster of computers. 

The recent availability of the [Amazon Elastic File System](https://aws.amazon.com/efs/) 
(EFS), a fully featured NFS based file system hosted on the AWS infrastructure represents 
a major step in this context, unlocking the deployment of scientific computing 
in the cloud and taking it to the next level. 

### Nextflow support for the cloud 

Nextflow could already be deployed in the cloud, either using tools such as 
[ElastiCluster](https://github.com/gc3-uzh-ch/elasticluster) or 
[CfnCluster](https://aws.amazon.com/hpc/cfncluster/), or by using custom deployment 
scripts. However the procedure was still cumbersome and, above all, it was not optimised 
to fully take advantage of cloud elasticity i.e. the ability to (re)shape the computing 
cluster dynamically as the computing needs change over time. 

For these reasons, we decided it was time to provide Nextflow with a first-class support 
for the cloud, integrating the Amazon EFS and implementing an optimised native cloud
scheduler, based on [Apache Ignite](https://ignite.apache.org/), with a full support for cluster 
auto-scaling and spot/preemptible instances. 

In practice this means that Nextflow can now spin-up and configure a fully featured computing 
cluster in the cloud with a single command, after that you need only to login to the master 
node and launch the pipeline execution as you would do in your on-premise cluster. 


### Demo !

Since a demo is worth a thousands words, I've record a short screencast showing how 
Nextflow can setup a cluster in the cloud and mount the Amazon EFS shared file system. 

<script type="text/javascript" src="https://asciinema.org/a/9vupd4d72ivaz6h56pajjjkop.js" id="asciicast-9vupd4d72ivaz6h56pajjjkop" async></script>

<p class="text-muted" style='font-size: 0.9em; position: relative; top:-15px' >
Note: in this screencast it has been cut the Ec2 instances startup delay. It required around 
5 minutes to launch them and setup the cluster.   
</p>
   
Let's recap the steps showed in the demo: 
 
* The user provides the cloud parameters (such as the VM image ID and the instance type) 
  in the `nextflow.config` file. 

* To configure the EFS file system you need to provide your EFS storage ID and the mount path
  by using the `sharedStorageId` and `sharedStorageMount` properties.  
  
* To use [EC2 Spot](https://aws.amazon.com/ec2/spot/) instances, just specify the price 
  you want to bid by using the `spotPrice` property.
  
* The AWS access and secret keys are provided by using the usual environment variables. 

* The `nextflow cloud create` launches the requested number of instances, configures the user and 
  access key, mounts the EFS storage and setups the Nextflow cluster automatically.
  Any Linux AMI can be used, it is only required that the [cloud-init](https://cloudinit.readthedocs.io/en/latest/) 
  package, a Java 7+ runtime and the Docker engine are present. 
  
* When the cluster is ready, you can SSH in the master node and launch the pipeline execution 
  as usual with the `nextflow run <pipeline name>` command.
  
* For the sake of this demo we are using [paraMSA](https://github.com/pditommaso/paraMSA), 
  a pipeline for generating multiple sequence alignments and bootstrap replicates developed 
  in our lab. 

* Nextflow automatically pulls the pipeline code from its GitHub repository when the 
  execution is launched. This repository includes also a dataset which is used by default.  
  [The many bioinformatic tools used by the pipeline](https://github.com/pditommaso/paraMSA#dependencies-) 
  are packaged using a Docker image, which is downloaded automatically on each computing node. 
  
* The pipeline results are uploaded automatically in the S3 bucket specified 
  by the `--output s3://cbcrg-eu/para-msa-results` command line option. 
  
* When the computation is completed, the cluster can be safely shutdown and the 
  EC2 instances terminated with the `nextflow cloud shutdown` command. 
  
### Try it yourself 

We are releasing the Nextflow integrated cloud support in the upcoming version `0.22.0`.  You 
can test it now by defining the following environment variable and running the Nextflow 
installer script as shown below: 

		export NXF_VER=0.22.0-RC1 
		curl get.nextflow.io | bash 
		
Bare in mind that Nextflow requires a Unix-like operating system and a Java runtime version 7+ 
(Windows 10 users which have installed the [Ubuntu subsystem](https://blogs.windows.com/buildingapps/2016/03/30/run-bash-on-ubuntu-on-windows/)
should be able to run it, at their risk..). 

Once you have installed it, you can follow the steps in the above demo. For your convenience 
we made publicly available the EC2 image `ami-43f49030` (EU Ireland region) used to record this 
screencast. 
		
Also make sure you have the following the following variables defined in your environment:

    AWS_ACCESS_KEY_ID="<your aws access key>"
    AWS_SECRET_ACCESS_KEY="<your aws secret key>"
    AWS_DEFAULT_REGION="<your aws region>"
    		
        
### Conclusion 

Nextflow provides state of the art support for cloud and containers technologies making 
it possibile to create computing clusters in the cloud and deploy computational workflows 
in a no-brainer way, with just two commands on your terminal. 

In an upcoming post I will describe the autoscaling capabilities implemented by the 
Nextflow scheduler that allows, along with the use of spot/preemptible instances, 
a cost effective solution for the execution of your pipeline in the cloud. 

#### Credits 

Thanks to [Evan Floden](https://github.com/skptic) for reviewing this post and for writing 
the [paraMSA](https://github.com/skptic/paraMSA/) pipeline.

  
 


