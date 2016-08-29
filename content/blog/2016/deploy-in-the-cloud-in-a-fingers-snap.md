title=Deploy your computational pipelines in the cloud at the snap-of-a-finger 
date=2016-08-30
type=post
tags=aws,cloud,pipelines,nextflow,genomic,docker
status=published
author=Paolo Di Tommaso 
icon=paolo.jpg
~~~~~~

*Learn how deploy and run a computational pipeline in the Amazon AWS cloud with ease 
thanks to Nextflow and Docker containers*

Nextflow is a framework to simplify the writing of parallel and distributed computational
pipelines in a portable and reproducible way across different computing platforms, from 
a single computer to the cloud, with minimal maintenance or interventions. 

Indeed, the very original idea, when this project started three years ago, was to 
implement a tool that would allow the researchers in 
[our lab](http://www.crg.eu/es/programmes-groups/comparative-bioinformatics) to smoothly migrate 
their data analysis applications in the cloud when needed - without having 
to change or adapt their code. 

However to date Nextflow has been used mostly to deploy computational workflows at on-premise 
computing clusters or HPC data-centers, such as 
[BSC MareNostrum](https://www.bsc.es/marenostrum-support-services/mn3), because these 
infrastructures are easier to use and provide, on average, a cheaper cost and better 
performance when compared to a cloud environments. 

A major obstacle to efficient deployment of scientific workflows in the cloud is the lack 
of this kind of infrastructure for a performant POSIX compatible shared file system. These
kinds of applications are usually made-up by putting together a collection of tools, scripts and 
systems command that need a reliable shared file system to share with each other the input and 
output of files as they are produced. 

The recent availability of the [Amazon Elastic File System](https://aws.amazon.com/efs/) 
(EFS), a fully featured NFS based file system hosted on the AWS infrastructure represents 
a major step in this context, unlocking the deployment of scientific computing 
in the cloud and taking it to the next level. 

### Nextflow support for the cloud 

Nextflow could already be deployed in the cloud, either using tools such as 
[ElastiCluster](https://github.com/gc3-uzh-ch/elasticluster) or [CfnCluster](https://aws.amazon.com/hpc/cfncluster/),
or by using custom deployment scripts. However the procedure was still cumbersome and,
above all, it was not optimised to fully take advantage of cloud elasticity i.e. 
the ability to shape the computing cluster dynamically as the computing needs change.

For these reasons, we decided it was time to provide in Nextflow a first-class support 
for the cloud, implementing a native support for Amazon EFS, a revamped cloud-optimised 
scheduler based on [Apache Ignite](https://ignite.apache.org/) and full-support for cluster 
auto-scaling. 

In practice this means that Nextflow can now spin and configure a fullly featured computing 
cluster in the cloud with a single command, after that you need only to login in the master 
node and launch the pipeline execution as you would do in your cluster. 


### Demo !

Since a demo is worth a thousands words, I've record a short screencast showing to use
Nextflow to setup a cluster in the cloud made-up of the 10 EC2 spot instances and using 
the Amazon EFS file system. When the cluster is ready, the only requirement is to SSH login 
to the cluster master node and launch the Nextflow pipeline execution as usual. 

<script type="text/javascript" src="https://asciinema.org/a/bq44gqp2sizxen4i05a6vomd6.js" id="asciicast-bq44gqp2sizxen4i05a6vomd6" async></script> 
   
Let's recap the steps showed in the demo 
 
* The user provides the cloud parameters (such as the VM image ID and the instance type) 
  in the `nextflow.config` file. 

* To configure the EFS file sytem you need to provide your storage ID and the mount path
  by using the  `sharedStorageId` and `sharedStorageMount` properties.  
  
* To use [EC2 Spot](https://aws.amazon.com/ec2/spot/) instances, just specify the price 
  you want to bid by using the `spotPrice` property.
  
* The AWS access and secret keys are provided by using the usual environment variables. 

* The `cloud launch` launches the requested number of instances, configures the user and 
  access keys, mounts the EFS storage and setsup the Nextflow cluster automatically.
  Any Linux AMI can be used, it only required that the [cloud-init](https://cloudinit.readthedocs.io/en/latest/) 
  package, a Java 7+ runtime and the Docker engine are present. 
  
* When the cluster is ready, you can SSH in the master node and launch the pipeline execution 
  as usual with the `nextflow run ..` command.
  
* For the sake of this demo we are using [pditommaso/paraMSA](https://github.com/pditommaso/paraMSA), 
  a pipeline for generating multiple sequence alignments and bootstrap replicates 
  developed by [Evan Floden](https://github.com/skptic). 

* Nextflow automatically puuls the pipeline code from its GitHub repository when the 
  execution is launched. This repository includes also a dataset which is used by default.  
  The many bioinformatic tools used by the pipeline are packaged using a Docker image, 
  which is downloaded automatically on each computing node. 
  
* The pipeline results are uploaded automatically in the S3 bucket specified 
  by the `--output s3://cbcrg-eu/para-msa-results` command line option. 
  
* When the computation is completed we can safely shutdown then cluster and terminate the 
  EC2 instances with the `nextflow cloud shutdown command` 
  
### Try it yourself 

We are releasing the Nextflow integrated cloud support and technological preview, and 
you can access it by using the Nextflow development development snapshot as shown below: 

		export NXF_VER=0.23.0-SNAPSHOT 
		curl get.nextflow.io | bash 
		
Bare in mind that Nextflow requires a Unix-like operating system and a Java runtime version 7 (or higher).

Once you have installed it, you can follow the steps in the above demo. For your convenience 
we made public available the EC2 image `ami-43f49030` (EU Ireland region) used in this demo. 
		
        
### Conclusion 

Nextflow provides state-of-art support for cloud and containers technologies making 
it possibile to create and configure computing clusters in the cloud. The deployment and execution 
of computational workflows become executible in a no-brainer way, with just two commands on
your terminal. 

In an upcoming post I will describe the autoscaling capabilities implemented by the 
Nextflow cloud native scheduler that allows, along with the use of 
spot/preemptible instances, a cost effective solution to execution of your pipeline in the cloud. 

  
