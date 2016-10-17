title=Enabling effective cloud computing 
date=2016-10-18
type=post
tags=aws,cloud,pipelines,nextflow,genomic,docker
status=published
author=Paolo Di Tommaso 
icon=paolo.jpg
~~~~~~



In the [previous post](/blog/2016/deploy-in-the-cloud-at-snap-of-a-finger.html) I've introduced the new cloud native support for AWS provided by Nextflow. 

This allows the creation of a computing in the cloud in a no-brainer way, enabling the 
deployment of complex computational pipelines in a few commands. 

This solution is characterised used a super-thin application stack which does not 
require any third party component in the cloud instances other than the Java VM and the 
Docker engine. 

![Nextflow cloud deployment](/img/cloud-deployment.png)


Each EC2 instance runs during he bootstrap a script that mounts the [EFS](https://aws.amazon.com/efs/) 
storage and downloads and launches the Nextflow cluster daemon. This daemon is self-configuring, 
it automatically discovers the other running instances and join them forming the computing cluster. 

### Going elastic 

The simplicity of this stack makes it possible to setup the cluster in the cloud a few minutes, 
little more the time required to spin up the EC2 instances, and above all independently the number 
launched instances because each of them configures itself in an autonomous manner. 

This makes also possible to add or remove instances as needed, realising the [long promised 
elastic scalability](http://www.nextplatform.com/2016/09/21/three-great-lies-cloud-computing/) 
of cloud computing.  

This ability is even more important for bioinformatic workloads which frequently are composed 
by tasks with very different computing requirements (eg. a few very long running tasks and 
many short-lived tasks) and crunch not homogenous dataset.

The Nextflow cloud scheduler features an elastic cluster that is able to resize itself 
to adapt to the actual computing needs at runtime, thus spinning up new EC2 instances 
when there are jobs waiting too long in the execution queue or terminating instances not 
used for a certain amount of time. 

In order to enable the cluster autoscaling you will need to specify the autoscale 
properties in the `nextflow.config` file. For example: 

```
cloud {
  imageId = 'ami-43f49030'
  instanceType = 'm4.xlarge'

  autoscale {
   	 enabled = true
     minInstances = 5
     maxInstances = 10 
  }
}
``` 

The above configuration enables the autoscaling features in such a way that the cluster 
will be made up of at least 5 nodes. If at any point some tasks will spend more than 5 minutes 
without being processed, a number of instances needed to fullfil the pending tasks - up to 
the specified limit ie. 10 - will be launched. On the other hand, if these instances are 
idle, they will be  terminated before reaching the 60 minutes usage boundary. The autoscaler 
launches instances by using the same AMI ID and type specified in the `cloud` configuration. 

However it is possible to define different instance attributes as shown below: 

```
cloud {
  imageId = 'ami-43f49030'
  instanceType = 'm4.large'

  autoscale {
   	 enabled = true
     maxInstances = 10 
  	 instanceType = 'm4.2xlarge'
  	 spotPrice = 0.05 
  }
}
``` 

The cluster is first created by using instance(s) of type `m4.large`. Then, when new 
computing nodes are required the autoscaler will launch instances of type `m4.2xlarge`. 
Also, since the `spotPrice` attribute is specified, [EC2 spot](https://aws.amazon.com/ec2/spot/) 
instances are launched, instead of regular on-demand ones, bidding for the price specified. 

### Conclusion 


xxx


