title=:birthday: Nextflow Turns 5! :birthday:
date=2019-03-29
type=post
tags=nextflow,kubernetes
status=published
author=Evan Floden
icon=evan.jpg
~~~~~~

Nextflow is growing up. Last week marked five years since the [first commit](https://github.com/nextflow-io/nextflow/commit/c080150321e5000a2c891e477bb582df07b7f75f) of the project which took place on the 22 of March 2013. Like a parent reflecting on their child attending school for the first time, we know reaching this point hasn’t been an entirely solo journey, despite Paolo's best efforts!

A lot has happened recently and we thought it was time to highlight some of the recent evolutions. We also take the opportunity to extend the warmest of thanks to all those who have contributed to the development of Nextflow as well as the fantastic community of users who consistently provide ideas, feedback and the occasional late night banter on the [Gitter channel](https://gitter.im/nextflow-io/nextflow).

 Here are a few neat developments churning out of the birthday cake mixer.

### nf-core
[nf-core](https://nf-core.github.io/) is a community effort to provide a home for high quality, production-ready, curated analysis pipelines built using Nextflow. The project has been initiated and is being led by [Phil Ewels](https://github.com/ewels) of [MultiQC](http://multiqc.info/) fame. The principle is that nf-core pipelines can be used out-of-the-box or as inspiration for something different.

As well as being a place for best-practise pipelines, other features of nf-core include the [cookie cutter template tool](https://github.com/nf-core/cookiecutter) which provides a fast way to create a dependable workflow using many of Nextflow’s sweet capabilities such as:

* *Outline:* Skeleton pipeline script
* *Data:* Reference Genome implementation (AWS iGenomes)
* *Configuration:* Robust configuration setup
* *Containers:* Skeleton files for Docker image generation
* *Reporting:* HTML email functionality and and HTML results output
* *Documentation:* Installation, Usage, Output, Troubleshooting etc
* *Continuous Integration:* Skeleton files for automated testing using Travis CI

There is also a Python package with helper tools for Nextflow.

You can find more information about the community via the projects [website](nf-core.github.io), [GitHub repository](https://github.com/nf-core), [Twitter account](https://twitter.com/nf_core) or join the dedicated [Gitter](https://gitter.im/nf-core/Lobby) chat.

<img alt='nf-core logo' width='760' src='/img/nf-core-logo.png' style='margin:1em auto'/>

### Kubernetes has landed

As of v0.28.0 Nextflow has now has support for kubebernetes. If you don’t know much about kubernetes, at its heart it is an open-source platform for the management and deployment of containers. Google led the initial design and it is now maintained by the Cloud Native Computing Foundation. I found the [The Illustrated Children's Guide to Kubernetes](https://www.youtube.com/watch?v=4ht22ReBjno) particularly useful in explaining the basic vocabulary and concepts.

Kubernetes looks be one of the key technologies for the application of containers in the cloud as well as for building Infrastructure as a Service (IaaS) and Platform and a Service (PaaS) applications. We have been approached by many users who wish to use Nextflow with kubernetes to be able to deploy workflows across both academic and commercial settings. With enterprise versions of kubernetes such as Red Hat's [OpenShift](https://www.openshift.com/), it was becoming apparent there was a need for native execution with Nextflow.

The new command `nextflow kuberun` launches the Nextflow driver as a “pod” which is then able to launch workflow tasks as other pods within a kubernetes cluster. You can read more in the documentation on kubernetes support for Nextflow [here](https://www.nextflow.io/docs/latest/kubernetes.html) and follow a basic walk-through example on how to get up and running at the [bottom](### Get started with kubernetes and nextflow) of this post.

<img alt='Nextflow and kubernetes' width='760' src='/img/nextflow-kubernetes.png' style='margin:1em auto'/>

### Improved reporting and notifications

Following the hackathon in September we wrote about of the addition of HTML trace reports  that allow for the generation HTML detailing resource usage (CPU time, memory, disk i/o etc).

Thanks to valuable feedback there has continued to be many improvements to the reportsm as tracked through the Nextflow GitHub issues page. Reports are now able to display [thousands of tasks](https://github.com/nextflow-io/nextflow/issues/547) and include extra information such as the [container engine used](https://github.com/nextflow-io/nextflow/issues/521). Tasks can be filtered and an [overall progress bar](https://github.com/nextflow-io/nextflow/issues/534) has been added.

You can explore a [real-world HTML report](/misc/nf-trace-report2.html) and more information on HTML reports can be found in the [docs](https://www.nextflow.io/docs/latest/tracing.html).

There has also been additions to workflow notifications. Currently these can be configured to automatically send a notification message when a workflow execution terminates. You can read more about how to setup notifications in the [docs](https://www.nextflow.io/docs/latest/mail.html?highlight=notification#workflow-notification).

### Syntax-tic!

Writing workflows no longer has to be done in monochrome. There is now syntax highlighting for Nextflow in the popular Atom Text Editor as well as in Visual Studio.

You can find the Atom plugin by searching for Nextflow in Atoms package installer or clicking [here](https://atom.io/packages/language-nextflow). The Visual Studio plugin can be downloaded [here](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow)

On a related note, Nextflow is now an official language on GitHub!

<img alt='Nextflow syntax highlighting with atom' width='760' src='/img/atom-syntax.png' style='margin:1em auto'/>

### Get started with kubernetes and nextflow

The following walk-through takes advantage of the recent addition of [kubernetes to docker](https://www.docker.com/kubernetes). As of writing, kubernetes support in Docker for Mac part of the "Edge" release but is planned to be incorporated in the main releases soon.

The following code runs is confirmed to run on MacOS (High Sierra) running Docker version 18.03.0-ce-rc4-mac57.

1. Open Docker and ensure it is running with kubernetes
<img alt='Docker for mac with kubernetes' width='760' src='/img/docker-for-mac-with-kubernetes.png' style='margin:1em auto'/>

2. Install Nextflow if necessary
```
curl -s https://get.nextflow.io | bash
```

3. Create a kubernetes persistent volume
```
kubectl create -f https://k8s.io/docs/tasks/configure-pod-container/task-pv-volume.yaml
```

4. Create a kubernetes persistent volume claim
```
kubectl create -f https://k8s.io/docs/tasks/configure-pod-container/task-pv-claim.yaml
```

5. Run any Nextflow pipeline on your new kubernetes clusters
```
nextflow kuberun nextflow-io/rnatoy -v task-pv-claim:/mnt/
```
