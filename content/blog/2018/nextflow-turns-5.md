title=Nextflow turns five! Happy birthday!
date=2019-04-3
type=post
tags=nextflow,kubernetes,nf-core
status=published
author=Evan Floden
icon=evan.jpg
~~~~~~

Nextflow is growing up. The past week marked five years since the [first commit](https://github.com/nextflow-io/nextflow/commit/c080150321e5000a2c891e477bb582df07b7f75f) of the project which took place on the 22 of March 2013. Like a parent reflecting on their child attending school for the first time, we know reaching this point hasn’t been an entirely solo journey, despite Paolo's best efforts!

A lot has happened recently and we thought it was time to highlight some of the recent evolutions. We also take the opportunity to extend the warmest of thanks to all those who have contributed to the development of Nextflow as well as the fantastic community of users who consistently provide ideas, feedback and the occasional late night banter on the [Gitter channel](https://gitter.im/nextflow-io/nextflow).

 Here are a few neat developments churning out of the birthday cake mix.

### NF-Core

[NF-Core](https://nf-core.github.io/) is a community effort to provide a home for high quality, production-ready, curated analysis pipelines built using Nextflow. The project has been initiated and is being led by [Phil Ewels](https://github.com/ewels) of [MultiQC](http://multiqc.info/) fame. The principle is that NF-Core pipelines can be used out-of-the-box or as inspiration for something different.

As well as being a place for best-practise pipelines, other features of NF-Core include the [cookie cutter template tool](https://github.com/nf-core/cookiecutter) which provides a fast way to create a dependable workflow using many of Nextflow’s sweet capabilities such as:

* *Outline:* Skeleton pipeline script.
* *Data:* Reference Genome implementation (AWS iGenomes).
* *Configuration:* Robust configuration setup.
* *Containers:* Skeleton files for Docker image generation.
* *Reporting:* HTML email functionality and and HTML results output.
* *Documentation:* Installation, Usage, Output, Troubleshooting, etc.
* *Continuous Integration:* Skeleton files for automated testing using Travis CI.

There is also a Python package with helper tools for Nextflow.

You can find more information about the community via the project [website](nf-core.github.io), [GitHub repository](https://github.com/nf-core), [Twitter account](https://twitter.com/nf_core) or join the dedicated [Gitter](https://gitter.im/nf-core/Lobby) chat.

<img alt='nf-core logo' width='560' src='/img/nf-core-logo.png' style='margin:1em auto'/>

### Kubernetes has landed

As of version 0.28.0 Nextflow now has support for Kubernetes. If you don’t know much about Kubernetes, at its heart it is an open-source platform for the management and deployment of containers at scale. Google led the initial design and it is now maintained by the Cloud Native Computing Foundation. I found the [The Illustrated Children's Guide to Kubernetes](https://www.youtube.com/watch?v=4ht22ReBjno) particularly useful in explaining the basic vocabulary and concepts.

Kubernetes looks be one of the key technologies for the application of containers in the cloud as well as for building Infrastructure as a Service (IaaS) and Platform and a Service (PaaS) applications. We have been approached by many users who wish to use Nextflow with Kubernetes to be able to deploy workflows across both academic and commercial settings. With enterprise versions of Kubernetes such as Red Hat's [OpenShift](https://www.openshift.com/), it was becoming apparent there was a need for native execution with Nextflow.

The new command `nextflow kuberun` launches the Nextflow driver as a *pod* which is then able to launch workflow tasks as other pods within a Kubernetes cluster. You can read more in the documentation on Kubernetes support for Nextflow [here](https://www.nextflow.io/docs/latest/kubernetes.html).

<img alt='Nextflow and Kubernetes' width='760' src='/img/nextflow-kubernetes.png' style='margin:1em auto'/>

### Improved reporting and notifications

Following the hackathon in September we wrote about the addition of HTML trace reports that allow for the generation HTML detailing resource usage (CPU time, memory, disk i/o etc).

Thanks to valuable feedback there has continued to be many improvements to the reports as tracked through the Nextflow GitHub issues page. Reports are now able to display [thousands of tasks](https://github.com/nextflow-io/nextflow/issues/547) and include extra information such as the [container engine used](https://github.com/nextflow-io/nextflow/issues/521). Tasks can be filtered and an [overall progress bar](https://github.com/nextflow-io/nextflow/issues/534) has been added.

You can explore a [real-world HTML report](/misc/nf-trace-report2.html) and more information on HTML reports can be found in the [documentation](https://www.nextflow.io/docs/latest/tracing.html).

There has also been additions to workflow notifications. Currently these can be configured to automatically send a notification email when a workflow execution terminates. You can read more about how to setup notifications in the [documentation](https://www.nextflow.io/docs/latest/mail.html?highlight=notification#workflow-notification).

### Syntax-tic!

Writing workflows no longer has to be done in monochrome. There is now syntax highlighting for Nextflow in the popular [Atom editor](https://atom.io) as well as in [Visual Studio Code](https://code.visualstudio.com).

<img alt='Nextflow syntax highlighting with Atom' width='360' src='/img/atom-min.png' style='margin:1em auto; cursor: pointer;' onclick="window.open(this.src)" />
<img alt='Nextflow syntax highlighting with VSCode' width='360' src='/img/vscode-min.png' style='margin:1em auto; cursor: pointer;' onclick="window.open(this.src)" />

You can find the Atom plugin by searching for Nextflow in Atoms package installer or clicking [here](https://atom.io/packages/language-nextflow). The Visual Studio plugin can be downloaded [here](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow).

On a related note, Nextflow is now an [official language on GitHub](https://github.com/search?q=language%3Anextflow)!

<img alt='GitHub nextflow syntax' width='760' src='/img/github-nf-syntax-min.png' style='margin:1em auto'/>

### Conclusion

Nextflow developments are progressing faster than ever and with the help of the community, there are a ton of great new features on the way. If you have any suggestions of your killer NF idea then please drop us a line, open an issue or even better, join in on the fun.

Over the coming months Nextflow will be reaching out with several training sessions across the US and Europe. We hope to see as many of you as possible on the road.
