title=The State of Nextflow and Kubernetes
date=2023-03-10
type=post
description=It seems that Kubernetes is everywhere we look! Nextflow users increasingly see K8s as a viable compute environment for their bioinformatic pipelines. In this article, Ben Sherman, Seqera Lab’s master of all things Kubernetes, details recent improvements to the Nextflow/Kubernetes integration.
image=img/the-state-of-nextflow-and-kubernetes.jpg
tags=nextflow, kubernetes
status=published
author=Ben Sherman
icon=ben.jpg
~~~~~~

Hi, my name is Ben, and I’m a software engineer at Seqera Labs. I joined Seqera in November 2021 after finishing my Ph.D. at Clemson University. I work on a number of things at Seqera, but my primary role is that of a Nextflow core contributor.

I have run Nextflow just about everywhere, from my laptop to my university cluster to the cloud and Kubernetes. I have written Nextlfow pipelines for bioinformatics and machine learning, and I even wrote a pipeline to run other Nextflow pipelines for my [dissertation research](https://github.com/bentsherman/tesseract). While I tried to avoid contributing code to Nextflow as a student (I had enough work already), now I get to work on it full-time!

Which brings me to the topic of this post: Nextflow and Kubernetes.

One of my first contributions was a “[best practices guide](https://github.com/seqeralabs/nf-k8s-best-practices)” for running Nextflow on Kubernetes. The guide has helped many people, but for me it provided a map for how to improve K8s support in Nextflow. You see, Nextflow was originally built for HPC, while Kubernetes and cloud batch executors were added later. While Nextflow’s extensible design makes adding features like new executors relatively easy, support for Kubernetes is still a bit spotty.

So, I set out to make Nextflow + K8s great! Over the past year, in collaboration with talented members of the Nextflow community, we have added all sorts of enhancements to the K8s executor. In this blog post, I’d like to show off all of these improvements in one place. So here we go!

## New features

### Submit tasks as Kubernetes Jobs

*New in version 22.05.0-edge.*

Nextflow submits tasks as Pods by default, which is sort of a bad practice. In Kubernetes, every Pod should be created through a controller (e.g., Deployment, Job, StatefulSet) so that Pod failures can be handled automatically. For Nextflow, the appropriate controller is a K8s Job. Using Jobs instead of Pods directly has greatly improved the stability of large Nextflow runs on Kubernetes, and will likely become the default behavior in a future version.

You can enable this feature with the following configuration option:

```groovy
k8s.computeResourceType = ‘Job’
```

Credit goes to @xhejtman from CERIT-SC for leading the charge on this one!

### Object storage as the work directory

*New in version 22.10.0.*

One of the most difficult aspects of using Nextflow with Kubernetes is that Nextflow needs a `PersistentVolumeClaim` (PVC) to store the shared work directory, which also means that Nextflow itself must run inside the Kubernetes cluster in order to access this storage. While the `kuberun` command attempts to automate this process, it has never been reliable enough for production usage.

At the Nextflow Summit in October 2022, we introduced [Fusion](https://seqera.io/fusion/), a file system driver that can mount S3 buckets as POSIX-like directories. The combination of Fusion and [Wave](https://seqera.io/wave/) (a just-in-time container provisioning service) enables you to have your work directory in S3-compatible storage. See the [Wave blog post](https://nextflow.io/blog/2022/rethinking-containers-for-cloud-native-pipelines.html) for an explanation of how it works — it’s pretty cool.

This functionality is useful in general, but it is especially useful for Kubernetes, because (1) you don’t need to provision your own PVC and (2) you can run Nextflow on Kubernetes without using `kuberun` or creating your own submitter Pod.

This feature currently supports AWS S3 on Elastic Kubernetes Service (EKS) clusters and Google Cloud Storage on Google Kubernetes Engine (GKE) clusters.

Check out [this article](https://seqera.io/blog/deploying-nextflow-on-amazon-eks/) over at the Seqera blog for an in-depth guide to running Nextflow (with Fusion) on Amazon EKS.

### No CPU limits by default

*New in version 22.11.0-edge.*

We have changed the default behavior of CPU requests for the K8s executor. Before, a single number in a Nextflow resource request (e.g., `cpus = 8`) was interpreted as both a “request” (lower bound) and a “limit” (upper bound) in the Pod definition. However, setting an explicit CPU limit in K8s is increasingly seen as an anti-pattern (see [this blog post](https://home.robusta.dev/blog/stop-using-cpu-limits) for an explanation). The bottom line is that it is better to specify a request without a limit, because that will ensure that each task has the CPU time it requested, while also allowing the task to use more CPU time if it is available. Unlike other resources like memory and disk, CPU time is compressible — it can be given and taken away without killing the application.

We have also updated the Docker integration in Nextflow to use [CPU shares](https://www.batey.info/cgroup-cpu-shares-for-docker.html), which is the mechanism used by [Kubernetes](https://www.batey.info/cgroup-cpu-shares-for-kubernetes.html) and [AWS Batch](https://docs.aws.amazon.com/AmazonECS/latest/developerguide/task_definition_parameters.html#container_definitions) under the hood to define expandable CPU requests. These changes make the behavior of CPU requests in Nextflow much more consistent across executors.

### CSI ephemeral volumes

*New in version 22.11.0-edge.*

In Kubernetes, volumes are used to provide storage and data (e.g., configuration and secrets) to Pods. Persistent volumes exist independently of Pods and can be mounted and unmounted over time, while ephemeral volumes are attached to a single Pod and are created and destroyed alongside it. While Nextflow can use any persistent volume through a `PersistentVolumeClaim`, ephemeral volume types are supported on a case-by-case basis. For example, `ConfigMaps` and `Secrets` are two ephemeral volume types that are already supported by Nextflow.

Nextflow now also supports [CSI ephemeral volumes](https://kubernetes.io/docs/concepts/storage/ephemeral-volumes/#csi-ephemeral-volumes). CSI stands for Container Storage Interface, and it is a standard used by Kubernetes to support third-party storage systems as volumes. The most common example of a CSI ephemeral volume is [Secrets Store](https://secrets-store-csi-driver.sigs.k8s.io/getting-started/usage.html), which is used to inject secrets from a remote vault such as [Hashicorp Vault](https://www.vaultproject.io/) or [Azure Key Vault](https://azure.microsoft.com/en-us/products/key-vault/).

*Note: CSI persistent volumes can already be used in Nextflow through a `PersistentVolumeClaim`.*

### Local disk storage for tasks

*New in version 22.11.0-edge.*

Nextflow uses a shared work directory to coordinate tasks. Each task receives its own subdirectory with the required input files, and each task is expected to write its output files to this directory. As a workflow scales to thousands of concurrent tasks, this shared storage becomes a major performance bottleneck. We are investigating a few different ways to overcome this challenge.

One of the tools we have to reduce I/O pressure on the shared work directory is to make tasks use local storage. For example, if a task takes input file A, produces an intermediate file B, then produces an output file C, the file B can be written to local storage instead of shared storage because it isn’t a required output file. Or, if the task writes an output file line by line instead of all at once at the end, it can stream the output to local storage first and then copy the file to shared storage.

While it is far from a comprehensive solution, local storage can reduce I/O congestion in some cases. Provisioning local storage for every task looks different on every platform, and in some cases it is not supported. Fortunately, Kubernetes provides a seamless interface for local storage, and now Nextflow supports it as well.

To provision local storage for tasks, you must (1) add an `emptyDir` volume to your Pod options, (2) request disk storage via the `disk` directive, and (3) direct tasks to use the local storage with the `scratch` directive. Here’s an example:

```groovy
process {
    disk = 10.GB
    pod = [ [emptyDir: [:], mountPath: '/scratch'] ]
    scratch = '/scratch'
}
```

As a bonus, you can also provision an `emptyDir` backed by memory:

```groovy
process {
    memory = 10.GB
    pod = [ [emptyDir: [medium: 'Memory'], mountPath: '/scratch'] ]
    scratch = '/scratch'
}
```

Nextflow maps the `disk` directive to the `[ephemeral-storage](https://kubernetes.io/docs/concepts/configuration/manage-resources-containers/#setting-requests-and-limits-for-local-ephemeral-storage)` resource request, which is provided by the `[emptyDir](https://kubernetes.io/docs/concepts/storage/volumes/#emptydir)` volume (another ephemeral volume type).

### Miscellaneous

Check the [release notes](https://github.com/nextflow-io/nextflow/releases) or the list of [K8s pull requests](https://github.com/nextflow-io/nextflow/pulls?q=is%3Apr+label%3Aplatform%2Fk8s) on Github to see what else has been added. Here are some notable improvements from the past year:

- Support Pod `affinity` ([640cbed4](https://github.com/nextflow-io/nextflow/commit/640cbed4813a34887d4dc10f87fa2e4aa524d055))
- Support Pod `automountServiceAccountToken` ([1b5908e4](https://github.com/nextflow-io/nextflow/commit/1b5908e4cbbb79f93be2889eec3acfa6242068a1))
- Support Pod `priorityClassName` ([51650f8c](https://github.com/nextflow-io/nextflow/commit/51650f8c411ba40f0966031035e7a47c036f542e))
- Support Pod `tolerations` ([7f7cdadc](https://github.com/nextflow-io/nextflow/commit/7f7cdadc6a36d0fb99ef125f6c6f89bfca8ca52e))
- Support `time` directive via `activeDeadlineSeconds` ([2b6f70a8](https://github.com/nextflow-io/nextflow/commit/2b6f70a8fa55b993fa48755f7a47ac9e1b584e48))
- Improved control over error conditions ([064f9bc4](https://github.com/nextflow-io/nextflow/commit/064f9bc4), [58be2128](https://github.com/nextflow-io/nextflow/commit/58be2128), [d86ddc36](https://github.com/nextflow-io/nextflow/commit/d86ddc36))
- Improved support for labels and queue annotation ([9951fcd9](https://github.com/nextflow-io/nextflow/commit/9951fcd9), [4df8c8d2](https://github.com/nextflow-io/nextflow/commit/4df8c8d2))
- Add support for AWS IAM role for Service Accounts ([62df42c3](https://github.com/nextflow-io/nextflow/commit/62df42c3), [c3364d0f](https://github.com/nextflow-io/nextflow/commit/c3364d0f), [b3d33e3b](https://github.com/nextflow-io/nextflow/commit/b3d33e3b))

## Beyond Kubernetes

We’ve added tons of value to Nextflow over the past year – not just in terms of Kubernetes support, but also in terms of performance, stability, and integrations with other technologies – and we aren’t stopping any time soon! We have greater ambitions still for Nextflow, and I for one am looking forward to what we will accomplish together. As always, keep an eye on this blog, as well as the [Nextflow GitHub](https://github.com/nextflow-io/nextflow) page, for the latest updates to Nextflow.
