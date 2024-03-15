---
title: Scaling with AWS Batch
date: 2017-11-08
type: post
tags: pipelines,nextflow,genomic,workflow,aws,batch
author: Paolo Di Tommaso
icon: paolo.jpg
---

The latest Nextflow release (0.26.0) includes built-in support for [AWS Batch](https://aws.amazon.com/batch/),
a managed computing service that allows the execution of containerised workloads
over the Amazon EC2 Container Service (ECS).

This feature allows the seamless deployment of Nextflow pipelines in the cloud by offloading
the process executions as managed Batch jobs. The service takes care to spin up the required
computing instances on-demand, scaling up and down the number and composition of the instances
to best accommodate the actual workload resource needs at any point in time.

AWS Batch shares with Nextflow the same vision regarding workflow containerisation
i.e. each compute task is executed in its own Docker container. This dramatically
simplifies the workflow deployment through the download of a few container images.
This common design background made the support for AWS Batch a natural extension for Nextflow.

### Batch in a nutshell

Batch is organised in _Compute Environments_, _Job queues_, _Job definitions_ and _Jobs_.

The _Compute Environment_ allows you to define the computing resources required for a specific workload (type).
You can specify the minimum and maximum number of CPUs that can be allocated,
the EC2 provisioning model (On-demand or Spot), the AMI to be used and the allowed instance types.

The _Job queue_ definition allows you to bind a specific task to one or more Compute Environments.

Then, the _Job definition_ is a template for one or more jobs in your workload. This is required
to specify the Docker image to be used in running a particular task along with other requirements
such as the container mount points, the number of CPUs, the amount of memory and the number of
retries in case of job failure.

Finally the _Job_ binds a Job definition to a specific Job queue
and allows you to specify the actual task command to be executed in the container.

The job input and output data management is delegated to the user. This means that if you
only use Batch API/tools you will need to take care to stage the input data from a S3 bucket
(or a different source) and upload the results to a persistent storage location.

This could turn out to be cumbersome in complex workflows with a large number of
tasks and above all it makes it difficult to deploy the same applications across different
infrastructure.

### How to use Batch with Nextflow

Nextflow streamlines the use of AWS Batch by smoothly integrating it in its workflow processing
model and enabling transparent interoperability with other systems.

To run Nextflow you will need to set-up in your AWS Batch account a [Compute Environment](http://docs.aws.amazon.com/batch/latest/userguide/compute_environments.html)
defining the required computing resources and associate it to a [Job Queue](http://docs.aws.amazon.com/batch/latest/userguide/job_queues.html).

Nextflow takes care to create the required _Job Definitions_ and _Job_ requests as needed.
This spares some Batch configurations steps.

In the `nextflow.config`, file specify the `awsbatch` executor, the Batch `queue` and
the container to be used in the usual manner. You may also need to specify the AWS region
and access credentials if they are not provided by other means. For example:

    process.executor = 'awsbatch'
    process.queue = 'my-batch-queue'
    process.container = your-org/your-docker:image
    aws.region = 'eu-west-1'
    aws.accessKey = 'xxx'
    aws.secretKey = 'yyy'

Each process can eventually use a different queue and Docker image (see Nextflow documentation for details).
The container image(s) must be published in a Docker registry that is accessible from the
instances run by AWS Batch eg. [Docker Hub](https://hub.docker.com/), [Quay](https://quay.io/)
or [ECS Container Registry](https://aws.amazon.com/ecr/).

The Nextflow process can be launched either in a local computer or a EC2 instance.
The latter is suggested for heavy or long running workloads.

Note that input data should be stored in the S3 storage. In the same manner
the pipeline execution must specify a S3 bucket as a working directory by using the `-w` command line option.

A final caveat about custom containers and computing AMI. Nextflow automatically stages input
data and shares tasks intermediate results by using the S3 bucket specified as a work directory.
For this reason it needs to use the `aws` command line tool which must be installed either
in your process container or be present in a custom AMI that can be mounted and accessed
by the Docker containers.

You may also need to create a custom AMI because the default image used by AWS Batch only
provides 22 GB of storage which may not be enough for real world analysis pipelines.

See the documentation to learn [how to create a custom AMI](/docs/latest/awscloud.html#custom-ami)
with larger storage and how to setup the AWS CLI tools.

### An example

In order to validate Nextflow integration with AWS Batch, we used a simple RNA-Seq pipeline.

This pipeline takes as input a metadata file from the Encode project corresponding to a [search
returning all human RNA-seq paired-end datasets](https://www.encodeproject.org/search/?type=Experiment&award.project=ENCODE&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&files.file_type=fastq&files.run_type=paired-ended&replicates.library.nucleic_acid_term_name=RNA&replicates.library.depleted_in_term_name=rRNA)
(the metadata file has been additionally filtered to retain only data having a SRA ID).

The pipeline automatically downloads the FASTQ files for each sample from the EBI ENA database,
it assesses the overall quality of sequencing data using FastQC and then runs [Salmon](https://combine-lab.github.io/salmon/)
to perform the quantification over the human transcript sequences. Finally all the QC and
quantification outputs are summarised using the [MultiQC](http://multiqc.info/) tool.

For the sake of this benchmark we used the first 38 samples out of the full 375 samples dataset.

The pipeline was executed both on AWS Batch cloud and in the CRG internal Univa cluster,
using [Singularity](/blog/2016/more-fun-containers-hpc.html) as containers runtime.

It's worth noting that with the exception of the two configuration changes detailed below,
we used exactly the same pipeline implementation at [this GitHub repository](https://github.com/nextflow-io/rnaseq-encode-nf).

The AWS deploy used the following configuration profile:

    aws.region = 'eu-west-1'
    aws.client.storageEncryption = 'AES256'
    process.queue = 'large'
    executor.name = 'awsbatch'
    executor.awscli = '/home/ec2-user/miniconda/bin/aws'

While for the cluster deployment the following configuration was used:

    executor = 'crg'
    singularity.enabled = true
    process.container = "docker://nextflow/rnaseq-nf"
    process.queue = 'cn-el7'
    process.time = '90 min'
    process.$quant.time = '4.5 h'

### Results

The AWS Batch Compute environment was configured to use a maximum of 132 CPUs as the number of CPUs
that were available in the queue for local cluster deployment.

The two executions ran in roughly the same time: 2 hours and 24 minutes when running in the
CRG cluster and 2 hours and 37 minutes when using AWS Batch.

It must be noted that 14 jobs failed in the Batch deployment, presumably because one or more spot
instances were retired. However Nextflow was able to re-schedule the failed jobs automatically
and the overall pipeline execution completed successfully, also showing the benefits of a truly
fault tolerant environment.

The overall cost for running the pipeline with AWS Batch was **$5.47** ($ 3.28 for EC2 instances,
$1.88 for EBS volume and $0.31 for S3 storage). This means that with ~ $55 we could have
performed the same analysis on the full Encode dataset.

It is more difficult to estimate the cost when using the internal cluster, because we don't
have access to such detailed cost accounting. However, as a user, we can estimate it roughly
comes out at $0.01 per CPU-Hour. The pipeline needed around 147 CPU-Hour to carry out the analysis,
hence with an estimated cost of **$1.47** just for the computation.

The execution report for the Batch execution is available at [this link](https://cdn.rawgit.com/nextflow-io/rnaseq-encode-nf/db303a81/benchmark/aws-batch/report.html)
and the one for cluster is available [here](https://cdn.rawgit.com/nextflow-io/rnaseq-encode-nf/db303a81/benchmark/crg-cluster/report.html).

### Conclusion

This post shows how Nextflow integrates smoothly with AWS Batch and how it can be used to
deploy and execute real world genomics pipeline in the cloud with ease.

The auto-scaling ability provided by AWS Batch along with the use of spot instances make
the use of the cloud even more cost effective. Running on a local cluster may still be cheaper,
even if it is non trivial to account for all the real costs of a HPC infrastructure.
However the cloud allows flexibility and scalability not possible with common on-premises clusters.

We also demonstrate how the same Nextflow pipeline can be _transparently_ deployed in two very
different computing infrastructure, using different containerisation technologies by simply
providing a separate configuration profile.

This approach enables the interoperability across different deployment sites, reduces
operational and maintenance costs and guarantees consistent results over time.

### Credits

This post is co-authored with [Francesco Strozzi](https://twitter.com/fstrozzi),
who also helped to write the pipeline used for the benchmark in this post and contributed
to and tested the AWS Batch integration. Thanks to [Emilio Palumbo](https://github.com/emi80)
that helped to set-up and configure the AWS Batch environment and [Evan Floden](https://gitter.im/skptic)
for the comments.
