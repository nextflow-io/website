title=Deploy Nextflow Pipelines with Google Cloud Batch! Extended and updated version
date=2023-02-01
type=post
description=We've talked about deploying Nextflow pipelines with Google Cloud Batch in the past. Join us again for an extended and updated version
image=img/gbatch_extended.jpg
tags=nextflow,google,cloud
status=published
author=Marcel Ribeiro-Dantas
icon=marcel.png
~~~~~~

[We have talked about Google Cloud Batch before](https://www.nextflow.io/blog/2022/deploy-nextflow-pipelines-with-google-cloud-batch.html). Not only that, we were proud to announce Nextflow support to Google Cloud Batch right after it was publicly released, back in July 2022. How amazing is that? But we didn't stop there! The [Nextflow official documentation](https://www.nextflow.io/docs/latest/google.html) also provides a lot of useful information on how to use Google Cloud Batch as the compute environment for your Nextflow pipelines. Having said that, feedback from the community is valuable, and we agreed that in addition to the documentation, teaching by example, and in a more informal language, can help many of our users. So, here is a tutorial on how to use the Batch service of the Google Cloud Platform with Nextflow 🥳

### Running an RNAseq pipeline with Google Cloud Batch

Welcome to our RNAseq tutorial using Nextflow and Google Cloud Batch! RNAseq is a powerful technique for studying gene expression and is widely used in a variety of fields, including genomics, transcriptomics, and epigenomics. In this tutorial, we will show you how to use Nextflow, a popular workflow management tool, to perform an RNAseq analysis on Google Cloud Batch, a scalable cloud-based computing platform.

Nextflow allows you to easily develop, execute, and scale complex pipelines on any infrastructure, including the cloud. Google Cloud Batch enables you to run batch workloads on Google Cloud Platform (GCP), with the ability to scale up or down as needed. Together, Nextflow and Google Cloud Batch provide a powerful and flexible solution for RNAseq analysis.

We will walk you through the entire process, from setting up your Google Cloud account and installing Nextflow to running an RNAseq pipeline and interpreting the results. By the end of this tutorial, you will have a solid understanding of how to use Nextflow and Google Cloud Batch for RNAseq analysis. So let's get started!

### Setting up Google Cloud CLI (gcloud)

In this tutorial, you will learn how to use the gcloud command-line interface to interact with the Google Cloud Platform and set up your Google Cloud account for use with Nextflow. If you do not already have gcloud installed, you can follow the instructions [here](https://cloud.google.com/sdk/docs/install) to install it. Once you have gcloud installed, run the command `gcloud init` to initialize the CLI. You will be prompted to choose an existing project to work on or create a new one. For the purpose of this tutorial, we will create a new project. Name your project "my-rnaseq-pipeline". There may be a lot of information displayed on the screen after running this command, but you can ignore it for now.

### Setting up Batch and Storage in Google Cloud Platform

#### Enable Google Batch

According to the [official Google documentation](https://cloud.google.com/batch/docs/get-started), "Batch is a fully managed service that lets you schedule, queue, and execute [batch processing](https://en.wikipedia.org/wiki/Batch_processing) workloads on Compute Engine virtual machine (VM) instances. Batch provisions resources and manages capacity on your behalf, allowing your batch workloads to run at scale."

The first step is to download the `beta` command group. You can do this by executing:

```bash
$ gcloud components install beta
```

Then, enable billing for this project. You will first need to get your account id with

```bash
$ gcloud beta billing accounts list
```

After that, you will see something like the following appear in your window:

```console
ACCOUNT_ID            NAME                OPEN  MASTER_ACCOUNT_ID
XXXXX-YYYYYY-ZZZZZZ  My Billing Account  True
```

If you get the error “Service Usage API has not been used in project 842841895214 before or it is disabled”, simply run the command again and it should work. Then copy the account id, and the project id and paste them into the command below. This will enable billing for your project id.

```bash
$ gcloud beta billing projects link PROJECT-ID --billing-account XXXXXX-YYYYYY-ZZZZZZ
```

Next, you must enable the Batch API, along with the Compute Engine and Cloud Logging APIs. You can do so with the following command:

```bash
$ gcloud services enable batch.googleapis.com compute.googleapis.com logging.googleapis.com
```

You should see a message similar to the one below:

```console
Operation "operations/acf.p2-AAAA-BBBBB-CCCC--DDDD" finished successfully.
```

#### Create a Service Account

In order to access the APIs we enabled, you need to [create a Service Account](https://cloud.google.com/iam/docs/creating-managing-service-accounts#iam-service-accounts-create-gcloud) and set the necessary IAM roles for the project. You can create the Service Account by executing:

```bash
$ gcloud iam service-accounts create rnaseq-pipeline-sa
```

After this, set appropriate roles for the project using the commands below:

```bash
$ gcloud projects add-iam-policy-binding my-rnaseq-pipeline \
--member="serviceAccount:rnaseq-pipeline-sa@my-rnaseq-pipeline.iam.gserviceaccount.com" \
--role="roles/iam.serviceAccountUser"

$ gcloud projects add-iam-policy-binding my-rnaseq-pipeline \
--member="serviceAccount:rnaseq-pipeline-sa@my-rnaseq-pipeline.iam.gserviceaccount.com" \
--role="roles/batch.jobsEditor"

$ gcloud projects add-iam-policy-binding my-rnaseq-pipeline \
--member="serviceAccount:rnaseq-pipeline-sa@my-rnaseq-pipeline.iam.gserviceaccount.com" \
--role="roles/logging.viewer"

$ gcloud projects add-iam-policy-binding my-rnaseq-pipeline \
--member="serviceAccount:rnaseq-pipeline-sa@my-rnaseq-pipeline.iam.gserviceaccount.com" \
--role="roles/storage.admin"
```

#### Create your Bucket

Now it's time to create your Storage bucket, where both your input, intermediate and output files will be hosted and accessed by the Google Batch virtual machines. Your bucket name must be globally unique (across regions). For the example below, the bucket is named rnaseq-pipeline-nextflow-bucket. However, as this name has now been used you have to create a bucket with a different name

```bash
$ gcloud storage buckets create gs://rnaseq-pipeline-bckt
```

Now it's time for Nextflow to join the party! 🥳

### Setting up Nextflow to make use of Batch and Storage

#### Write the configuration file

Here you will set up a simple RNAseq pipeline with Nextflow to be run entirely on Google Cloud Platform (GCP) directly from your local machine.

Start by creating a folder for your project on your local machine, such as “rnaseq-example”. It's important to mention that you can also go fully cloud and use a Virtual Machine for everything we will do here locally.

Inside the folder that you created for the project, create a file named `nextflow.config` with the following content (remember to replace PROJECT-ID with the project id you created above):

```Nextflow
workDir = 'gs://rnaseq-pipeline-bckt/scratch'

process {
  executor = 'google-batch'
  container = 'nextflow/rnaseq-nf'
  errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
  maxRetries = 5
}

google {
  project = 'PROJECT-ID'
  location = 'us-central1'
  batch.spot = true
}
```

The `workDir` option tells Nextflow to use the bucket you created as the work directory. Nextflow will use this directory to stage our input data and store intermediate and final data. Nextflow does not allow you to use the root directory of a bucket as the work directory -- it must be a subdirectory instead. Using a subdirectory is also just a good practice.

The `process` scope tells Nextflow to run all the processes (steps) of your pipeline on Google Batch and to use the `nextflow/rnaseq-nf` Docker image hosted on DockerHub (default) for all processes. Also, the error strategy will automatically retry any failed tasks with exit code 14, which is the exit code for spot instances that were reclaimed.

The `google` scope is specific to Google Cloud. You need to provide the project id (don't provide the project name, it won't work!), and a Google Cloud location (leave it as above if you're not sure of what to put). In the example above, spot instances are also requested (more info about spot instances [here](https://www.nextflow.io/docs/latest/google.html#spot-instances)), which are cheaper instances that, as a drawback, can be reclaimed at any time if resources are needed by the cloud provider. Based on what we have seen so far, the `nextflow.config` file should contain "rnaseq-nxf" as the project id.

#### Write your Nextflow pipeline

With that done, it's time to write your Nextflow pipeline. The example below is a proof of concept RNAseq pipeline. More information about this pipeline and where this came from can be found in the public training material developed by Seqera Labs at [https://training.seqera.io/](https://training.seqera.io/). Copy the script below and save it as `main.nf` in the same folder as your nextflow.config file:

```Nextflow
/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"
log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

process QUANTIFICATION {
    tag "Salmon on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path salmon_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}

process FASTQC {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    index_ch = INDEX(params.transcriptome_file)
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(quant_ch.mix(fastqc_ch).collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
```

#### Download the input data

Create a folder named `data` in your project folder and download the example files (also taken from the Seqera training material).

```bash
$ mkdir data; cd data
$ wget https://raw.githubusercontent.com/seqeralabs/nf-training-public/master/nf-training/data/ggal/gut_1.fq
$ wget https://raw.githubusercontent.com/seqeralabs/nf-training-public/master/nf-training/data/ggal/gut_2.fq
$ wget https://raw.githubusercontent.com/seqeralabs/nf-training-public/master/nf-training/data/ggal/transcriptome.fa
$ cd ..
```

Files can also be uploaded directly to Google Cloud Storage, or downloaded directly from GitHub, depending on your own preferences.

Use the command below to authenticate with Google Cloud Platform. Nextflow will use this account by default when you run a pipeline.

```bash
$ gcloud auth application-default login
```

#### Launch the pipeline!

With that done, you’re now ready to run the example Nextflow pipeline:

```bash
$ nextflow run main.nf
```

While the pipeline stores everything in the bucket, our example pipeline will also download the final outputs to a local directory called `results`, because of how the `publishDir` directive was specified in the `main.nf` script. If you want to avoid the egress cost associated with downloading data from a bucket, you can change the `publishDir` to another bucket directory, e.g. `gs://rnaseq-pipeline-bckt/results`.

In your terminal, you should see something like this:

![Nextflow ongoing run on Google Cloud Batch](/img/ongoing-nxf-gbatch.png)

You can check the status of your jobs on Google Batch by opening another terminal and running the following command:

```bash
$ gcloud batch jobs list
```

By the end of it, if everything worked well, you should see something like:

![Nextflow run on Google Cloud Batch finished](/img/nxf-gbatch-finished.png)

And that's all, folks! 😆

You will find more information about Nextflow on Google Batch in [this blog post](https://www.nextflow.io/blog/2022/deploy-nextflow-pipelines-with-google-cloud-batch.html) and the [official Nextflow documentation](https://www.nextflow.io/docs/latest/google.html).

Special thanks to Hatem Nawar, Chris Hakkaart, and Ben Sherman for providing valuable feedback to this document.
