title=Deploy Nextflow Pipelines with Google Cloud Batch!
date=2022-07-13
type=post
description=Deploy Nextflow pipelines at scale with new the Google Cloud Batch compute service.
image=img/google_cloud_batch_nextflow-min.png
tags=nextflow,google,cloud
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~

A key feature of Nextflow is the ability to abstract the implementation of data analysis pipelines so they can be deployed in a portable manner across execution platforms. 

As of today, Nextflow supports a rich variety of HPC schedulers and all major cloud providers. Our goal is to support new services as they emerge to enable Nextflow users to take advantage of the latest technology and deploy pipelines on the compute environments that best fit their requirements. 

For this reason, we are delighted to announce that Nextflow now supports [Google Cloud Batch](https://cloud.google.com/batch), a new fully managed batch service just announced for beta availability by Google Cloud.

### A New On-Ramp to the Google Cloud

Google Cloud Batch is a comprehensive cloud service suitable for multiple use cases, including HPC, AI/ML, and data processing. While it is similar to the Google Cloud Life Sciences API, used by many Nextflow users today, Google Cloud Batch offers a broader set of capabilities. As with Google Cloud Life Sciences, Google Cloud Batch automatically provisions resources, manages capacity, and allows batch workloads to run at scale. It offers several advantages, including:

* The ability to re-use VMs across jobs steps to reduce overhead and boost performance.
* Granular control over task execution, compute, and storage resources.
* Infrastructure, application, and task-level logging.
* Improved task parallelization, including support for multi-node MPI jobs, with support for array jobs, and subtasks.
* Improved support for spot instances, which provides a significant cost saving when compared to regular instance.
* Streamlined data handling and provisioning.

A nice feature of Google Cloud Batch API, that fits nicely with Nextflow, is its built-in support for data ingestion from Google Cloud Storage buckets. A batch job can *mount* a storage bucket and make it directly accessible to a container running a Nextflow task. This feature makes data ingestion and sharing resulting data sets more efficient and reliable than other solutions. 

### Getting started with Google Cloud Batch

Support for the Google Cloud Batch requires the latest release of Nextflow from the edge channel (version `22.07.1-edge` or later). If you don't already have it, you can install this release using these commands:

```
export NXF_EDGE=1 
curl get.nextflow.io | bash
./nextflow -self-update
```

Make sure your Google account is allowed to access the Google Cloud Batch service by checking the [API & Service](https://console.cloud.google.com/apis/dashboard) dashboard. 

Credentials for accessing the service are picked up by Nextflow from your environment using the usual [Google Application Default Credentials](https://github.com/googleapis/google-auth-library-java#google-auth-library-oauth2-http) mechanism. That is, either via the `GOOGLE_APPLICATION_CREDENTIALS` environment variable, or by using the following command to set up the environment: 

```
gcloud auth application-default login
```

After authenticating yourself to Google Cloud, create a `nextflow.config` file and specify `google-batch` as the Nextflow executor. You will also need to specify the Google Cloud project where execution will occur and the Google Cloud Storage working directory for pipeline execution. 

```
cat <<EOT > nextflow.config
process.executor = 'google-batch'
workDir = 'gs://YOUR-GOOGLE-BUCKET/scratch'
google.project = 'YOUR GOOGLE PROJECT ID'
EOT
```

In the above snippet replace `<YOUR_GOOGLE_BUCKET>` with a Google Storage bucket of your choice where to store the pipeline output data and `<YOUR_GOOGLE_PROJECT_ID>` with your Google project Id where the computation will be deployed. 

With this information, you are ready to start. You can verify that the integration is working by running the Nextflow “hello” pipeline as shown below:

```
nextflow run https://github.com/nextflow-io/hello
```


### Migrating Google Cloud Life Sciences pipelines to Google Cloud Batch 

Google Cloud Life Sciences users can easily migrate their pipelines to Google Cloud Batch by making just a few edits to their pipeline configuration settings. Simply replace the `google-lifesciences` executor with `google-batch`. 

For each setting having the prefix `google.lifeScience.`, there is a corresponding  `google.batch.` setting. Simply update these configuration settings to reflect the new service.


The usual process directives such as: [cpus](https://www.nextflow.io/docs/latest/process.html#cpus), [memory](https://www.nextflow.io/docs/latest/process.html#memory), [time](https://www.nextflow.io/docs/latest/process.html#time), [machineType](https://www.nextflow.io/docs/latest/process.html#machinetype) are natively supported by Google Cloud Batch, and should not be modified. 


### 100% Open, Built to Scale

The Google Cloud Batch executor for Nextflow is offered as an open source contribution to the Nextflow project. The integration was developed by Google in collaboration with [Seqera Labs](https://seqera.io/). This is a validation of Google Cloud’s ongoing commitment to open source software (OSS) and a testament to the health and vibrancy of the Nextflow project. We wish to thank the entire Google Cloud Batch team, and Shamel Jacobs in particular, for their support of this effort. 

### Conclusion 

Support for Google Cloud Batch further expands the wide range of computing platforms supported by Nextflow. It empowers Nextflow users to easily access cost-effective resources, and take full advantage of the rich capabilities of the Google Cloud. Above all, it enables researchers to easily scale and collaborate, improving their productivity, and resulting in better research outcomes.
