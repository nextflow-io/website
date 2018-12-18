title=Bringing Nextflow to Google Cloud Platform with WuXi NextCODE
date=2018-12-18
type=post
tags=nextflow,wuxinextcode,google,cloud
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~

<div class='text-muted' style='margin-bottom:2em'>
<i>This is a guest post authored by Halli Bjornsson, Head of Product Development Operations at WuXi NextCODE and Jonathan Sheffi, Product Manager, Biomedical Data at Google Cloud.
</i>
</div>


Google Cloud and WuXi NextCODE are dedicated to advancing the state of the art in biomedical informatics, especially through open source, which allows developers to collaborate broadly and deeply.

WuXi NextCODE is itself a user of Nextflow, and Google Cloud has many customers that use Nextflow. Together, we’ve collaborated to deliver Google Cloud Platform (GCP) support for Nextflow using the [Google Pipelines API](https://cloud.google.com/genomics/pipelines). Pipelines API is a managed computing service that allows the execution of containerized workloads on GCP.

<div class="row">
  <div class="column">
  <img src='/img/google-cloud.svg' style='width:80%; padding:1.7em; '>
  </div>
  <div class="column">
  <img src='/img/wuxi-nextcode.jpeg' style='width:80%; padding:1em'>
  </div>
</div>
<style>
.column {
  float: left;
  width: 50%;
}

.row:after {
  content: "";
  display: table;
  clear: both;
}
</style>

Nextflow now provides built-in support for Google Pipelines API which allows the seamless deployment of a Nextflow pipeline in the cloud, offloading the process executions as pipelines running on Google's scalable infrastructure with a few commands. This makes it even easier for customers and partners like WuXi NextCODE to process biomedical data using Google Cloud.


### Get started! 

This feature is currently available in the Nextflow edge channel. Follow these steps to get started:

* Install Nextflow from the edge channel exporting the variables shown below and then running the usual Nextflow installer Bash snippet:

    ```
    export NXF_VER=18.12.0-edge
    export NXF_MODE=google
    curl https://get.nextflow.io | bash
    ```

* [Enable the Google Genomics API for your GCP projects](https://console.cloud.google.com/flows/enableapi?apiid=genomics.googleapis.com,compute.googleapis.com,storage-api.googleapis.com). 

* [Download and set credentials for your Genomics API-enabled project](https://cloud.google.com/docs/authentication/production#obtaining_and_providing_service_account_credentials_manually).

* Change your `nextflow.config` file to use the Google Pipelines executor and specify the required config values for it as [described in the documentation](/docs/edge/google.html#google-pipelines).

* Finally, run your script with Nextflow like usual, specifying a Google Storage bucket as the pipeline work directory with the `-work-dir` option. For example:

    ```
    nextflow run rnaseq-nf -work-dir gs://your-bucket/scratch
    ```

<br>
You can find more detailed info about available configuration settings and deployment options at [this link](/docs/edge/google.html).

We’re thrilled to make this contribution available to the Nextflow community!