title=Introducing Nextflow for Azure Batch
date=2021-02-22
type=post
tags=nextflow,azure
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~


When the Nextflow project was created, one of the main drivers was to enable reproducible data pipelines that could be deployed across a wide range of execution platforms with minimal effort as well as to empower users to scale their data analysis while facilitating the migration to the cloud.

Throughout the years, the computing services provided by cloud vendors have evolved in a spectacular manner. Eight years ago, the model was focused on launching virtual machines in the cloud, then came containers and then the idea of serverless computing which changed everything again. However, the power of the Nextflow abstraction consists of hiding the complexity of the underlying platform. Through the concept of executors, emerging technologies and new platforms can be easily adapted with no changes required to user pipelines.

With this in mind, we could not be more excited to announce that over the past months we have been working with Microsoft to implement built-in support for [Azure Batch](https://azure.microsoft.com/en-us/services/batch/) into Nextflow. Today we are delighted to make it available to all users as a beta release. 


### How does it work

Azure Batch is a cloud-based computing service that allows the execution of highly scalable, container based, workloads in the Azure cloud.

The support for Nextflow comes in the form of a plugin which implements a new executor, not surprisingly named `azurebatch`, which offloads the execution of the pipeline jobs to corresponding Azure Batch jobs. 

Each job run consists in practical terms of a container execution which ships the job dependencies and carries out the job computation. As usual, each job is assigned a unique working directory allocated into a [Azure Blob](https://azure.microsoft.com/en-us/services/storage/blobs/) container.

### Let's get started!

The support for Azure Batch requires the latest release of Nextflow from the *edge* channel (version 21.02-edge or later). If you don't have this, you can install it using these commands:

```
export NXF_EDGE=1 
curl get.nextflow.io | bash
```

Note for Windows users, as Nextflow is *nix based tool you will need to run it using the [Windows subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10). Also make sure Java 8 or later is installed in the Linux environment.

Once Nextflow is installed, to run your data pipelines with Azure Batch, you will need to create an Azure Batch account in the region of your choice using the Azure Portal. In a similar manner, you will need an Azure Blob container.

With the Azure Batch and Blob storage container configured, your `nextflow.config` file should be set up similar to the example below:

```
plugins {
  id 'nf-azure'
}

process {
  executor = 'azurebatch'
}

azure {
  batch {
    location = 'westeurope'
    accountName = '<YOUR BATCH ACCOUNT NAME>'
    accountKey = '<YOUR BATCH ACCOUNT KEY>'
    autoPoolMode = true
  }
  storage {
    accountName = "<YOUR STORAGE ACCOUNT NAME>"
    accountKey = "<YOUR STORAGE ACCOUNT KEY>"
  }
}
```

Using this configuration snippet, Nextflow will automatically create the virtual machine pool(s) required to deploy the pipeline execution in the Azure Batch service.

Now you will be able to launch the pipeline execution using the following command: 

```
nextflow run <pipeline name> -w az://my-container/work
```

Replace `<pipeline name>` with a pipeline name e.g. nextflow-io/rnaseq-nf and `my-container` with a blob container in the storage account as defined in the above configuration.

For more details regarding the Nextflow configuration setting for Azure Batch 
refers to the Nextflow documentation at [this link](/docs/edge/azure.html). 

### Conclusion

The support for Azure Batch further expands the wide range of computing platforms supported by Nextflow and empowers Nextflow users to deploy their data pipelines in the cloud provider of their choice. Above all, it allows researchers to scale, collaborate and share their work without being locked into a specific platform.

We thank Microsoft, and in particular  [Jer-Ming Chia](https://www.linkedin.com/in/jermingchia/) who works in the HPC and AI team for having supported and sponsored this open source contribution to the Nextflow framework.
