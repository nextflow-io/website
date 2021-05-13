title=5 Netxflow tips for HPC users
date=2021-05-13
type=post
tags=nextflow,hpc
status=published
author=Kevin Sayer
icon=kevin.jpg
~~~~~~

Nextflow is a powerful tool for developing scientific workflows for use on HPC systems. It provides a simple solution to deploy parallelized workloads at scale using an elegant reactive/function programming model in a portable manner. 

It supports the most popular workload managers such as GridEngine, Slurm, LSF and PBS, among other out-of-the-box executors, and comes with sensible defaults for each. However, each HPC system is a complex machine with its own characteristics and constraints. For this reason you should always consult your system admin before running a new piece of software or a compute intensive pipeline that spawns a large number of jobs.  

In this series of posts, we will be sharing the top tips we have learnt along the way that should help you get your results faster while keeping in the good books of your sysadmins.


### 1. Don't forget the executor  

Nextflow, by default, spawns parallel task executions in the computer in which it is running. This is generally useful for development purposes, however when using an HPC system you should specify the executor matching your system. This instructs Nextflow to submit pipeline tasks as jobs into your HPC workload manager. This can be done adding the following setting to the `nextflow.config` file in the launching directory, for example:

```
executor.slurm = 'slurm' 
```

With the above setting Nextflow will submit the job executions to your Slurm cluster spawning a `sbatch` command for each job in your pipeline. Find the executor matching your system at [this link](https://www.nextflow.io/docs/latest/executor.html).

Even better, to prevent the undesired use of the local executor in a specific environment, define the *default* executor to be used by Nextflow using the following system variable: 

```
export NXF_EXECUTOR=slurm
```

### 2. Nextflow as a job AKA don’t run in the login node! 

Quite surely your sysadmin already warned you that the login should only be used to submit job execution and not run compute intensive tasks. 
When running a Nextflow pipeline, the driver application submits and monitors the job executions of your cluster (provided you have correctly specified the executor as said at the point 1), and therefore it should not run compute intensive tasks. 

However it's never a good practice to launch a long running job in the login node, and therefore a good practice consists of running Nextflow itself as a cluster job. This can be done by wrapping the `nextflow run` command in a shell script and submitting it as any other job. An average pipeline may require 2 CPUs and 2 GB of resources allocation. 

Note: the queue where the Nextflow driver job is submitted should allow the spawning of the pipeline jobs to carry out the pipeline execution.

### 3. Specify queueSize

The `queueSize` directive is part of the executor configuration in the `nextflow.config` file, and defines how many processes are queued at a given time. By default Nextflow will submit up to  100 jobs at a time for execution. Increase or decrease this setting depending your HPC system quota and throughput, for example:
 
```
executor {
    name = 'slurm'
    queueSize = 50
}

```

### 4. Specify the max heap size 

The Nextflow runtime runs on top of the Java virtual machine which, by design, tries to allocate as much memory as is available. This is not a good practice in HPC systems which are made to share compute resources across many users and applications. 
To avoid this, specify the max amount of memory that can be used by the Java VM using the -Xms and -Xmx Java flags. These can be specified using the `NXF_OPTS` environment variable. 

For example: 

```
export NXF_OPTS="-Xms500M -Xmx2G" 
```

The above setting instructs Nextflow to use allocate a Java heap in the range of 500 MB and 2 GB of RAM.

### 5. Specify the submit rate limit 

Nextflow tries to submit the job executions as quickly as possible, which is generally not a problem. However, in some HPC systems the submission throughput is constrained or it should be limited to avoid degrading the overall system performance.

To prevent this problem you can use `submitRateLimit` to control the Nextflow job submission throughput. This directive is part of the `executor` configuration scope, and defines the number of tasks that can be submitted per a unit of time. The default for the `submitRateLimit` is unlimited. 

You can specify the `submitRateLimit` like this:

```
executor {
	submitRateLimit = ‘10 sec’
}
```

You  can also more explicitly specify it as a rate of # processes / time unit: 

```
executor {
	submitRateLimit = ‘10/2min’
}
```

### Conclusion

Nextflow aims to give you control over every aspect of your workflow. These options allow you to shape how Nextflow communicates with your HPC system. This can make workflows more robust while avoiding overloading the executor. Some systems have hard limits, which if you do not take into account, will stop any jobs from being scheduled. 

Stay tuned for part two where we will discuss background executions, retry strategies, maxForks and other tips. 

