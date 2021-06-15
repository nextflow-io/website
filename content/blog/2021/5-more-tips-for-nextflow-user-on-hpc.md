title=Five more tips for Nextflow user on HPC
date=2021-06-15
type=post
tags=nextflow,hpc
status=published
author=Kevin Sayers
icon=kevin.jpg
~~~~~~

In May we blogged about [Five Nextflow Tips for HPC Users](/blog/2021/5_tips_for_hpc_users.html) and now we continue the series with five additional tips for deploying Nextflow with on HPC batch schedulers.

### 1. Use the scratch directive

To allow the pipeline tasks to share data with each other, Nextflow requires a shared file system path as a working directory. When using this model, a common recommendation is to use the node's local scratch storage as the job working directory to avoid unnecessary use of the network shared file system and achieve better performance.

Nextflow implements this best-practice which can be enabled by adding the following setting in your `nextflow.config` file.

```
process.scratch = true
```

When using this option, Nextflow:
* Creates a unique directory in the computing node's local `/tmp` or the path assigned by your cluster via the `TMPDIR` environment variable.
* Creates a [symlink](https://en.wikipedia.org/wiki/Symbolic_link) for each input file required by the job execution.
* Runs the job in the local scratch path.
Copies the job output files into the job shared work directory assigned by Nextflow.

### 2. Use -bg option to launch the execution in the background 

In some circumstances, you may need to run your Nextflow pipeline in the background without losing the execution output. In this scenario use the `-bg` command line option as shown below. 

```
nextflow run <pipeline> -bg > my-file.log
``` 

This can be very useful when launching the execution from an SSH connected terminal and ensures that any connection issues don't stop the pipeline. You can use `ps` and `kill` to find and stop the execution.

### 3. Disable interactive logging 

Nextflow has rich terminal logging which uses ANSI escape codes to update the pipeline execution counters interactively. However, this is not very useful when submitting the pipeline execution as a cluster job or in the background. In this case, disable the rich ANSI logging using the command line option `-ansi-log false` or the environment variable `NXF_ANSI_LOG=false`.

### 4. Cluster native options

Nextlow has portable directives for common resource requests such as [cpus](https://www.nextflow.io/docs/latest/process.html#cpus), [memory](https://www.nextflow.io/docs/latest/process.html#memory) and [disk](https://www.nextflow.io/docs/latest/process.html#disk) allocation. 

These directives allow you to specify the request for a certain number of computing resources e.g CPUs, memory, or disk and Nextflow converts these values to the native setting of the target execution platform specified in the pipeline configuration.

However, there can be settings that are only available on some specific cluster technology or vendors. 

The [clusterOptions](https://www.nextflow.io/docs/latest/process.html#clusterOptions) directive allows you to specify any option of your resource manager for which there isn't direct support in Nextflow. 


### 5. Retry failing jobs increasing resource allocation 

A common scenario is that instances of the same process may require different computing resources. For example, requesting an amount of memory that is too low for some processes will result in those tasks failing. You could specify a higher limit which would accommodate the task with the highest memory utilization, but you then run the risk of decreasing your job’s execution priority. 

Nextflow provides a mechanism that allows you to modify the amount of computing resources requested in the case of a process failure and attempt to re-execute it using a higher limit. For example:

```
process foo {

  memory { 2.GB * task.attempt }
  time { 1.hour * task.attempt }

  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  script:
  """
  your_job_command --here
  """
}
```

In the above example the memory and execution time limits are defined dynamically. The first time the process is executed the task.attempt is set to 1, thus it will request 2 GB of memory and one hour of maximum execution time.

If the task execution fails, reporting an exit status in the range between 137 and 140, the task is re-submitted (otherwise it terminates immediately). This time the value of task.attempt is 2, thus increasing the amount of the memory to four GB and the time to 2 hours, and so on.

NOTE: These exit statuses are not standard and can change depending on the resource manager you are using. Consult your cluster administrator or scheduler administration guide for details on the exit statuses used by your cluster in similar error conditions. 


### Conclusion

Nextflow aims to give you control over every aspect of your workflow. These Nextflow options allow you to shape how Nextflow submits your processes to your executor, that can make your workflow more robust by avoiding the overloading of the executor. Some systems have hard limits which if you do not take into account, no processes will be executed. Being aware of these configuration values and how to use them is incredibly helpful when working with larger workflows. 
