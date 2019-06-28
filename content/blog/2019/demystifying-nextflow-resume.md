title=Demystifying Nextflow resume
date=2019-06-24
type=post
tags=nextflow,resume
status=published
author=Evan Floden
icon=evan.jpg
~~~~~~

*This two-part blog aims to help users understand Nextflow’s powerful caching mechanism. Part one describes how it works whilst part two will focus on execution provenance and troubleshooting. You can read part two [here](blog/2019/troubleshooting-nextflow-resume.html)*

Task execution caching and checkpointing is an essential feature of any modern workflow manager and Nextflow provides an automated caching mechanism with every workflow execution. When using the `-resume` flag, successfully completed tasks are skipped and the previously cached results are used in downstream tasks. But understanding the specifics of how it works and debugging situations when the behaviour is not as expected is a common source of frustration.

The mechanism works by assigning a unique ID to each task. This unique ID is used to create a separate execution directory, called the working directory, where the tasks are executed and the results stored. A task’s unique ID is generated as a 128-bit hash number obtained from a composition of the task’s: 

* Inputs values
* Input files
* Command line string
* Container ID
* Conda environment
* Environment modules
* Any executed scripts in the bin directory


### How does resume work?
The `-resume` command line option allows for the continuation of a workflow execution. It can be used in its most basic form with:

```
$ nextflow run nextflow-io/hello -resume
```

In practice, every execution starts from the beginning. However, when using resume, before launching a task, Nextflow uses the unique ID to check if:

* the working directory exists 
* it contains a valid command exit status
* it contains the expected output files.

If these conditions are satisfied, the task execution is skipped and the previously computed outputs are applied. When a task requires recomputation, ie. the conditions above are not fulfilled, the downstream tasks are automatically invalidated.

### The working directory
By default, the task work directories are created in the directory from where the pipeline is launched. This is often a scratch storage area that can be cleaned up once the computation is completed. A different location for the execution work directory can be specified using the command line option `-w` e.g.

```
$ nextflow run <script> -w /some/scratch/dir
```

Note that if you delete or move the pipeline work directory, this will prevent to use the resume feature in subsequent runs. 

Also note that the pipeline work directory is intended to be used as a temporary scratch area. The final 
workflow outputs are expected to be stored in a different location specified using the [`publishDir` directive](https://www.nextflow.io/docs/latest/process.html#publishdir). 

### How is the hash calculated on input files?
The hash provides a convenient way for Nextflow to determine if a task requires recomputation. For each input file, the hash code is computed with:

* The complete file path
* The file size
* The last modified timestamp

Therefore, even just performing a touch on a file will invalidate the task execution.

### How to ensure resume works as expected?
It is good practice to organize each experiment in its own folder. An experiment’s input parameters can be specified using a Nextflow config file which also makes it simple to track and replicate an experiment over time. Note that you should avoid launching two (or more) Nextflow instances in the same directory concurrently.

The nextflow log command lists the executions run in the current folder:

<style>
pre {
    white-space: pre;
    overflow-x: auto;
}
</style>
<pre>
$ nextflow log

TIMESTAMP            DURATION  RUN NAME          STATUS  REVISION ID  SESSION ID                            COMMAND                                    
2019-05-06 12:07:32  1.2s      focused_carson    ERR     a9012339ce   7363b3f0-09ac-495b-a947-28cf430d0b85  nextflow run hello                         
2019-05-06 12:08:33  21.1s     mighty_boyd       OK      a9012339ce   7363b3f0-09ac-495b-a947-28cf430d0b85  nextflow run rnaseq-nf -with-docker        
2019-05-06 12:31:15  1.2s      insane_celsius    ERR     b9aefc67b4   4dc656d2-c410-44c8-bc32-7dd0ea87bebf  nextflow run rnaseq-nf                     
2019-05-06 12:31:24  17s       stupefied_euclid  OK      b9aefc67b4   4dc656d2-c410-44c8-bc32-7dd0ea87bebf  nextflow run rnaseq-nf -resume -with-docker
</pre>

You can use the resume command with either the session ID or the run name to recover a specific execution. For example:

```
$ nextflow run rnaseq-nf -resume mighty_boyd
```

or 

```
nextflow run naseq-nf -resume 4dc656d2-c410-44c8-bc32-7dd0ea87bebf
```


*Stay tuned for part two where we will discuss resume in more detail with respect to provenance and troubleshooting techniques!*