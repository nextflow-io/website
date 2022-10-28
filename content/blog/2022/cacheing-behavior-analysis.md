title=Analysing cacheing behavior of pipelines
date=2022-10-28
type=post
description=A guide to analysis of dump-hashes for troubleshooting cacheing behavior.
image=img/delta-hash-dump-diff.png
tags=nextflow,cache
status=published
author=Abhinav Sharma
icon=abhinav.jpg
~~~~~~

<!-- NOTE: Talk about the hash-computation funnel?  -->
<!-- The following code snippet from the Nextflow code base  -->
<!-- https://github.com/nextflow-io/nextflow/blob/665dae3a5ec63b8a990d691a3441a36468847e1f/plugins/nf-azure/src/main/nextflow/cloud/azure/config/AzPoolOpts.groovy#L91-L109 -->

Resumability of analysis (i.e. cacheing) is one of the core strengths of Nextflow and allows us to avoid re-running the processes which have not been affected by any new change in the pipeline code by simply appending `-resume` to the Nextflow `run` command. Depending upon the nature of the change introduced sometimes it becomes necessary to understand why exactly a specific process was rerun.

We have previously written about Nextflow's [resume functionality](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html) as well as highlighted some [troubleshooting strategies](https://www.nextflow.io/blog/2019/troubleshooting-nextflow-resume.html) to gain more insights on the cacheing behavior. 

In this blog post we will take a more hands on approach and highlight the strategies which we can use to understand what is causing a particular process (or processes) to rerun and therefore not rely upon the cache from previous runs of the pipeline. We will use the standard proof-of-concept [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) and introduce a minor change in one of the process definitions and investigate how it affects the overall cacheing behavior when compared to the initial execution of the pipeline.

## Pre-requisites for the analysis

This blog post makes use of the baseline Nextflow (and Java) setup along with Docker to provide a self-contained example.

- Nextflow

```terminal
>_ nextflow info

  Version: 22.10.0 build 5826
  Created: 13-10-2022 05:44 UTC (07:44 SAST)
  System: Mac OS X 12.6
  Runtime: Groovy 3.0.13 on OpenJDK 64-Bit Server VM 17+35-LTS
  Encoding: UTF-8 (UTF-8)

```
- Java-17

```terminal
>_ java -version

  openjdk version "17" 2021-09-14 LTS
  OpenJDK Runtime Environment Corretto-17.0.0.35.1 (build 17+35-LTS)
  OpenJDK 64-Bit Server VM Corretto-17.0.0.35.1 (build 17+35-LTS, mixed mode, sharing)


```

- Docker

```terminal
>_ docker version

  Client:
   Version:           20.10.8

  Server: Docker Engine - Community
   Engine:
    Version:          20.10.8

```

### Local setup for the test

We clone the standard proof-of-concept [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf
) pipeline locally, so that we can introduce a change and observe its effect on cacheing behavior of subsequent runs.

```
>_ git clone https://github.com/nextflow-io/rnaseq-nf
>_ cd rnaseq-nf
```

### Pipeline flowchart

The flowchart below would help in understanding the design of the pipeline and the dependencies between the various tasks.

![rnaseq-nf](/img/rnaseq-nf.base.png)

### Logs from initial (fresh) run

To analyze the logs for a subsequent run, we need to have the initial hashes for the processes in the pipeline which would be captured in `fresh_run.log` and would be used later on as "ground-truth" for the analysis.

**TIP:** We rely upon the [`-log` option](https://www.nextflow.io/docs/latest/cli.html#execution-logs) in the `nextflow` CLI to use a custom log file name instead of the default `.nextflow.log`.

```
>_ nextflow -log fresh_run.log run ./main.nf -profile docker -dump-hashes

...
...
executor >  local (4)
[f2/160b79] process > RNASEQ:INDEX (ggal_1_48850000_49020000) [100%] 1 of 1 ✔
[82/e3c268] process > RNASEQ:FASTQC (FASTQC on ggal_gut)      [100%] 1 of 1 ✔
[48/e4ae12] process > RNASEQ:QUANT (ggal_gut)                 [100%] 1 of 1 ✔
[fa/e0a8a6] process > MULTIQC                                 [100%] 1 of 1 ✔

...
...

```

### Introduce a change within the `modules/fastqc.nf`

After the initial run of the pipeline, we introduce a change in the `fastqc.nf` module and hard-code the number of threads which should be used to run the `FASTQC` process via Nextflow's [`cpus` directive](https://www.nextflow.io/docs/latest/process.html#cpus).


Here's the output of `git diff` on the contents of `modules/fastqc/main.nf` file after the introduction the change.

```diff
>_ git diff

diff --git a/modules/fastqc/main.nf b/modules/fastqc/main.nf
index c235b25..bdff6dc 100644
--- a/modules/fastqc/main.nf
+++ b/modules/fastqc/main.nf
@@ -4,6 +4,7 @@ process FASTQC {
     tag "FASTQC on $sample_id"
     conda 'bioconda::fastqc=0.11.9'
     publishDir params.outdir, mode:'copy'
+    cpus 2

     input:
     tuple val(sample_id), path(reads)
@@ -13,6 +14,6 @@ process FASTQC {

     script:
     """
-    fastqc.sh "$sample_id" "$reads"
+    fastqc.sh "$sample_id" "$reads" -t ${task.cpus}
     """
 }


```

### Logs from the follow up run

In this run, we will add the `-resume` option and instruct Nextflow to rely upon the cached results from previous run and therefore only run the parts of the pipeline which have changed. We also instruct Nextflow to print the hashes in the `resumed_run.log` file.

```
>_ nextflow -log resumed_run.log run ./main.nf -profile docker -dump-hashes -resume

...
...
executor >  local
[f2/160b79] process > RNASEQ:INDEX (ggal_1_48850000_49020000) [100%] 1 of 1, cached: 1 ✔
[48/d8d964] process > RNASEQ:FASTQC (FASTQC on ggal_gut)      [100%] 1 of 1 ✔
[48/e4ae12] process > RNASEQ:QUANT (ggal_gut)                 [100%] 1 of 1, cached: 1 ✔
[c6/0c6b68] process > MULTIQC                                 [100%] 1 of 1 ✔
...
...


```

## Analysis of cache hashes

In the resumed run, we noticed that the `RNASEQ:FASTQC (FASTQC on ggal_gut)` and `MULTIQC` processes were rerun while the others were cached. To understand the reason why both of these processes were triggered again, we can examine the hashes generated by the processes from the logs of the `fresh_run` and `resumed_run`.

For the analysis, we need to keep in mind that:

1. The time-stamps would of course be different, therefore we can safely ignore the timestamps and narrow down the `grep` pattern to the Nextflow `TaskProcessor` class.

2. The order of the log entries isn't fixed due to the nature of the underlying parallel computation dataflow model which is used by Nextflow. For example, `FASTQC` in `fresh_run.log` isn't the first logged processes in `resumed_run.log`

**NOTE**: The use of classic Unix tools like `cut` / `grep` / `sort` in the commands address the points above.

### Find the process level hashes

- `fresh_run.log`

```
rnaseq-nf  🍣 master 🅒 base
>_ cat ./fresh_run.log | grep 'INFO.*TaskProcessor.*cache hash' | cut -d '-' -f 3 | cut -d ';' -f 1 | sort | tee  ./fresh_run.tasks.log

 [MULTIQC] cache hash: 167d7b39f7efdfc49b6ff773f081daef
 [RNASEQ:FASTQC (FASTQC on ggal_gut)] cache hash: 47e8c58d92dbaafba3c2ccc4f89f53a4
 [RNASEQ:INDEX (ggal_1_48850000_49020000)] cache hash: ac8be293e1d57f3616cdd0adce34af6f
 [RNASEQ:QUANT (ggal_gut)] cache hash: d8b88e3979ff9fe4bf64b4e1bfaf4038

```

- `resumed_run.log`

```
>_ cat ./resumed_run.log | grep 'INFO.*TaskProcessor.*cache hash' | cut -d '-' -f 3 | cut -d ';' -f 1 | sort | tee  ./resumed_run.tasks.log

 [MULTIQC] cache hash: d3f200c56cf00b223282f12f06ae8586
 [RNASEQ:FASTQC (FASTQC on ggal_gut)] cache hash: 92478eeb3b0ff210ebe5a4f3d99aed2d
 [RNASEQ:INDEX (ggal_1_48850000_49020000)] cache hash: ac8be293e1d57f3616cdd0adce34af6f
 [RNASEQ:QUANT (ggal_gut)] cache hash: d8b88e3979ff9fe4bf64b4e1bfaf4038

```

### Inference from process top-level hashes

Computing a hash is a multi-step process and various factors contribute to it, such as the inputs of the process, platform, time-stamps of the input files etc. All of these factors are used to compute the final process-level hash, therefore the change we introduced in the task level CPUs directive as well as the script section of the `FASTQC` process triggered a re-computation of hashes.

Upon comparing these hashes, we can confirm that the difference in computed process-level hashes validates the observed behavior.

```diff
>_ diff ./fresh_run.tasks.log ./resumed_run.tasks.log

1,2c1,2
<  [MULTIQC] cache hash: 167d7b39f7efdfc49b6ff773f081daef
<  [RNASEQ:FASTQC (FASTQC on ggal_gut)] cache hash: 47e8c58d92dbaafba3c2ccc4f89f53a4
---
>  [MULTIQC] cache hash: d3f200c56cf00b223282f12f06ae8586
>  [RNASEQ:FASTQC (FASTQC on ggal_gut)] cache hash: 92478eeb3b0ff210ebe5a4f3d99aed2d

```

Note that even if we only introduced changes in `FASTQC`, the `MULTIQC` process was re-run since it relies upon the output of the FASTQC process, as is higlighted in the pipeline flowchart after we introduced the change in `FASTQC` process and how that affects the downstream `MULTIQC` process.

![rnaseq-nf after modification](/img/rnaseq-nf.modified.png)


### Understanding why FASTQC was re-run

For this process, we inspect the relevant `FASTQC` section *along* with the contributing factors such as hashes for its `input`, `script` and other directives.

Here's what the full `FASTQC` process hashes look like from within the `fresh_run.log` file

```
INFO  nextflow.processor.TaskProcessor - [RNASEQ:FASTQC (FASTQC on ggal_gut)] cache hash: 47e8c58d92dbaafba3c2ccc4f89f53a4; mode: STANDARD; entries: 
  6c11792ef32523f19ac69f13766387cf [java.util.UUID] 4aade6ae-67d3-4e34-9c67-3ecb32a5e4fa 
  195c7faea83c75f2340eb710d8486d2a [java.lang.String] RNASEQ:FASTQC 
  43e5a23fc27129f92a6c010823d8909b [java.lang.String] """
    fastqc.sh "$sample_id" "$reads"
    """
 
  8e58c0cec3bde124d5d932c7f1579395 [java.lang.String] quay.io/nextflow/rnaseq-nf:v1.1 
  7ec7cbd71ff757f5fcdbaa760c9ce6de [java.lang.String] sample_id 
  16b4905b1545252eb7cbfe7b2a20d03d [java.lang.String] ggal_gut 
  553096c532e666fb42214fdf0520fe4a [java.lang.String] reads 
  3f0553d854f79285fa4db56cf13d62b2 [nextflow.util.ArrayBag] [FileHolder(sourceObj:/home/abhi18av/rnaseq-nf/data/ggal/ggal_gut_1.fq, storePath:/home/abhi18av/rnaseq-nf/data/ggal/ggal_gut_1.fq, stageName:ggal_gut_1.fq), FileHolder(sourceObj:/home/abhi18av/rnaseq-nf/data/ggal/ggal_gut_2.fq, storePath:/home/abhi18av/rnaseq-nf/data/ggal/ggal_gut_2.fq, stageName:ggal_gut_2.fq)] 
  4f9d4b0d22865056c37fb6d9c2a04a67 [java.lang.String] $ 
  16fe7483905cce7a85670e43e4678877 [java.lang.Boolean] true 
  4a7103d7455de7dceca72e305e2d9382 [sun.nio.fs.UnixPath] /home/abhi18av/rnaseq-nf/bin/fastqc.sh 

```

When we compare the hashes between `fresh_run.log` and `resumed_run.log`, we see the following diff. 

TIP: To get a word level diff you can rely upon [`delta`](https://dandavison.github.io/delta/introduction.html) diff tool

- From the hash entries we can see that the content of the script has changed, `$task.cpus` has been added
- This caused another entry in the `resumed_run.log` regarding the content of the process level directive `cpus` as well as the updated script


![diff in hashes](/img/delta-hash-dump-diff.png)


### Understanding why `MULTIQC` was re-run

Now, we apply the same analysis technique for the `MULTIQC` process in both log files

```diff
- [MULTIQC] cache hash: 79f6f63ff816a3aecb24ff7a0499b669; mode: STANDARD; entries: 
+ [MULTIQC] cache hash: 7b05eb097cc35b00610f5ccb66773918; mode: STANDARD; entries: 
34c33
-   fac1d0c29da59a72bd617f4ecda550b1 [nextflow.util.ArrayBag] [FileHolder(sourceObj:/home/abhi18av/rnaseq-nf/work/74/08e7523f43026dd7ac80ba7e57df95/ggal_gut, storePath:/home/abhi18av/rnaseq-nf/work/74/08e7523f43026dd7ac80ba7e57df95/ggal_gut, stageName:ggal_gut), FileHolder(sourceObj:/home/abhi18av/rnaseq-nf/work/c0/6d41e0a2d9ee3f0083ef9834a41015/fastqc_ggal_gut_logs, storePath:/home/abhi18av/rnaseq-nf/work/c0/6d41e0a2d9ee3f0083ef9834a41015/fastqc_ggal_gut_logs, stageName:fastqc_ggal_gut_logs)] 
+   e66b49bd490b43c99088a93058c21c0d [nextflow.util.ArrayBag] [FileHolder(sourceObj:/home/abhi18av/rnaseq-nf/work/74/08e7523f43026dd7ac80ba7e57df95/ggal_gut, storePath:/home/abhi18av/rnaseq-nf/work/74/08e7523f43026dd7ac80ba7e57df95/ggal_gut, stageName:ggal_gut), FileHolder(sourceObj:/home/abhi18av/rnaseq-nf/work/1e/687e9b5d894b83e5e0a38b0b8defbb/fastqc_ggal_gut_logs, storePath:/home/abhi18av/rnaseq-nf/work/1e/687e9b5d894b83e5e0a38b0b8defbb/fastqc_ggal_gut_logs, stageName:fastqc_ggal_gut_logs)] 
41d39

```

As we can see above, one of the input file paths has changed as a result of `FASTQC` being rerun and therefore, `MULTIQC` had to be rerun as well.


## Conclusion

Debugging the cacheing behavior of a pipeline can be tricky, however a systematic analysis approach can help in investigation what exactly is causing a particular process to be rerun again. 

One good takeway for analyzing larger datasets is to consider using the `-dump-hashes` parameter by default for any pipeline execution, which would save us the trouble of running the pipeline again just to obtain the hashes in the log file.

Call to action: we invite the community to building tooling for a better cache-debugging experience for Nextflow, perhaps an nf-cache plugin? 
Stay tuned for an upcoming blog post on how to extend and add new functionality within Nextflow using plugins.
