title=Analysing cacheing behavior of pipelines
date=2022-09-26
type=post
description=A guide to analysis of dump-hashes for troubleshooting cacheing behavior.
image=img/FIXME.png
tags=nextflow,cache
status=published
author=Abhinav Sharma
icon=abhinav.jpg

~~~~~~

TODO: Add mermaid diagram to explain dependencies

## Experimental setup

- Nextflow
- Java-11
- Docker

Here's the `nextflow` setup I relied upon

```
+  >_ nextflow info
  Version: 22.09.4-edge build 5793
  Created: 19-09-2022 19:09 UTC (21:09 SAST)
  System: Mac OS X 11.5
  Runtime: Groovy 3.0.13 on OpenJDK 64-Bit Server VM 17-ea+12
  Encoding: UTF-8 (UTF-8)

```

### Local setup for the test

```
+  >_ git clone https://github.com/nextflow-io/rnaseq-nf
+  >_ cd rnaseq-nf
```

### Logs from initial base run

`base_run.log`

```
+  >_ nextflow -log base_run.log run ./main.nf -profile docker -dump-hashes

...
...
process > RNASEQ:INDEX (ggal_1_48850000_49020000) [100%] 1 of 1 ✔
process > RNASEQ:FASTQC (FASTQC on ggal_gut)      [100%] 1 of 1 ✔
process > RNASEQ:QUANT (ggal_gut)                 [100%] 1 of 1 ✔
process > MULTIQC                                 [100%] 1 of 1 ✔

...
...

```

### Introduce a change within the `modules/fastqc.nf`

```diff
diff --git a/modules/fastqc.nf b/modules/fastqc.nf
index 571519b..afc19cc 100644
--- a/modules/fastqc.nf
+++ b/modules/fastqc.nf
@@ -12,6 +12,6 @@ process FASTQC {

     script:
     """
-    fastqc.sh "$sample_id" "$reads"
+    fastqc.sh "$sample_id" "$reads" $task.cpus
     """
 }


```

### Logs from the follow up run

`resumed_run.log`

```
+  >_ nextflow -log resumed_run.log run ./main.nf -profile docker -dump-hashes -resume

...
...
process > RNASEQ:INDEX (ggal_1_48850000_49020000) [100%] 1 of 1, cached: 1 ✔
process > RNASEQ:FASTQC (FASTQC on ggal_gut)      [100%] 1 of 1 ✔
process > RNASEQ:QUANT (ggal_gut)                 [100%] 1 of 1, cached: 1 ✔
process > MULTIQC                                 [100%] 1 of 1 ✔

...
...


```

## Analysis

In the resumed run, we noticed that the `RNASEQ:FASTQC (FASTQC on ggal_gut)` and `MULTIQC` processes were rerun while the others were cached. To understand this we can examine the hashes

Some things to keep in mind

- The time-stamps would of course be different, therefore we have narrowed down the `grep` pattern to the nextflow `TaskProcessor` class.

- The order of the log entries isn't fixed due to the parallel nature of computation. For example, `FASTQC` in `base_run.log` isn't the first logged processes in `resumed_run.log`

**NOTE**: The use of classic Unix tools like `cut` / `grep` / `sort` in the commands address the points above.

### Find the process level hashes

- `base_run.log`

```
rnaseq-nf  🍣 master 🅒 base
+  >_ cat ./base_run.log | grep 'INFO.*TaskProcessor.*cache hash' | cut -d '-' -f 3 | cut -d ';' -f 1 | sort | tee  ./base_run.tasks.log
 [MULTIQC] cache hash: 79f6f63ff816a3aecb24ff7a0499b669
 [RNASEQ:FASTQC (FASTQC on ggal_gut)] cache hash: 0c18fe122c17bae488b068321c77cac8
 [RNASEQ:INDEX (ggal_1_48850000_49020000)] cache hash: 44b19219efaa91c3d297de5e35a4bbf1
 [RNASEQ:QUANT (ggal_gut)] cache hash: e30f5220838891c1fbc6d95211d51fe4

```

- `resumed_run.log`

```
rnaseq-nf  🍣 master 🅒 base
+  >_ cat ./resumed_run.log | grep 'INFO.*TaskProcessor.*cache hash' | cut -d '-' -f 3 | cut -d ';' -f 1 | sort | tee  ./resumed_run.tasks.l
og
 [MULTIQC] cache hash: 7b05eb097cc35b00610f5ccb66773918
 [RNASEQ:FASTQC (FASTQC on ggal_gut)] cache hash: 2778779a59946ba7b9047ba0fe1bf5a5
 [RNASEQ:INDEX (ggal_1_48850000_49020000)] cache hash: 44b19219efaa91c3d297de5e35a4bbf1
 [RNASEQ:QUANT (ggal_gut)] cache hash: e30f5220838891c1fbc6d95211d51fe4


```

### Inference from top-level process hash log entries

The difference in hashes validates the observed behavior.

*NOTE* The MULTIQC process was rerun since a direct upstream process had to be rerun due to a change in the `script` section.

```diff
-  [MULTIQC] cache hash: 79f6f63ff816a3aecb24ff7a0499b669
+  [MULTIQC] cache hash: 7b05eb097cc35b00610f5ccb66773918

-  [RNASEQ:FASTQC (FASTQC on ggal_gut)] cache hash: 0c18fe122c17bae488b068321c77cac8
+  [RNASEQ:FASTQC (FASTQC on ggal_gut)] cache hash: 2778779a59946ba7b9047ba0fe1bf5a5


```

### Drilling down to individual process

For this process, we inspect the relevant `FASTQC` section *along* with the hashes for it's input, script and container directives.

Here's the `FASTQC` process hashes from the `base_run.log` file

- FASTQC

```
[RNASEQ:FASTQC (FASTQC on ggal_gut)] cache hash: 0c18fe122c17bae488b068321c77cac8; mode: STANDARD; entries: 
  061a0fb0c8f3cecfe11b33689ff21b3d [java.util.UUID] 09247d22-bbf9-4563-9b87-235e2e3dda22 
  195c7faea83c75f2340eb710d8486d2a [java.lang.String] RNASEQ:FASTQC 
  43e5a23fc27129f92a6c010823d8909b [java.lang.String] """
    fastqc.sh "$sample_id" "$reads"
    """
 
  8e58c0cec3bde124d5d932c7f1579395 [java.lang.String] quay.io/nextflow/rnaseq-nf:v1.1 
  7ec7cbd71ff757f5fcdbaa760c9ce6de [java.lang.String] sample_id 
  16b4905b1545252eb7cbfe7b2a20d03d [java.lang.String] ggal_gut 
  553096c532e666fb42214fdf0520fe4a [java.lang.String] reads 
  f60c89cf6646d94fbef280fdfbb7a57c [nextflow.util.ArrayBag] [FileHolder(sourceObj:/home/abhi18av/rnaseq-nf/data/ggal/ggal_gut_1.fq, storePath:/home/abhi18av/rnaseq-nf/data/ggal/ggal_gut_1.fq, stageName:ggal_gut_1.fq), FileHolder(sourceObj:/home/abhi18av/rnaseq-nf/data/ggal/ggal_gut_2.fq, storePath:/home/abhi18av/rnaseq-nf/data/ggal/ggal_gut_2.fq, stageName:ggal_gut_2.fq)] 
  4f9d4b0d22865056c37fb6d9c2a04a67 [java.lang.String] $ 
  16fe7483905cce7a85670e43e4678877 [java.lang.Boolean] true 
  87e9380394e9c557974f5161a6475c95 [sun.nio.fs.UnixPath] /home/abhi18av/rnaseq-nf/bin/fastqc.sh 
```

When we compare these to those in the `resumed_run.log`, we see the following diff

```diff
- [RNASEQ:FASTQC (FASTQC on ggal_gut)] cache hash: 0c18fe122c17bae488b068321c77cac8; mode: STANDARD; entries: 
+ [RNASEQ:FASTQC (FASTQC on ggal_gut)] cache hash: 2778779a59946ba7b9047ba0fe1bf5a5; mode: STANDARD; entries: 

-   43e5a23fc27129f92a6c010823d8909b [java.lang.String] """
-     fastqc.sh "$sample_id" "$reads"
+   ce7882599c370487014c1a2d8fea0e83 [java.lang.String] """
+     fastqc.sh "$sample_id" "$reads" $task.cpus
15a15
+   192fa17dc8adc8300832ec8c9b257b8f [java.util.HashMap$EntrySet] [task.cpus=null] 
18d17

```

**Inferences**

- From the hash entries we can see that the content of the script has changed, `$task.cpus` has been added
- This caused another entry in the `resumed_run.log` regarding the content off the process level directive `cpus` - in this example, `null` as we didn't specify that anywhere in the process or config.

- `MULTIQC`

Now, we drill down the diff for `MULTIQC` process in both logs

```diff
- [MULTIQC] cache hash: 79f6f63ff816a3aecb24ff7a0499b669; mode: STANDARD; entries: 
+ [MULTIQC] cache hash: 7b05eb097cc35b00610f5ccb66773918; mode: STANDARD; entries: 
34c33
-   fac1d0c29da59a72bd617f4ecda550b1 [nextflow.util.ArrayBag] [FileHolder(sourceObj:/home/abhi18av/rnaseq-nf/work/74/08e7523f43026dd7ac80ba7e57df95/ggal_gut, storePath:/home/abhi18av/rnaseq-nf/work/74/08e7523f43026dd7ac80ba7e57df95/ggal_gut, stageName:ggal_gut), FileHolder(sourceObj:/home/abhi18av/rnaseq-nf/work/c0/6d41e0a2d9ee3f0083ef9834a41015/fastqc_ggal_gut_logs, storePath:/home/abhi18av/rnaseq-nf/work/c0/6d41e0a2d9ee3f0083ef9834a41015/fastqc_ggal_gut_logs, stageName:fastqc_ggal_gut_logs)] 
+   e66b49bd490b43c99088a93058c21c0d [nextflow.util.ArrayBag] [FileHolder(sourceObj:/home/abhi18av/rnaseq-nf/work/74/08e7523f43026dd7ac80ba7e57df95/ggal_gut, storePath:/home/abhi18av/rnaseq-nf/work/74/08e7523f43026dd7ac80ba7e57df95/ggal_gut, stageName:ggal_gut), FileHolder(sourceObj:/home/abhi18av/rnaseq-nf/work/1e/687e9b5d894b83e5e0a38b0b8defbb/fastqc_ggal_gut_logs, storePath:/home/abhi18av/rnaseq-nf/work/1e/687e9b5d894b83e5e0a38b0b8defbb/fastqc_ggal_gut_logs, stageName:fastqc_ggal_gut_logs)] 
41d39

```

**Inference**

In this case the process had to be rerun since one of the input files was now being sourced from a different directory, as a result of `FASTQC` being rerun.
