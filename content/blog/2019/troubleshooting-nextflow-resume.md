title=Troubleshooting Nextflow resume
date=2019-07-01
type=post
tags=nextflow,resume
status=published
author=Evan Floden
icon=evan.jpg
~~~~~~

*This two-part blog aims to help users understand Nextflow’s powerful caching mechanism. Part one describes how it works whilst part two will focus on execution provenance and troubleshooting. You can read part one [here](blog/2019/demystifying-nextflow-resume.html)*

### Troubleshooting resume

If your workflow execution is not resumed as expected, there exists several strategies to debug the problem.

#### Modified input file(s)

Make sure that there has been no change in your input files. Don’t forget the unique task hash is computed by taking into account the complete file path, the last modified timestamp and the file size. If any of these change, the workflow will be re-executed,  even if the input content is the same. 

#### A process modifying one or more inputs

A process should never alter input files. When this happens, the future execution of tasks will be invalidated for the same reason explained in the previous point.

#### Inconsistent input file attributes

Some shared file system, such as NFS, may report inconsistent file timestamp i.e. a different timestamp for the same file even if it has not been modified. There is an option to use the [lenient mode of caching](https://www.nextflow.io/docs/latest/process.html#cache) to avoid this problem.

#### Race condition in a global variable

Nextflow does its best to simplify parallel programming and to prevent race conditions and the access of shared resources. One of the few cases in which a race condition may arise is when using a global variable with two (or more) operators. For example:

```
Channel
    .from(1,2,3)
    .map { it -> X=it; X+=2 }
    .println { "ch1 = $it" }

Channel
    .from(1,2,3)
    .map { it -> X=it; X*=2 }
    .println { "ch2 = $it" }
```

The problem with this snippet is that the X variable in the closure definition is defined in the global scope. Since operators are executed in parallel, the X value can, therefore, be overwritten by the other map invocation.

The correct implementation requires the use of the def keyword to declare the variable local.

```
Channel
    .from(1,2,3)
    .map { it -> def X=it; X+=2 }
    .println { "ch1 = $it" }

Channel
    .from(1,2,3)
    .map { it -> def X=it; X*=2 }
    .println { "ch2 = $it" }
```

#### Non-deterministic input channels

While dataflow channel ordering is guaranteed i.e. data is read in the same order in which it’s written in the channel, when a process declares as input two or more channels, each of which is the output of a different process, the overall input ordering is not consistent across different executions.

Consider the following snippet:

```
process foo {
  input: set val(pair), file(reads) from reads_ch
  output: set val(pair), file('*.bam') into bam_ch
  """
  your_command --here
  """
}

process bar {
  input: set val(pair), file(reads) from reads_ch
  output: set val(pair), file('*.bai') into bai_ch
  """
  other_command --here
  """
}

process gather {
  input:
  set val(pair), file(bam) from bam_ch
  set val(pair), file(bai) from bai_ch
  """
  merge_command $bam $bai
  """
}
```

The inputs declared in the gather process can be delivered in any order as the execution order of the process foo and bar is not deterministic due to parallel executions.

Therefore, the input of the third process needs to be synchronized using the join operator or a similar approach. The third process should be written as:


```
process gather {
  input:
  set val(pair), file(bam), file(bai) from bam_ch.join(bai_ch)
  """
  merge_command $bam $bai
  """
}
```

#### Still in trouble?

These are most frequent causes of problems with the Nextflow resume mechanism. If you are still not able to resolve 
your problem, identify the first process not resuming correctly, then run your script twice using `-dump-hashes`. You can then compare the resulting `.nextflow.log` files (the first will be named `.nextflow.log.1`). 

Unfortunately, the information reported by `-dump-hashes` can be quite cryptic, however, with the help of a good diff tool it is possible to compare the two log files to identify the reason for the cache to be invalidated.  

#### The golden rule

Never try to debug this kind of problem with production data! This issue can be annoying, but when it happens
it ahould be able to be replicated in a consistent manner with any data.

Therefore, we always suggest Nextflow developers include in their pipeline project 
a small synthetic dataset to easily execute and test the complete pipeline execution in a few seconds. 
This is the golden rule for debugginh and troubleshooting execution problems avoids getting stuck with production data.

#### Resume by default?
Given the majority of users always apply resume, we recently discussed having resume applied by the default. 

Is there any situation where you do not use resume? Would a flag specifying `-no-cache` be enough to satisfy these use cases? 

We want to hear your thoughts on this. Help steer Nextflow development and vote in the twitter poll below.

<blockquote class="twitter-tweet" data-partner="tweetdeck"><p lang="en" dir="ltr">Should -resume⏯️ be the default when launching a Nextflow pipeline?</p>&mdash; Nextflow (@nextflowio) <a href="https://twitter.com/nextflowio/status/1145599932268785665?ref_src=twsrc%5Etfw">July 1, 2019</a></blockquote>
<script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

*In the following post of this series, we will show how to produce a provenance report using a built-in Nextflow command.*
