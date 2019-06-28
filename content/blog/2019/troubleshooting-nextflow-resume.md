title=Troubleshooting Nextflow resume
date=2019-07-01
type=post
tags=nextflow,resume
status=published
author=Evan Floden
icon=evan.jpg
~~~~~~

*This two-part blog aims to help users understand Nextflow’s powerful caching mechanism. Part one describes how it works whilst part two will focus on execution provenance and troubleshooting. You can read part one [here](blog/2019/demystifying-nextflow-resume.html)*

#### Troubleshooting resume
If your workflow execution is not resumed as expected, there exists several strategies to debug the problem.

* **Input file changed:** Make sure that there has been no change in your input files. Don’t forget the unique task hash is computed by taking into account the complete file path, the last modified timestamp and the file size. If any of these change, the workflow will be re-executed,  even if the input content is the same. 

* **A process modifies an input:** A process should never alter input files. When this happens, the future execution of tasks will be invalidated for the same reason explained in the previous point.

* **Inconsistent file attributes:** Some shared file system, such as NFS, may report inconsistent file timestamp i.e. a different timestamp for the same file even if it has not been modified. There is an option to use the [lenient mode of caching](https://www.nextflow.io/docs/latest/process.html#cache) to avoid this problem.

* **A race condition in a global variable:** Nextflow is designed to simplify parallel programming without taking care of race conditions and the access of shared resources. One of the few cases in which a race condition may arise is when using a global variable with two (or more) operators. For example:

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

* **Non-deterministic input channels:** While dataflow channel ordering is guaranteed i.e. data is read in the same order in which it’s written in the channel, when a process declares as input two or more channels, each of which is the output of a different process, the overall input ordering is not consistent across different executions.

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

#### Execution provenance

When provided with a run name or session ID, the log command can return useful information about a pipeline execution. This can be composed to track the provenance of a workflow result.

When supplying a run name or session ID, the log command lists all the work directories used to compute the final result. For example:

```
$ nextflow log tiny_fermat

/data/.../work/7b/3753ff13b1fa5348d2d9b6f512153a
/data/.../work/c1/56a36d8f498c99ac6cba31e85b3e0c
/data/.../work/f7/659c65ef60582d9713252bcfbcc310
/data/.../work/82/ba67e3175bd9e6479d4310e5a92f99
/data/.../work/e5/2816b9d4e7b402bfdd6597c2c2403d
/data/.../work/3b/3485d00b0115f89e4c202eacf82eba
```

Using the option `-f` (fields) it’s possible to specify which metadata should be printed by the log command. For example:

```
$ nextflow log tiny_fermat -f 'process,exit,hash,duration'

index	0	7b/3753ff	2s
fastqc	0	c1/56a36d	9.3s
fastqc	0	f7/659c65	9.1s
quant	0	82/ba67e3	2.7s
quant	0	e5/2816b9	3.2s
multiqc	0	3b/3485d0	6.3s
```

The complete list of available fields can be retrieved with the command:

```
$ nextflow log -l
```

The option `-F` allows the specification of filtering criteria to print only a subset of tasks. For example:

```
$ nextflow log tiny_fermat -F 'process =~ /fastqc/'

/data/.../work/c1/56a36d8f498c99ac6cba31e85b3e0c
/data/.../work/f7/659c65ef60582d9713252bcfbcc310
```

This can be useful to locate specific tasks work directories.

Finally, the `-t` option allows for the creation of a basic custom HTML provenance report that can be generated by providing a template file, in any format of your choice. For example:

```
<div>
<h2>${name}</h2>
<div>
Script:
<pre>${script}</pre>
</div>

<ul>
    <li>Exit: ${exit}</li>
    <li>Status: ${status}</li>
    <li>Work dir: ${workdir}</li>
    <li>Container: ${container}</li>
</ul>
</div>
```

By saving the above snippet in a file named template.html, you can run the following command:

```
$ nextflow log tiny_fermat -t template.html > provenance.html
```

#### Resume by default?
Given the majority of users always apply resume, we recently discussed having resume applied by the default. 

Is there any situation where you do not use resume? Would a flag specifying ‘-no-cache’ be enough to satisfy these use cases? 

We want to hear your thoughts on this. Help steer Nextflow development and vote on the twitter poll.