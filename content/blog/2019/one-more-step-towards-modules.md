title=One more step towards Nextflow modules
date=2019-05-22
type=post
tags=nextflow,release,modules,dsl2
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~

The ability to create components, libraries or module files has been 
among the most requested feature ever over the years. 

For this reason, today we are very happy to announce that a preview implementation 
of the [modules feature](https://github.com/nextflow-io/nextflow/issues/984) has been merged 
on master branch of the project and included in the
[19.05.0-edge](https://github.com/nextflow-io/nextflow/releases/tag/v19.05.0-edge) release. 

The implementation of this feature has opened the possibility for many fantastic improvements to Nextflow and its syntax. We are extremely excited as it results in a radical new way of writing Nextflow applications! So much so, that we are referring to these changes as DSL 2.

#### Enabling DSL 2 syntax

Since this is still a preview technology and, above all, to not break 
any existing applications, to enable the new syntax you will need to add 
the following line at the beginning of your workflow script: 

```
nextflow.preview.dsl=2
```

#### Module files 

A module file simply consists of one or more `process` definitions, written with the usual syntax. The *only* difference is that the `from` and `into` clauses in the `input:` and `output:` definition blocks has to be omitted. For example: 

```
process INDEX {
  input:
    file transcriptome 
  output:
    file 'index' 
  script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}
```

The above snippet defines a process component that can be imported in the main 
application script using the `include` statement shown below.  

Also, module files can declare optional parameters using the usual `params` idiom, 
as it can be done in any standard script file.

This approach, which is consistent with the current Nextflow syntax, makes very easy to migrate existing code to the new modules system, reducing it to a mere copy & pasting exercise in most cases. 

You can see a complete module file [here](https://github.com/nextflow-io/rnaseq-nf/blob/66ebeea/modules/rnaseq.nf). 

### Module inclusion

A module file can be included into a Nextflow script using the `include` statement. 
With this it becomes possible to reference any process defined in the module using the usual syntax for a function invocation, and specifying the expected input channels as they were function arguments. 

```
nextflow.preview.dsl=2
include 'modules/rnaseq'

read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )
transcriptome_file = file( params.transcriptome )

INDEX( transcriptome_file )
FASTQC( read_pairs_ch )
QUANT( INDEX.out, read_pairs_ch )
MULTIQC( QUANT.out.mix(FASTQC.out).collect(), multiqc_file )
```

Notably, each process defines its own namespace in the script scope which allows the access of the process output channel(s) using the `.out` attribute. This can be used then as any other Nextflow channel variable in your pipeline script.

The `include` statement gives also the possibility to include only a [specific process](https://www.nextflow.io/docs/edge/dsl2.html#selective-inclusion)
or to include a process with a different [name alias](https://www.nextflow.io/docs/edge/dsl2.html#module-aliases). 

### Channel smart forking 

One of the most important changes of the new syntax is that any channel can be read as many 
times as you need removing the requirement to duplicate them using the `into` operator. 

For example, in the above snippet, the `read_pairs_ch` channel has been used twice, as input both for the `FASTQC` and `QUANT` processes. Nextflow forks it behind the scene for you. 

This makes the writing of workflow scripts much more fluent, readable and ... fun! No more channel names proliferation!


### Nextflow pipes!

Finally, maybe our favourite one. The new DSL introduces the `|` (pipe) operator which allows for the composition 
of Nextflow channels, processes and operators together seamlessly in a much more expressive way.  

Consider the following example: 

```
process align {
  input:
    file seq
  output:
    file 'result'

  """
    t_coffee -in=${seq} -out result
  """
}

Channel.fromPath(params.in) | splitFasta | align | view { it.text }
```

In the last line, the `fromPath` channel is piped to the [`splitFasta`](https://www.nextflow.io/docs/latest/operator.html#splitfasta) operator whose result is used as input by 
the `align` process. Then the output is finally printed by the [`view`](https://www.nextflow.io/docs/latest/operator.html#view) 
operator.

This syntax finally realizes the Nextflow vision of empowering developers to write 
complex data analysis applications with a simple but powerful language that mimics 
the expressiveness of the Unix pipe model but at the same time makes it possible to 
handle complex data structures and patterns as is required for highly 
parallelised and distributed computational workflows.  

#### Conclusion

This wave of improvements brings a radically new experience when it comes to 
writing Nextflow workflows. We are releasing it as a preview technology to allow 
users to try, test, provide their feedback and give us the possibility 
stabilise it. 

We are also working to other important enhancements that will be included soon, 
such as remote modules, sub-workflows composition, simplified file path 
wrangling and more. Stay tuned!
