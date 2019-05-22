title=A step more towards Nextflow DSL 2
date=2019-05-21
type=post
tags=nextflow,edge,dsl2
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~

The ability to create components libraries or module files has been 
one of the most requested feature ever over the years! 

For this reason, I'm very happy to announce that finally a  preview implementation has been merged on master branch of the  project and included in the
[19.05.0-edge](https://github.com/nextflow-io/nextflow/releases/tag/v19.05.0-edge) release. 

The module implementation opened the possibility for some many improvements to Nextflow and its syntax that results in a radical new way of writing Nextflow application about which we are extremely excited! For this reason we are referring to these changes as DLS 2.

#### Enabling DSL 2 syntax

Since this is still a preview technology and, above all, to not break 
any existing application, to enable the new syntax you will need to add 
the following line at the beginning of your workflow script: 

```
nextflow.preview.dsl=2
```

#### Module files 

A module file consists just of one or more Nextflow `process` definitions written with the usual syntax. The *only* difference is that the `from` and `into` clauses in the `input:` and `output:` definition blocks should not be provided any more. For example: 

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

Also, module files can declare optional parameters using the usual Nextflow `params` idiom 
as any other script file.

This approach, which is consistent with the current Nextflow syntax, makes very easy to migrate existing code to the new modules system, reducing it to a mere copy & pasting activity 
most of the time. 

See a complete module file [here](https://github.com/nextflow-io/rnaseq-nf/blob/66ebeea/modules/rnaseq.nf). 

### Module inclusion

A module script can be included from another Nextflow script using the `include` statement. Then it's possible to reference any process defined by the module using the usual syntax for 
a function invocation and specifying the expected input channels as the function arguments. 

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

Notably, each process defines its own namespace in the script scope that 
allows the access of the process output channel(s) using the `.out` attribute 
and that can be used as any other Nextflow channel variable. 

See the complete example [here](https://github.com/nextflow-io/rnaseq-nf/blob/0bde483/main.nf#L63-L71)

The `include` statement allows also to include only a [specific process](https://www.nextflow.io/docs/edge/dsl2.html#selective-inclusion)
or to include with a different [name alias](https://www.nextflow.io/docs/edge/dsl2.html#module-aliases). 

### Channel smart forking 

One of the most important changes of the new syntax is that any channel can be read as many times as you need without the need to duplicate it using an `into` operator. 

For example, in the above snippet, the `read_pairs_ch` channel has been used twice, as input both for the `FASTQC` and `QUANT` processes. Nextflow forks it behind the scene for you. 

This makes the writing of workflow scripts much more fluent, readable and ... fun! No more channel names proliferation!


### Nextflow pipes!

Finally, maybe my favourite one. The new DSL introduces the `|` (pipe) operator which allows composing Nextflow channels, processes and operators together seamlessly in a much more fluent and expressive way.  

Consider the following example: 

```
process align {
  input:
    file seq
  output:
    file 'result'

  """
    t_coffee -in=${seq} > 
  """
}

Channel.fromPath(params.in) | splitFasta | align | view { it.text }
```

In the last line the `fromPath` channel is piped to the `splitFasta` 
operator which result is used as input by the `align` process, finally 
the output is printed by the `view` operator.

This syntax finally realizes the Nextflow vision to empower developers  to write complex data analysis applications with a simple but powerful syntax that mimics the expressiveness of the Unix pipe model but at the same time make it possible to handle complex data structures and patterns as required to carry out highly parallelised distributed computational workflows.   

#### Conclusion

This wave of improvements brings a radically new experience when it comes to writing Nextlfow data analysis applications. We realising it as a preview technology to allow users to try and test and stabilise it. We are also working to other important enhancements that will be included soon, such as remote modules, sub-workflows composition, simplified file paths wrangling and more. Stay tuned!
