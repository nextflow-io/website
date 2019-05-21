title=A step more towards Nextflow DSL 2
date=2019-05-21
type=post
tags=nextflow,edge,dsl2
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~

The ability to create components libraries or module files has been 
one of the most requested feature ever! 

For this reason I'm very exiting to announce that preview version 
has been merged on master and released in the 
[19.05.0-edge](https://github.com/nextflow-io/nextflow/releases/tag/v19.05.0-edge) version. 

In reality there are some many improvements that this represents a complete
new way to write Nextflow workflows

Since this is still a preview technology and above all to enable to not
break any existing pipelines, to enable this new syntax to need to add 
the following line in your workflow script: 

```
nextflow.preview.dsl=2
```

#### Module files 

A module file consist just of one or more Nextflow process declarations
with the usual syntax. The *only* different is that the `from` and `into` 
clauses in the `input:` and `output:` definitions should not provided 
any more. For example: 

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

A module file can declare optional parameters using the usual Nextflow `params` idiom 
as any other script file.

See a complete module file [here](https://github.com/nextflow-io/rnaseq-nf/blob/66ebeea/modules/rnaseq.nf). 

This makes very easy to migrate existing code to the new modules system. 

### Module includes

A module script can be included from another Nextflow script using the include keyword. Then it's possible to use the included processes composing them each other.

```
nextflow.preview.dsl=2
include 'modules/rnaseq' params(params)

read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )
transcriptome_file = file(params.transcriptome)

INDEX(transcriptome_file)
QUANT(INDEX.out, read_pairs_ch)
```

See the complete example [here](https://github.com/nextflow-io/rnaseq-nf/blob/0bde483/main.nf#L63-L71)

### Channel smart forking 

One of the most exciting change of the new syntax is that any channel 
can be read as many times as you need without the need to duplicate it 
using a `into` operator and making the writing of workflow script much 
fluent and readable. No more channel names proliferation!

```

```

### Nextflow pipes!

The new DSL introduces the `|` (pipe) operator 
which allows to pipe and combine Nextflow processes and operators 
together in a much more fluent and expressive way.  

Consider the following example: 

```
process ampa {

    input:
    file seq

    output:
    file result

    // The BASH script to be executed - for each - sequence
    """
    AMPA.pl -in=${seq} -noplot -rf=result -df=data
    """

}

Channel.fromPath(params.in) | splitFasta | ampa | view { it.text }
```

In the last line the `fromPath` channel is piped to the `splitFast` 
operator which in turn is piped with the `ampa` process and 
finally the output is printed by the `view` operator.





