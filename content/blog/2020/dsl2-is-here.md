title=Nextflow DSL 2 is here!
date=2020-07-24
type=post
tags=nextflow,release,modules,dsl2
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~

We are thrilled to announce the stable release of Nextflow DSL 2 as part of the latest 20.07.1 version! 

Nextflow DSL 2 represents a major evolution of the Nextflow language and makes it possible to scale and modularise your data analysis pipeline while continuing to use the Dataflow programming paradigm that characterises the Nextflow processing model.  

We spent more than one year collecting user feedback and making sure that DSL 2 would naturally fit the programming experience Nextflow developers are used to. 


#### DLS 2 in a nutshell

Backward compatibility is a paramount value, for this reason the changes introduced in the syntax have been minimal and above all, guarantee the support of all existing applications. DSL 2 will be an opt-in feature for at least the next 12 to 18 months. After this transitory period, we plan to make it the default Nextflow execution mode. 

As of today, to use DSL 2 in your Nextflow pipeline, you are required to use the following declaration at the top of your script: 

```
nextflow.enable.dsl=2
```

Note that the previous `nextflow.preview` directive is still available, however, when using the above declaration the use of the final syntax is enforced.

#### Nextflow modules

A module file is nothing more than a Nextflow script containing one or more `process` definitions that can be imported from another Nextflow script. 

The only difference when compared with legacy syntax is that the process is not bound with specific input and output channels, as was previously required using the `from` and `into` keywords respectively. Consider this example of the new syntax:

```
process INDEX {
  input:
    path transcriptome 
  output:
    path 'index' 
  script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}
```

This allows the definition of workflow processes that can be included from any other script and invoked as a custom function within the new `workflow` scope. This effectively allows for the composition of the pipeline logic and enables reuse of workflow components. We anticipate this to improve both the speed that users can develop new pipelines, and the robustness of these pipelines through the use of validated modules. 

Any process input can be provided as a function argument using the usual channel semantics familiar to Nextflow developers. Moreover process outputs can either be assigned to a variable or accessed using the implicit `.out` attribute in the scope implicitly defined by the process name itself. See the example below:

```
include { INDEX; FASTQC; QUANT; MULTIQC } from './some/module/script.nf' 

read_pairs_ch = channel.fromFilePairs( params.reads)

workflow {
  INDEX( params.transcriptome )
  FASTQC( read_pairs_ch )
  QUANT( INDEX.out, read_pairs_ch )
  MULTIQC( QUANT.out.mix(FASTQC.out).collect(), multiqc_file )
}
```


Also enhanced is the ability to use channels as inputs multiple times without the need to duplicate them (previously done with the special into operator) which makes the resulting pipeline code more concise, fluent and therefore readable!


#### Sub-workflows

Notably, the DSL 2 syntax allows for the definition of reusable processes as well as sub-workflow libraries. The only requirement is to provide a `workflow` name that will be used to reference and declare the corresponding inputs and outputs using the new `take` and `emit` keywords. For example:

```
workflow RNASEQ {
  take:
    transcriptome
    read_pairs_ch
 
  main: 
    INDEX(transcriptome)
    FASTQC(read_pairs_ch)
    QUANT(INDEX.out, read_pairs_ch)

  emit: 
     QUANT.out.mix(FASTQC.out).collect()
}
```

Now named sub-workflows can be used in the same way as processes, allowing you to easily include and reuse multi-step workflows as part of larger workflows. Find more details [here](/docs/latest/dsl2.html).

#### More syntax sugar


Another exciting feature of Nextflow DSL 2 is the ability to compose built-in operators, pipeline processes and sub-workflows with the pipe (|) operator! For example the last line in the above example could be written as: 

```
emit:
 QUANT.out | mix(FASTQC.out) | collect 
```

This syntax finally realizes the Nextflow vision of empowering developers to write complex data analysis applications with a simple but powerful language that mimics the expressiveness of the Unix pipe model but at the same time makes it possible to handle complex data structures and patterns as is required for highly parallelised and distributed computational workflows.

Another change is the introduction of `channel` as an alternative name as a synonym of `Channel` type identifier and therefore allows the use of `channel.fromPath` instead of `Channel.fromPath` and so on. This is a small syntax sugar to keep the capitazionation  consistent with the rest of the language.

Moreover, several process inputs and outputs syntax shortcuts were removed when using the final version of DSL 2 to make it more predictable. For example, with DSL1, in a tuple input or output declaration the component type could be omitted, for example: 

```
input: 
  tuple foo, 'bar'
```

The `foo` identifier was implicitly considered an input value declaration instead the string `'bar'` was considered a shortcut for `file('bar')`. However, this was a bit confusing especially for new users and therefore using DSL 2, the fully qualified version must be used: 

```
input: 
  tuple val(foo), path('bar')
```

You can find more detailed migration notes at [this link](/docs/latest/dsl2.html#dsl2-migration-notes).  


#### What's next

As always, reaching an important project milestone can be viewed as a major success, but at the same time the starting point for challenges and developments. Having a modularization mechanism opens new needs and possibilities. The first one of which will be focused on the ability to test and validate process modules independently using a unit-testing style approach. This will definitely help to make the resulting pipelines more resilient. 

Another important area for the development of the Nextflow language will be the ability to better formalise pipeline inputs and outputs and further decouple for the process declaration. Nextflow currently strongly relies on the `publishDir` constructor for the generation of the workflow outputs. 

However in the new *module* world, this approach results in `publishDir` being tied to a single process definition. The plan is instead to extend this concept in a more general and abstract manner, so that it will be possible to capture and redirect the result of any process and sub-workflow based on semantic annotations instead of hardcoding it at the task level.

### Conclusion 

We are extremely excited about today's release. This was a long awaited advancement and therefore we are very happy to make it available for general availability to all Nextflow users. We greatly appreciate all of the community feedback and ideas over the past year which have shaped DSL 2. 

We are confident this represents a big step forward for the project and will enable the writing of a more scalable and complex data analysis pipeline and above all, a more enjoyable experience. 
