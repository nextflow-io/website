title=Evolution of the Nextflow runtime
date=2022-03-24
type=post
description=Software development is a constantly evolving, this post will summarise the major changes in the evolution of the framework over the next 12 to 18 months.
image=img/nextflow-moving-to-slack.jpg
tags=nextflow,dsl2
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~

Software development is a constantly evolving process that requires continuous adaptation to keep pace with new technologies, user needs, and trends. Likewise, changes are needed in order to introduce new capabilities and guarantee a sustainable development process.

Nextflow is no exception. This post will summarise the major changes in the evolution of the framework over the next 12 to 18 months. 

### Java baseline version 

Nextflow runs on top of Java (or, more precisely, the Java virtual machine). So far, Java 8 has been the minimal version required to run Nextflow. However, this version was released 8 years ago and is going to reach its end-of-life status at the end of [this month](https://endoflife.date/java). For this reason, as of version 22.01.x-edge and the upcoming stable release 22.04.0, Nextflow will require Java version 11 or later for its execution. This also allows the introduction of new capabilities provided by the modern Java runtime. 

Tip: If you are confused about how to install or upgrade Java on your computer, consider using [Sdkman](https://sdkman.io/). It’s a one-liner install tool that allows easy management of Java versions.  

### DSL2 as default syntax 

Nextflow DSL2 has been introduced nearly [2 years ago](https://www.nextflow.io/blog/2020/dsl2-is-here.html) (how time flies!) and definitely represented a major milestone for the project. Established pipeline collections such as those in [nf-core](https://nf-co.re/pipelines) have migrated their pipelines to DSL2 syntax. 

This is a confirmation that the DSL2 syntax represents a natural evolution for the project and is not considered to be just an experimental or alternative syntax. 

For this reason, as for Nextflow version 22.03.0-edge and the upcoming 22.04.0 stable release, DSL2 syntax is going to be the **default** syntax version used by Nextflow, if not otherwise specified. 

In practical terms, this means it will no longer be necessary to add the declaration `nextflow.enable.dsl = 2` at the top of your script or use the command line option `-dsl2 ` to enable the use of this syntax. 

If you still want to continue to use DSL1 for your pipeline scripts, you will need to add the declaration  `nextflow.enable.dsl = 1` at the top of your pipeline script or use the command line option `-dsl1`. 

To make this transition as smooth as possible, we have also added the possibility to declare the DSL version in the Nextflow configuration file, using the same syntax shown above. 

Finally, if you wish to keep the current DSL behaviour and not make any changes in your pipeline scripts, the following variable can be defined in your system environment: 

```
export NXF_DEFAULT_DSL=1 
```

### DSL1 end-of-life phase  

Maintaining two separate DSL implementations in the programming environment is not sustainable and above all, does not make much sense. For this reason, along with making DSL2 the default Nextflow syntax, DSL1 will enter into a 12-month end-of-life phase, at the end of which it will be removed. Therefore version 22.04.x and 22.10.x will be the last stable versions providing the ability to run DSL1 scripts. 

This is required to keep evolving the framework and to create a more solid implementation of Nextflow grammar. Maintaining compatibility with the legacy syntax implementation and data structures makes this a very difficult task.

Bear in mind, this does **not** mean it will not be possible to use DSL1 starting from 2023. All existing Nextflow runtimes will continue to be available, and it will be possible to for any legacy pipeline to run using the required version available from the GitHub [releases page](https://github.com/nextflow-io/nextflow/releases), or by specifying the version using the NXF_VER variable, e.g. 

```
NXF_VER=21.10.6 nextflow run <my/dsl1/pipeline>     
```


### New configuration format

The configuration file is a key component of the Nextflow framework since it allows workflow developers to decouple the pipeline logic from the execution parameters and infrastructure deployment settings.

The current Nextflow configuration file mechanism is extremely powerful, but it also has some serious drawbacks due to its *dynamic* nature that makes it very hard to keep stable and maintainable over time.

For this reason, we are planning to re-engineer the current configuration component and replace it with a better configuration component with two major goals: 1) continue to provide a rich and human-readable configuration system (so, no YAML or JSON), 2) have a well-defined syntax with a solid foundation that guarantees predictable configurations, simpler troubleshooting and more sustainable maintenance. 

Currently, the most likely options are [Hashicorp HCL](https://github.com/hashicorp/hcl) (as used by Terraform and other Hashicorp tools) and [Lightbend HOCON](https://github.com/lightbend/config). You can read more about this feature at [this link](https://github.com/nextflow-io/nextflow/issues/2723).


### Ignite executor deprecation 

The executor for Apache Ignite was an early attempt to provide Nextflow with a self-contained, distributed cluster for the deployment of pipelines into HPC environments. However, it had very little adoption over the years, which was not balanced by the increasing complexity of its maintenance. 

For this reason, it was decided to deprecate it and remove it from the default Nextflow distribution. The module is still available in the form of a separate project plugin and available at [this link](https://github.com/nextflow-io/nf-ignite), however, it will not be actively maintained. 

 
### Conclusion 

This post is focused on the most fundamental changes we are planning to make in the following months. 

With the adoption of Java 11, the full migration of DSL1 to DSL2 and the re-engineering of the configuration system, our purpose is to consolidate the Nextflow technology and lay the foundation for all the new exciting developments and features on which we are working on. Stay tuned for future blogs about each of them in upcoming posts.
