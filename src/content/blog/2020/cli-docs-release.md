---
title: The Nextflow CLI - tricks and treats!
date: 2020-10-22
type: post
tags: nextflow,docs
author: Abhinav Sharma
icon: abhinav.jpg
---

For most developers, the command line is synonymous with agility. While tools such as [Nextflow Tower](https://tower.nf) are opening up the ecosystem to a whole new set of users, the Nextflow CLI remains a bedrock for pipeline development. The CLI in Nextflow has been the core interface since the beginning; however, its full functionality was never extensively documented. Today we are excited to release the first iteration of the CLI documentation available on the [Nextflow website](https://www.nextflow.io/docs/edge/cli.html).

And given Halloween is just around the corner, in this blog post we'll take a look at 5 CLI tricks and examples which will make your life easier in designing, executing and debugging data pipelines. We are also giving away 5 limited-edition Nextflow hoodies and sticker packs so you can code in style this Halloween season!


### 1. Invoke a remote pipeline execution with the latest revision

Nextflow facilitates easy collaboration and re-use of existing pipelines in multiple ways. One of the simplest ways to do this is to use the URL of the Git repository.

```
$ nextflow run https://www.github.com/nextflow-io/hello
```

When executing a pipeline using the run command, it first checks to see if it has been previously downloaded in the ~/.nextflow/assets directory, and if so, Nextflow uses this to execute the pipeline. If the pipeline is not already cached, Nextflow will download it, store it in the `$HOME/.nextflow/` directory and then launch the execution.

How can we make sure that we always run the latest code from the remote pipeline? We simply need to add the `-latest` option to the run command, and Nextflow takes care of the rest.

```
$ nextflow run nextflow-io/hello -latest
```

### 2. Query work directories for a specific execution

For every invocation of Nextflow, all the metadata about an execution is stored including task directories, completion status and time etc. We can use the `nextflow log` command to generate a summary of this information for a specific run.

To see a list of work directories associated with a particular execution (for example, `tiny_leavitt`), use:

```
$ nextflow log tiny_leavitt
```

To filter out specific process-level information from the logs of any execution, we simply need to use the fields (-f) option and specify the fields.

```
$ nextflow log tiny_leavitt –f 'process, hash, status, duration'
```

The hash is the name of the work directory where the process was executed; therefore, the location of a process work directory would be something like `work/74/68ff183`.

The log command also has other child options including `-before` and `-after` to help with the chronological inspection of logs.


### 3. Top-level configuration

Nextflow emphasizes customization of pipelines and exposes multiple options to facilitate this. The configuration is applied to multiple Nextflow commands and is therefore a top-level option. In practice, this means specifying configuration options *before* the command.

Nextflow CLI provides two kinds of config overrides - the soft override and the hard override.

The top-level soft override "-c" option allows us to change the previous config in an additive manner, overriding only the fields included the configuration file.

```
$ nextflow -c my.config run nextflow-io/hello
```

On the other hand, the hard override `-C` completely replaces and ignores any additional configurations.

    $ nextflow –C my.config nextflow-io/hello

Moreover, we can also use the config command to inspect the final inferred configuration and view any profiles.

```
$ nextflow config -show-profiles
```

### 4. Passing in an input parameter file

Nextflow is designed to work across both research and production settings. In production especially, specifying multiple parameters for the pipeline on the command line becomes cumbersome. In these cases, environment variables or config files are commonly used which contain all input files, options and metadata. Love them or hate them, YAML and JSON are the standard formats for human and machines, respectively.

The Nextflow run option `-params-file` can be used to pass in a file containing parameters in either format.

```
$ nextflow run nextflow-io/rnaseq -params-file run_42.yaml
```

The YAML file could contain the following.

```
reads      : "s3://gatk-data/run_42/reads/*_R{1,2}_*.fastq.gz"
bwa_index  : "$baseDir/index/*.bwa-index.tar.gz"
paired_end : true
penalty    : 12
```

### 5. Specific workflow entry points

The recently released [DSL2](https://www.nextflow.io/blog/2020/dsl2-is-here.html) adds powerful modularity to Nextflow and enables scripts to contain multiple workflows. By default, the unnamed workflow is assumed to be the main entry point for the script, however, with numerous named workflows, the entry point can be customized by using the `entry` child-option of the run command.

    $ nextflow run main.nf -entry workflow1

This allows users to run a specific sub-workflow or a section of their entire workflow script. For more information, refer to the [implicit workflow](https://www.nextflow.io/docs/latest/dsl2.html#implicit-workflow) section of the documentation.

Additionally, as of version 20.09.1-edge, you can specify the script in a project to run other than `main.nf` using the command line option
`-main-script`.

    $ nextflow run http://github.com/my/pipeline -main-script my-analysis.nf


### Bonus trick! Web dashboard launched from the CLI

The tricks above highlight the functionality of the Nextflow CLI. However, for long-running workflows, monitoring becomes all the more crucial. With Nextflow Tower, we can invoke any Nextflow pipeline execution from the CLI and use the integrated dashboard to follow the workflow execution wherever we are. Sign-in to [Tower](https://tower.nf) using your GitHub credentials, obtain your token from the Getting Started page and export them into your terminal, `~/.bashrc` or include them in your `nextflow.config`.

```
$ export TOWER_ACCESS_TOKEN=my-secret-tower-key
$ export NXF_VER=20.07.1
```

Next simply add the "-with-tower" child-option to any Nextflow run command. A URL with the monitoring dashboard will appear.

```
$ nextflow run nextflow-io/hello -with-tower
```

### Nextflow Giveaway

If you want to look stylish while you put the above tips into practice, or simply like free stuff, we are giving away five of our latest Nextflow hoodie and sticker packs. Retweet or like the Nextflow tweet about this article and we will draw and notify the winners on October 31st!


### About the Author

[Abhinav Sharma](https://www.linkedin.com/in/abhi18av/) is a Bioinformatics Engineer at [Seqera Labs](https://www.seqera.io) interested in Data Science and Cloud Engineering. He enjoys working on all things Genomics, Bioinformatics and Nextflow.


### Acknowledgements

Shout out to [Kevin Sayers](https://github.com/KevinSayers) and [Alexander Peltzer](https://github.com/apeltzer) for their earlier efforts in documenting the CLI and which inspired this work.


*The latest CLI docs can be found in the edge release docs at [https://www.nextflow.io/docs/latest/cli.html](https://www.nextflow.io/docs/latest/cli.html).*
