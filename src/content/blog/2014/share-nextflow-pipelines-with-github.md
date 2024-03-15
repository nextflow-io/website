---
title: Share Nextflow pipelines with GitHub
date: 2014-08-07
type: post
tags: git,github,reproducibility
author: Paolo Di Tommaso
icon: paolo.jpg
---

The [GitHub](https://github.com) code repository and collaboration platform is widely
used between researchers to publish their work and to collaborate on projects source code.

Even more interestingly a few months ago [GitHub announced improved support for researchers](https://github.com/blog/1840-improving-github-for-science)
making it possible to get a Digital Object Identifier (DOI) for any GitHub repository archive.

With a DOI for your GitHub repository archive your code becomes formally citable
in scientific publications.

### Why use GitHub with Nextflow?

The latest Nextflow release (0.9.0) seamlessly integrates with GitHub.
This feature allows you to manage your code in a more consistent manner, or use other
people's Nextflow pipelines, published through GitHub, in a quick and transparent manner.

### How it works

The idea is very simple, when you launch a script execution with Nextflow, it will look for
a file with the pipeline name you've specified. If that file does not exist,
it will look for a public repository with the same name on GitHub. If it is found, the
repository is automatically downloaded to your computer and the code executed. This repository
is stored in the Nextflow home directory, by default `$HOME/.nextflow`, thus it will be reused
for any further execution.

You can try this feature out, having Nextflow (version 0.9.0 or higher) installed in your computer,
by simply entering the following command in your shell terminal:

    nextflow run nextflow-io/hello

The first time you execute this command Nextflow will download the pipeline
at the following GitHub repository `https://github.com/nextflow-io/hello`,
as you don't already have it in your computer. It will then execute it producing the expected output.

In order for a GitHub repository to be used as a Nextflow project, it must
contain at least one file named `main.nf` that defines your Nextflow pipeline script.

### Run a specific revision

Any Git branch, tag or commit ID in the GitHub repository can be used to specify a revision,
that you want to execute, when running your pipeline by adding the `-r` option to the run command line.
So for example you could enter:

    nextflow run nextflow-io/hello -r mybranch

or

    nextflow run nextflow-io/hello -r v1.1

This can be very useful when comparing different versions of your project.
It also guarantees consistent results in your pipeline as your source code evolves.

### Commands to manage pipelines

The following commands allows you to perform some basic operations that can be used to manage your pipelines.
Anyway Nextflow is not meant to replace functionalities provided by the [Git](http://git-scm.com/) tool,
you may still need it to create new repositories or commit changes, etc.

#### List available pipelines

The `ls` command allows you to list all the pipelines you have downloaded in
your computer. For example:

    nextflow ls

This prints a list similar to the following one:

    cbcrg/piper-nf
    nextflow-io/hello

#### Show pipeline information

By using the `info` command you can show information from a downloaded pipeline. For example:

    $ nextflow info hello

This command prints:

     repo name  : nextflow-io/hello
     home page  : http://github.com/nextflow-io/hello
     local path : $HOME/.nextflow/assets/nextflow-io/hello
     main script: main.nf
     revisions  :
     * master (default)
       mybranch
       v1.1 [t]
       v1.2 [t]

Starting from the top it shows: 1) the repository name; 2) the project home page; 3) the local folder where the pipeline has been downloaded; 4) the script that is executed
when launched; 5) the list of available revisions i.e. branches + tags. Tags are marked with
a `[t]` on the right, the current checked-out revision is marked with a `*` on the left.

#### Pull or update a pipeline

The `pull` command allows you to download a pipeline from a GitHub repository or to update
it if that repository has already been downloaded. For example:

    nextflow pull nextflow-io/examples

Downloaded pipelines are stored in the folder `$HOME/.nextflow/assets` in your computer.

#### Clone a pipeline into a folder

The `clone` command allows you to copy a Nextflow pipeline project to a directory of your choice. For example:

    nextflow clone nextflow-io/hello target-dir

If the destination directory is omitted the specified pipeline is cloned to a directory
with the same name as the pipeline _base_ name (e.g. `hello`) in the current folder.

The clone command can be used to inspect or modify the source code of a pipeline. You can
eventually commit and push back your changes by using the usual Git/GitHub workflow.

#### Drop an installed pipeline

Downloaded pipelines can be deleted by using the `drop` command, as shown below:

    nextflow drop nextflow-io/hello

### Limitations and known problems

- <s>GitHub private repositories currently are not supported</s> Support for private GitHub repositories has been introduced with version 0.10.0.
- <s>Symlinks committed in a Git repository are not resolved correctly
  when downloaded/cloned by Nextflow</s> Symlinks are resolved correctly when using Nextflow version 0.11.0 (or higher).
