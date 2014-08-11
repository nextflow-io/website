title=Share Nextflow pipelines with Github
date=2014-08-07
type=post
tags=git,github,reproducibility
status=published
author=Paolo Di Tommaso
~~~~~~

The [GitHub](https://github.com) code repository and collaboration platform is widely
used between researchers to publish their work and to collaborate on projects source code. 

Even more interestingly a few months ago [GitHub announced improved support for researchers](https://github.com/blog/1840-improving-github-for-science)
making it possible to get a Digital Object Identifier (DOI) for any GitHub repository archive. 

With a DOI for your GitHub repository archive your code becomes formally citable
in scientific publications.

### Why use GitHub with Nextflow

The latest Nextflow release (0.9.0) provides seamless integration with the GitHub platform.
This feature allows you to manage in a more consistent manner your code, or use other people's
Nextflow pipelines published through GitHub in a quick and transparent manner. 


### How it works

The idea is very simple, when you launch a script execution with Nextflow, it will look for
a file with the pipeline name you've specified. If that file does not exist,
it will look for a public repository with the same name on GitHub. If it is found, the 
repository is automatically downloaded to your computer and the code executed. This repository 
is stored in the Nextflow home directory, by default `$HOME/.nextflow`, thus it will be reused
for any further execution.

Having Nextflow installed in your computer (version 0.9.0 or higher), you can try this
feature out by simply entering the following command in your shell terminal:

    nextflow run nextflow-io/hello 
    
    
The first time you execute it, since you don't have pipeline with that name in your computer,
Nextflow will download it from the GitHub repository at this URL `https://github.com/nextflow-io/hello`,
will clone it and execute it, producing the expected output.

The only requirement for a GitHub repository to be managed as a Nextflow project is that it must 
contain at least one file named `main.nf` that defines your Nextflow pipeline script.

### Run a specific revision 

Any Git branch, tag or commit id in the GitHub repository can be used to specify a revision that 
you want to execute when running your pipeline by adding the `-r` option to the run command line. 
So for example you could enter:

    nextflow run nextflow-io-/hello -r mybranch   
  
or 

    nextflow run nextflow-io-/hello -r v1.1
  
  
This can be very useful in order to compare different versions of your project 
and to guarantee consistent results of your pipeline along the evolution of your source code.
  

### Commands to manage pipelines

The following commands provide some basic operations that can be used to manage your pipelines. Anyway 
take in consideration Nextflow is not meant to replace functionalities provided 
by [Git](http://git-scm.com/) tool, you still may need it to create new repository, commit changes, etc.  


#### List available pipelines 

The `ls` command allows you to list all the pipelines you have downloaded in
your computer. For example: 

    $ ./nextflow ls

It prints a list similar to the following one: 

    cbcrg/piper-nf
    nextflow-io/hello

#### Show pipeline information 

By using the `info` command you can show some information of a downloaded pipeline. For example:

    $ nextflow info hello

It prints:    

     repo name  : nextflow-io/hello
     home page  : http://github.com/nextflow-io/hello
     local path : $HOME/.nextflow/assets/nextflow-io/hello
     main script: main.nf
     revisions  : 
     * master (default)
       mybranch
       v1.1 [t]
       v1.2 [t]   
 
Starting from the top it shows: 1) the repository name; 2) the project home page;
3) the local folder where the pipeline has been downloaded; 4) the script that is executed 
when launched; 5) the list of available revisions i.e. branches + tags. Tags are marked with 
a `[t]` of the right, the current revision checked-out is marked with a `*` on the left.
  

#### Pull or updated a pipeline 

The `pull` command allows you to download a pipeline from a GitHub repository or to update 
it if that repository has already been downloaded. For example: 

    nextflow pull nextflow-io/pull
    
Pipeline repositories are stored in the folder `$HOME/.nextflow/assets` in your computer.


#### Clone a pipeline into to a folder 

The `clone` command allows you to copy a Nextflow pipeline project to directory of your choice. For example

    nextflow clone nextflow-io/hello target-dir 
    
If the destination directory is omitted the specified pipeline is cloned to a directory 
with the same name as the pipeline simple name (e.g. `hello`) in the current folder. 

The clone command can be useful to inspect or modify the source code of a pipeline. You can 
eventually commit and push back your changes by using the usual Git/GitHub workflow.  

#### Drop an installed pipeline  

Download pipelines can be deleted by using the `drop` command, as shown below: 

    nextflow drop nextflow-io/hello


### Limitations and know problems 

* GitHub private repository currently are not supported
* Symlinks committed in the Git repository are not resolved correctly 
  when downloaded/cloned by Nextflow 





 


  