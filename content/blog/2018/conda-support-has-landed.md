title=Conda support has landed!
date=2018-06-5
type=post
tags=nextflow,conda,bioconda
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~


Nextflow aims to ease the development of large scale, reproducible, workflows allowing  
developers to focus on the main application logic and to rely on best community tools and 
best practices.

For this reason we are vey excited to announce that the latest Nextflow version (`0.30.0`) finally 
provides built-in support for [Conda](https://conda.io/docs/). 

Conda is a popular package manager that simplifies the installation of software packages 
and the configuration of complex software environments. Above all, it provides access to large 
tool and software package collections maintained by domain specific communities such as 
as [Bioconda](https://bioconda.github.io) and [BioBuild](https://biobuilds.org/). 

The native integration with Nextflow allows researchers to develop workflow application 
in a rapid and easy repeatable manner, reusing community tools, whilst taking advantage of the 
configuration flexibility, portability and scalability provided by Nextflow.


### How it works 

Nextflow  automatically creates and activates the Conda environment(s) given the dependencies
specified by each process.

Dependencies are specified by using the [conda](/docs/latest/process.html#conda) directive, 
providing either the names of the required Conda packages, the path of a Conda environment yaml 
file or the path of an existing Conda environment directory.

Conda environments are stored on the file system. By default Nextflow instructs Conda to save
the required environments in the pipeline work directory. You can specify the directory where the 
Conda environments are stored using the ``conda.cacheDir`` configuration property.


#### Use Conda package names

The simplest way to use one or more Conda packages consists in specifying their names using the ``conda`` directive. 
Multiple package names can be specified by separating them with a blank space. For example:

```
process foo {
    conda "bwa samtools multiqc"

    """
    your_command --here
    """
}
```

Using the above definition a Conda environment that includes BWA, Samtools and MultiQC tools 
is created and activated when the process is executed.

The usual Conda package syntax and naming conventions can be used. The version of a package can be
specified after the package name as shown here ``bwa=0.7.15``.

The name of the channel where a package is located can be specified prefixing the package with
the channel name as shown here ``bioconda::bwa=0.7.15``.

#### Use Conda environment files 

When working in a project requiring a large number of dependencies it can be more convenient 
to consolidate all required tools using a Conda environment file. This is a file that
lists the required packages and channels structured using the YAML format. For example:

```
name: my-env
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - star=2.5.4a
  - bwa=0.7.15
```

The path of the environment file can be specified using the ``conda`` directive:

```
process foo {
  conda '/some/path/my-env.yaml'

  '''
  your_command --here
  '''
}
```

Note: the environment file name **must** end with a ``.yml`` or ``.yaml`` suffix otherwise 
it won't be properly recognized. Also relative paths are resolved against the workflow 
launching directory. 

The suggested approach is to store the the Conda environment file in your project root directory 
and reference it in the `nextflow.config` directory using the ``baseDir`` variable as shown below: 

```
process.conda = "$baseDir/my-env.yaml"
```

This guarantees that the environment paths is correctly resolved independently the execution path. 

See the [documentation](/docs/latest/conda.html) for more details how configure and 
use Conda environments in your Nextflow workflow.  

### Bonus! 

This release includes also a better support for [Biocontainers](https://biocontainers.pro/). So far, 
Nextflow users were able to use container images provided by the Biocontainers community, however
it was not possible to collect process metrics and runtime statistics due the usage of a legacy 
version of the `ps` system tool in those containers that is not compatible with the one 
expected by Nextflow. 

The latest version of Nextflow does not require any more the `ps` tool to fetch execution metrics 
and runtime statics, therefore these information are collected and correctly reported when using Biocontainers 
images.

### Conclusion 

We are very excited by this new feature bringing the ability to use popular Conda tool collections
such as Bioconda directly into Nextflow workflow applications. 

Nextflow developers have now yet another option to transparently manage the dependencies in their
workflows along with [Environment Modules](/docs/latest/process.html#module) and [containers](/docs/latest/docker.html) 
[technology](/docs/latest/singularity.html), giving them great configuration
flexibility. 

The resulting workflow applications can easily be reconfigured and deployed across a range of different 
platforms choosing the best technology according the requirements of the target system.  

