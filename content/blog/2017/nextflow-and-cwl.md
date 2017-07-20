title=Nextflow and the Common Workflow Language
date=2017-07-19
type=post
tags=nextflow,workflow,reproducibility,cwl
status=published
author=Kevin Sayers  
icon=kevin.jpg
~~~~~~

The Common Workflow Language ([CWL](http://www.commonwl.org/)) is a specification for defining 
workflows in a declarative manner. It has been implemented to varying degrees 
by different software packages. Nextflow and CWL share a common goal of enabling portable 
reproducible workflows. 

We are currently investigating the automatic conversion of CWL workflows into Nextflow scripts
to increase the portability of workflows. This work is being developed as 
the [cwl2nxf](https://github.com/nextflow-io/cwl2nxf) project, currently in early prototype stage. 

Our first phase of the project was to determine mappings of CWL to Nextflow and familiarize 
ourselves with how the current implementation of the converter supports a number of CWL specific 
features. 

### Mapping CWL to Nextflow

Inputs in the CWL workflow file are initially parsed as *channels* or other Nextflow input types. 
Each step specified in the workflow is then parsed independently. At the time of writing 
subworkflows are not supported, each step must be a CWL `CommandLineTool` file. 

The image below shows an example of the major components in the CWL files and then post-conversion (click to zoom). 

<a href='/img/cwl2nxf-min.png' target='_blank'>
<img alt='Nextflow cloud deployment' width='760' src='/img/cwl2nxf-min.png' style='margin:1em auto'/>
</a>

CWL and Nextflow share a similar structure of defining inputs and outputs as shown above. 

A notable difference between the two is how tasks are defined. CWL requires either a separate
file for each task or a verbose subworkflow. CWL also requires the explicit mapping of each command
line option for an executed tool. This is done using tags to indicate the position, prefix, etc for
each command line option.

In Nextflow a task command is defined as a separated component in the `process` definition and 
it is ultimately a multiline string which is interpreted by a command script by the underlying 
system. Input parameters can be used in the command string with a simple variable interpolation 
mechanism. This is beneficial as it simplifies porting existing BASH scripts to Nextflow 
with minimal refactoring. 

These examples highlight some of the differences between the two approaches, and the difficulties 
converting complex use cases such as scatter, CWL expressions, and conditional command line inclusion. 


### Current status 

The cwl2nxf is a Groovy based tool with a limited conversion ability. It parses the 
YAML documents and maps the various CWL objects to Nextflow. Conversion examples are 
provided as part of the repository along with documentation for each example specifying the mapping. 

This project was initially focused on developing an understanding of how to translate CWL to Nextflow. 
A number of CWL specific features such as scatter, secondary files and simple JavaScript expressions 
were analyzed and implemented. 


The GitHub repository includes instructions on how to build cwl2nxf and an example usage. 
The tool can be executed as either just a parser printing the converted CWL to stdout, 
or by specifying an output file which will generate the Nextflow script file and if necessary 
a config file. 

The tool takes in a CWL workflow file and the YAML inputs file. It does not currently work 
with a standalone `CommandLineTool`. The following example show how to run it: 

```
java -jar build/libs/cwl2nxf-*.jar rnatoy.cwl samp.yaml
```

See the GitHub [repository](https://github.com/nextflow-io/cwl2nxf) for further details. 

### Conclusion

We are continuing to explore ways to improve the interoperability of Nextflow with CWL. 
Although still an early prototype, the cwl2nxf tool provides some level of conversion of CWL to Nextflow. 


We are also planning to explore [CWL Avro](https://github.com/common-workflow-language/cwlavro), 
which may provide a more efficient way to parse and handle CWL objects for conversion to Nextflow.

Additionally, a number of workflows in the GitHub repository have been implemented in both 
CWL and Nextflow which can be used as a comparison of the two languages. 

The Nextflow team will be presenting a short talk and participating in the Codefest at [BOSC 2017](https://www.open-bio.org/wiki/BOSC_2017). 

We are interested in hearing from the community regarding CWL conversion, and would like 
to encourage anyone interested to contribute to the cwl2nxf project. 
