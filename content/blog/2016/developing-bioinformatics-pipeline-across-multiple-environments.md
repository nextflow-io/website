title=Developing a bioinformatics pipeline across multiple environments
date=2016-02-04
type=post
tags=bioinformatics,reproducibility,pipelines,nextflow,genomic,hpc
status=published
author=Evan Floden
icon=evan.jpg
~~~~~~

As a new bioinformatics student with little formal computer science training, there are 
few things that scare me more than PhD committee meetings and having to run my code in a 
completely different operating environment. 

Recently my work landed me in the middle of the phylogenetic tree jungle and the computational 
requirements of my project far outgrew the resources that were available on our institute’s 
[Univa Grid Engine](https://en.wikipedia.org/wiki/Univa_Grid_Engine) based cluster. Luckily for me, 
an opportunity arose to participate in a joint program at the MareNostrum HPC at the 
[Barcelona Supercomputing Centre](http://www.bsc.es) (BSC).

As one of the top 100 supercomputers in the world, the [MareNostrum III](https://www.bsc.es/marenostrum-support-services) 
dwarfs our cluster and consists of nearly 50'000 processors. However it soon became apparent 
that with great power comes great responsibility and in the case of the BSC, great restrictions. 
These include no internet access, restrictive wall times for jobs, longer queues, 
fewer pre-installed binaries and an older version of bash. Faced with the possibility of 
having to rewrite my 16 bodged scripts for another queuing system I turned to Nextflow.

Straight off the bat I was able to reduce all my previous scripts to a single Nextflow script. 
Admittedly, the original code was not great, but the data processing model made me feel confident 
in what I was doing and I was able to reduce the volume of code to 25% of its initial amount 
whilst making huge improvements in the readability. The real benefits however came from the portability.

I was able to write the project on my laptop (Macbook Air), continuously test it on my local 
desktop machine (Linux) and then perform more realistic heavy lifting runs on the cluster, 
all managed from a single GitHub repository. The BSC uses the [Load Sharing Facility](https://en.wikipedia.org/wiki/Platform_LSF) 
(LSF) platform with longer queue times, but a large number of CPUs. My project on the other 
hand had datasets that require over 100'000 tasks, but the tasks processes themselves run 
for a matter of seconds or minutes. We were able to marry these two competing interests 
deploying Nextflow in a [distributed execution manner that resemble the one of an MPI application](/blog/2015/mpi-like-execution-with-nextflow.html).

In this configuration, the queuing system allocates the Nextflow requested resources and 
using the embedded [Apache Ignite](https://ignite.apache.org/) clustering engine, Nextflow handles 
the submission of processes to the individual nodes. 

Here is some examples of how to run the same Nextflow project over multiple platforms.

#### Local

If I wished to launch a job locally I can run it with the command:

    nextflow run myproject.nf
    
#### Univa Grid Engine (UGE)

For the UGE I simply needed to specify the following in the `nextflow.config` file:  

    process {
            executor='uge'
            queue='my_queue'
    }  
    
    
And then launch the pipeline execution as we did before:
    
    nextflow run myproject.nf     
    
    
#### Load Sharing Facility (LSF)

For running the same pipeline in the MareNostrum HPC enviroment, taking advantage of the MPI 
standard to deploy my workload, I first created a wrapper script (for example `bsc-wrapper.sh`) 
declaring the resources that I want to reserve for the pipeline execution:

    #!/bin/bash
    #BSUB -oo logs/output_%J.out
    #BSUB -eo logs/output_%J.err
    #BSUB -J myProject
    #BSUB -q bsc_ls
    #BSUB -W 2:00
    #BSUB -x
    #BSUB -n 512
    #BSUB -R "span[ptile=16]"
    export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
    mpirun --pernode bin/nextflow run concMSA.nf -with-mpi


And then can execute it using `bsub` as shown below:

    bsub < bsc-wrapper.sh

By running Nextflow in this way and given the wrapper above, a single `bsub` job will run 
on 512 cores in 32 computing nodes (512/16 = 32) with a maximum wall time of 2 hours. 
Thousands of Nextflow processes can be spawned during this and the execution can be monitored 
in the standard manner from a single Nextflow output and error files. If any errors occur 
the execution can of course to continued with [`-resume` command line option](/docs/latest/getstarted.html?highlight=resume#modify-and-resume).

### Conclusion 

Nextflow provides a simplified way to develop across multiple platforms and removes 
much of the overhead associated with running niche, user developed pipelines in an HPC 
environment.    