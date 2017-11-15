title=Running CAW with Singularity
date=2017-11-16
type=post
tags=pipelines,nextflow,genomic,workflow,singularity,cancer
status=published
author=Maxime Garcia
icon=maxime.jpg
~~~~~~
## Running CAW with Singularity

To begin, I'll explain what our workflow is about, how the cluster we're using is set up, why we wanted to switch to Singularity and finally how we managed it.

### Our workflow

![alt text](https://raw.githubusercontent.com/SciLifeLab/CAW/master/doc/images/CAW_logo.png "Cancer Analysis Workflow logo")

First of all, a few word about the [Cancer Analysis Workflow](http://opensource.scilifelab.se/projects/caw/) (CAW for short).
It is developed in collaboration with two infrastructures within [Science for Life Laboratory](https://www.scilifelab.se/): [National Genomics Infrastructure](https://ngisweden.scilifelab.se/) (NGI), in The Stockholm [Genomics Applications Development Facility](https://www.scilifelab.se/facilities/ngi-stockholm/) to be precise and [National Bioinformatics Infrastructure Sweden](https://www.nbis.se/) (NBIS).

CAW is a Nextflow workflow designed to analyze tumor/normal pairs.
Adapted from [GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/) for the preprocessing of fastq files, then various variant calling tools are used to look for somatic SNVs and small indels ([MuTect1](https://github.com/broadinstitute/mutect/), [MuTect2](https://github.com/broadgsa/gatk-protected/), [Strelka](https://github.com/Illumina/strelka/), [Freebayes](https://github.com/ekg/freebayes/)), ([GATK HaplotyeCaller](https://github.com/broadgsa/gatk-protected/)), for structural variants([Manta](https://github.com/Illumina/manta/)) and for CNVs ([ASCAT](https://github.com/Crick-CancerGenomics/ascat/)).
Annotation tools ([snpEff](http://snpeff.sourceforge.net/), [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)) are also used, and finally [MultiQC](http://multiqc.info/) for handling reports.

We are currently working on a manuscript, so no publication yet, but if you want further information, you're welcome to look at (or even contribute to) our [github repository](https://github.com/SciLifeLab/CAW/), or come talk with us on our [gitter channel](https://gitter.im/SciLifeLab/CAW/).
Now that the details have been quickly explained, let's talk about our set up.

### UPPMAX clusters

We're based in Sweden, and [Uppsala Multidisciplinary Center for Advanced Computational Science](https://uppmax.uu.se/) (UPPMAX) provides Computational infrastructures for all Swedish researchers (you do have to apply though).
Since we're analyzing sensitive data, we are using secure clusters (with a two factor authentication), set up by UPPMAX: [SNIC-SENS](https://www.uppmax.uu.se/projects-and-collaborations/snic-sens/).

In my case, since we're still in the development part, I am mainly using the research cluster [Bianca](https://www.uppmax.uu.se/resources/systems/the-bianca-cluster/).
So I can only transfer files and data in one specific repository using SFTP.

UPPMAX is providing computing resources for Swedish researchers for all scientific domains, so getting software updates can occasionally take some times.
But they are using [environment modules](http://modules.sourceforge.net/) which allow several versions of different tools, which is good for reproducibility, and from my own perspective quite easy to use, but there is still some room left for improvement regarding portability.

### Why using containers?

The idea of using containers, for improved portability and reproducibility, and more up to date tools, came naturally to us, as it is easily managed within Nextflow.
Of course we cannot use [Docker](https://www.docker.com/) on our secure cluster, that's why we wanted to run CAW with [Singularity](http://singularity.lbl.gov/) images.

### How was the switch made?

We were already using Docker containers for our CI testing with Travis, and since we use many tools, I took the approach of making (almost) a container for each process.
Because it can take some time, it's repetitive and I~~'m lazy~~ like to automatize everything, I made a simple NF [script](https://github.com/SciLifeLab/CAW/blob/master/buildContainers.nf) to build and push all docker containers.
Basically it's just `build` and `pull` for all containers, with some configuration possibilities.

```groovy
docker build -t ${repository}/${container}:${tag} ${baseDir}/containers/${container}/.

docker push ${repository}/${container}:${tag}
```

Since Singularity can directly pull images from Docker-hub, I slightly modified my building script to pull all containers from Docker-hub to have local Singularity images.

```groovy
singularity pull --name ${container}-${tag}.img docker://${repository}/${container}:${tag}
```

After that it's just a matter of moving all containers to the secure cluster we're using, and using the right configuration file in the profile.
I'll spare you the details of the SFTP transfer.
This is what the configuration file for such Singularity images looks like: [`singularity-path.config`](https://github.com/SciLifeLab/CAW/blob/master/configuration/singularity-path.config)
```groovy
/*
vim: syntax=groovy
-*- mode: groovy;-*-
 * -------------------------------------------------
 * Nextflow config file for CAW project
 * -------------------------------------------------
 * Paths to Singularity images for every process
 * No image will be pulled automatically
 * Need to transfer and set up images before
 * -------------------------------------------------
 */

singularity {
  enabled = true
  runOptions = "--bind /scratch"
}

params {
  containerPath='containers'
  tag='1.2.3'
}

process {
  $ConcatVCF.container      = "${params.containerPath}/caw-${params.tag}.img"
  $RunMultiQC.container     = "${params.containerPath}/multiqc-${params.tag}.img"
  $IndelRealigner.container = "${params.containerPath}/gatk-${params.tag}.img"
  // I'm not putting the whole file here
  // you probably already got the point
}
```

And it almost ran perfectly on the first try, except a process failing due to a typo on a container name...

### Conlusion

This was done a couple of months ago.
We are now using Singularity containers in almost all of our Nextflow pipelines developed at NGI.
Even if we do enjoy the improved control, we must not forgot that:
> With great power comes great responsibility!


### Credits

Thanks to [Rickard Hammarén](https://twitter.com/RickardHammarn) for comments.
