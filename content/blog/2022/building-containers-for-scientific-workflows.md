title= Building Containers for Scientific Workflows
date=2022-10-20
type=post
description= Follow our step-by-step guide and learn how to build and deploy your own scientific containers with Nextflow.
image=img/inside-scientific-containers.jpg
tags=nextflow,cloud
status=published
author=Evan Floden
icon=evan.jpg
~~~~~~

## About containers and Nextflow
If you have worked in bioinformatics for any time, you have almost certainly come across containers. As a reminder, containers use OS-level process isolation features that enable software to be delivered in lightweight, portable packages across operating environments. While containers have existed in various forms for decades (FreeBSD Jails, Linux VServer, LXC, and Solaris Containers), their usage took off with the release of Docker in 2013. Docker provided a container runtime, a registry, and comprehensive tools for building, managing and deploying containers. In other words, Docker helped make containers mainstream. In addition to Docker, other container runtimes are widely used in scientific computing. These include Singularity, Shifter and Podman between the others.

While a deep knowledge of containers is not required to run Nextflow, having a basic understanding can make pipelines easier to understand and troubleshoot. This article provides an overview of containers and how to use them. It also provides a step-by-step guide to building your own scientific container and explains how to use custom containers in a Nextflow pipeline.

This guide covers the following:

* The case for containers
* Anatomy of a scientific container
* Building your own scientific container
* Using a container in a Nextflow pipeline

If you wish to follow along, you will need a Linux or Mac OS X host with Docker and Nextflow installed. Windows users can also use the Windows Subsystem for Linux (WSL). If you are unfamiliar with how to use Nextflow and Docker with WSL on Windows 10 or 11, see the article [Setting up a Nextflow environment on Windows 10](https://nextflow.io/blog/2021/setup-nextflow-on-windows.html) in the Nextflow blog. OS X users unfamiliar with Docker can reference [How to Use Docker on OS X: The Missing Guide](https://www.viget.com/articles/how-to-use-docker-on-os-x-the-missing-guide/).

## The case for containers
One of the most challenging aspects of ensuring reproducible scientific workflows is setting up and maintaining an application environment. There are literally hundreds of different open-source tools and frameworks in life sciences. Each tool has specific requirements, including pre-requisite software packages, runtimes, and libraries.

Before containers, all the software tools referenced scientific workflows needed to be pre-installed on execution hosts. A typical execution environment was an HPC cluster with a shared file system. You can imagine the complexity of getting all these different tools to co-exist and keeping cluster nodes synchronized. Even minor changes such as applying OS updates, changing middleware, or upgrading a particular tool could cause functionality to break in unpredictable ways. Worse, different pipelines were often tested with different versions of the same tool. This meant that organizations frequently needed to support multiple versions of each tool. 

Containers solve this problem by packaging applications in a fashion that makes them readily portable across compute environments. Rather than having every application installed and working, compute environments only need a container runtime. Applications are deployed to hosts as needed in portable containers. Containers have revolutionized scientific workflows and made application environments easier to deploy and maintain. Containers have also made it much easier to run workflows across diverse computing environments, including various private and public clouds.

## Anatomy of a scientific container
Before we build our own container, looking at a sample container is helpful. The _nextflow/rnaseq-nf_ available on [quay.io](https://quay.io/) (Red Hat’s container registry) is a good example. This container supports the [rnaseq-NF](https://github.com/nextflow-io/rnaseq-nf) proof of concept pipeline published on GitHub. This pipeline has the advantage that it is simple and easy to understand. For knowledgeable readers, a more comprehensive RNA sequencing pipeline ([nf-core/rnaseq](https://nf-co.re/rnaseq)) is available from nf-core.

The nextflow/rnaseq-nf image contains all the applications required by the _RNAseq-NF_ pipeline. This container makes it possible to run the pipeline on any compute environment without worrying about installing applications by encapsulating application logic.

With Docker installed on our local host, we can download (pull) version 1.1 of the container image using the following command:

```
ubuntu@my-host:~$ docker pull quay.io/nextflow/rnaseq-nf:v1.1
v1.1: Pulling from nextflow/rnaseq-nf
33847f680f63: Pull complete
ff810a0db00f: Pull complete
2cb7a358a8ff: Pull complete
977dd9c56199: Pull complete
b6c8d6f41857: Pull complete
0ec885e8834d: Pull complete
Digest: sha256:d6f56ed0eae171fabd324bf582dd5c49c6462662c80a7e69632c57043b6af143
Status: Downloaded newer image for quay.io/nextflow/rnaseq-nf:v1.1
quay.io/nextflow/rnaseq-nf:v1.1
```

After downloading the image, we can run `docker images` to verify that the image is available on our local host: 

```
ubuntu@my-host:~$ docker images
REPOSITORY                   TAG       IMAGE ID       CREATED         SIZE
rancher/rancher              latest    1a0da26e37fa   7 weeks ago     1.49GB
quay.io/nextflow/rnaseq-nf   v1.1      0a20f94e1ea6   8 months ago    2.23GB
hello-world                  latest    feb5d9fea6a5   11 months ago   13.3kB
```

Note that the _nextflow/rnaseq-nf_ container image requires just 2.23GB of disk space. While this is not tiny, it is much smaller than a typical VM image. Because containers are relatively small, they can be loaded and started in seconds.

Depending on your operating system, Docker physically stores downloaded images in `/var/lib/docker`. You can use the `docker info` command to get details about your local Docker settings.

Docker images have a tag that represents the image version number. Versioning container images makes them more maintainable and helps ensure that updates to a container do not accidentally break functionality. Different versions of Nextflow pipelines can be deployed with specific versions of Docker images that are tested and known to work correctly.

## Layers in Docker
When we downloaded the _nextflow/rnaseq-nf_ container above, six separate files were downloaded. This is because Docker images are logical constructs comprised of multiple layers. Each layer corresponds to an instruction in the Dockerfile used to build the image. These layers can be viewed as a series of _diffs_ reflecting what changed from the previous image as new instructions were applied. The final logical image is arrived at by applying each layers in sequence.

You can inspect the layers that comprise a Docker image by running the command below. The jq utility formats the JSON output of the Docker command easier to read.

```
ubuntu@my-host:~$ docker image inspect quay.io/nextflow/rnaseq-nf:v1.1 -f '{{json .RootFS.Layers}}' | jq
[
  "sha256:814bff7343242acfd20a2c841e041dd57c50f0cf844d4abd2329f78b992197f4",
  "sha256:6805e8a21a44da6325009e08c7f5cc9110a9d5259167487ba4c1cd86b14304cd",
  "sha256:afe1551304551dc76e3e5e045674c7a007b42ef9e1c5b98b7efc09976f940ebd",
  "sha256:b6d782ee61cfd1d90eed0c7b7788772582acfc9f4d3d46a92171faa0a0b55c73",
  "sha256:05c13a59978d498fde7fa5805440af469eb7cb73d903c3c1f13c0d17c907a070",
  "sha256:0d696a48560463bce793b6e4264095248e3fd384a08650873bf4e186d6c69da6"
]
ubuntu@ip-172-31-
```

Storing images in layers is helpful because it reduces the bandwidth requirements and the amount of storage required on disk. Suppose an update to a Dockerfile results in a new image version. In that case, Docker only needs to download the changed image layers. Also, new images derived from a common base image can share intermediate layers. These common layers are stored only once on the host operating system. 

We can see how the _nextflow/rnaseq-nf_ container was built using the `docker history` command as shown:

```
ubuntu@my-host:~$ docker history quay.io/nextflow/rnaseq-nf
IMAGE          CREATED       CREATED BY                                      SIZE      COMMENT
7ed5de31bd4d   2 years ago   /bin/sh -c apt-get install -y procps            1.62MB
<missing>      2 years ago   /bin/sh -c conda env update -n root -f conda…   2.1GB
<missing>      2 years ago   /bin/sh -c #(nop) COPY file:ae80bdd367bc9c5f…   153B
<missing>      2 years ago   /bin/sh -c apt-get -y install ttf-dejavu        10.7MB
<missing>      2 years ago   /bin/sh -c #(nop)  MAINTAINER Paolo Di Tomma…   0B
<missing>      2 years ago   /bin/sh -c #(nop)  CMD ["/bin/bash"]            0B
<missing>      2 years ago   /bin/sh -c wget --quiet https://repo.anacond…   131MB
<missing>      2 years ago   /bin/sh -c apt-get update --fix-missing &&  …   210MB
<missing>      2 years ago   /bin/sh -c #(nop)  ENV PATH=/opt/conda/bin:/…   0B
<missing>      2 years ago   /bin/sh -c #(nop)  ENV LANG=C.UTF-8 LC_ALL=C…   0B
<missing>      2 years ago   /bin/sh -c #(nop)  CMD ["bash"]                 0B
<missing>      2 years ago   /bin/sh -c #(nop) ADD file:1901172d265456090…   69.2MB
```

To gain additional insight, we can run the command above with the `--no-trunc` option:

```
ubuntu@my-host:~ $ docker history --no-trunc quay.io/nextflow/rnaseq-nf
```

## Running a container interactively
The _nextflow/rnaseq-nf_ image downloaded above has the following tools pre-installed:

* [Salmon](https://combine-lab.github.io/salmon/) - a tool for quantifying the expression of transcripts using RNA-seq data
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  - a tool that provides an overview of basic quality control metrics for raw next-generation sequencing data
* [MultiQC](https://multiqc.info/)  - a tool to create a single report visualizing output from multiple tools

After pulling the container, we can run it interactively using the `-it` switches on the `docker run` command line. The `-i` switch tells Docker to run the container interactively, and `-t` allocates a tty so that we can interact with the container from the command line.

We start a bash shell in the container and execute a series of commands to verify the applications' versions and ensure they run. We can exit the container by typing `exit` from the shell.

```
ubuntu@my-host:/var/lib$ docker run -it quay.io/nextflow/rnaseq-nf:v1.1 bash
(base) root@9766471494d1:/# which salmon
/opt/conda/bin/salmon
(base) root@8bfb622d25ad:/# salmon -v
salmon 1.0.0
(base) root@8bfb622d25ad:/# fastqc -v
FastQC v0.11.9
(base) root@9766471494d1:/# multiqc –version
multiqc, version 1.11
```

## Building your own scientific container
When building a container, we start with a base image. While opinions vary on which base image to use, it is a good idea to pick a small base image. The less software contained in the image, the faster it will load. Smaller images also have fewer security vulnerabilities because they present a smaller "attack surface." Most container registries provided a recommended set of trusted base images. 

While your containerized application may not be Python based, the following article provides some helpful tips for selecting a base image depending on your requirements - [The best Docker base image for your Python application](https://pythonspeed.com/articles/base-image-python-docker-images/).

For our example, we have chosen _Ubuntu 22.04_ as our base image. This is a Debian-based image maintained by Debian's developers. It is up to date and has a relatively small 78MB footprint.

To ensure the container meets your needs, you may want to download the base image and explore it interactively. This will help you get a feel for what is already installed and what additional software you need to add. You can also experiment by installing software and validating the installation procedure.

```
$ docker pull ubuntu:22.04
$ docker run -it ubuntu:22.04 bash
```
Docker images are constructed using a Dockerfile. A Dockerfile is a text document that contains a set of commands used to assemble an image. The Dockerfile starts with a reference to a base image, and includes a set of commands that will be run in sequence to construct the image.

Create a file called “Dockerfile” as shown using your favorite text editor:

```
FROM ubuntu:22.04
MAINTAINER GJS
RUN apt-get update && apt-get install -y curl
RUN apt-get install -y python3 python3-pip
RUN curl -sSL https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz | tar xz \
&& mv /salmon-*/bin/* /usr/bin/ \
&& mv /salmon-*/lib/* /usr/lib/
```

Docker instructions include FROM, COPY, ADD, USER, LABEL, WORKDIR, and others. A complete reference to Dockerfile instructions is available in the [Dockerfile reference](https://docs.docker.com/engine/reference/builder/) in the Docker documentation.

In our example, we have built a container with just a single life sciences application – salmon. The Dockerfile starts with our base image and applies a series of Dockerfile instructions. 

First, we download an updated list of packages and install curl because curl is not contained in the base image and is required to install salmon. The `-y` switch is required because the apt-get command needs to run without the user answering ‘yes’ to questions presented during the installation.

Next, we download salmon using curl, uncompress and extract it, and move the salmon binary and library files to directories included in our default PATH so that command(s) can be found.

Docker provides [best practices](https://docs.docker.com/develop/develop-images/dockerfile_best-practices/) for building containers using Dockerfile in the Docker documentation. For example, multiple commands are frequently chained together using && as shown in our example, to minimize the number of discrete commands and associated layers in the resulting Docker image.

Next, we can build an image based on the Dockerfile. In this example, we tag our image with the name _my-bio-container_. The `-t` switch is used to specify an image name and optional tag in the format _name:tag_.

```
ubuntu@my-host:~$ docker build -t my-bio-container .
Sending build context to Docker daemon  242.3MB
Step 1/5 : FROM ubuntu:22.04
 ---> df5de72bdb3b
Step 2/5 : MAINTAINER GJS
 ---> Using cache
 ---> 4f85cd946666
Step 3/5 : RUN apt-get update && apt-get install -y curl
 ---> Using cache
 ---> 0b4faff6c7c9
Step 4/5 : RUN apt-get install -y python3 python3-pip
 ---> Using cache
 ---> 91852b732883
Step 5/5 : RUN curl -sSL https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz | tar xz && mv /salmon-*/bin/* /usr/bin/ && mv /salmon-*/lib/* /usr/lib/
 ---> Running in 1f5d7f3b3f48
Removing intermediate container 1f5d7f3b3f48
 ---> 79a95ad1aa8d
Successfully built 79a95ad1aa8d
Successfully tagged my-bio-container:latest
```
The `docker build` command shows instructions in the Dockerfile being applied in sequence to arrive at the _my-bio-container_ image.

We can run _docker images_ to verify that the new container is available and ready to use:

```
ubuntu@my-host:~$ docker images my-bio-container
REPOSITORY         TAG       IMAGE ID       CREATED         SIZE
my-bio-container   latest    79a95ad1aa8d   2 minutes ago   684MB
```

We can see how the container was constructed using the docker history command:

```
ubuntu@my-host:~$ docker history my-bio-container:latest
IMAGE          CREATED              CREATED BY                                      SIZE      COMMENT
79a95ad1aa8d   2 minutes ago	   /bin/sh -c curl -sSL https://github.com/COMB…   219MB
91852b732883   14 minutes ago       /bin/sh -c apt-get install -y python3 python…   344MB
0b4faff6c7c9   14 minutes ago       /bin/sh -c apt-get update && apt-get install…   42.8MB
4f85cd946666   14 minutes ago       /bin/sh -c #(nop)  MAINTAINER GJS               0B
df5de72bdb3b   2 weeks ago          /bin/sh -c #(nop)  CMD ["bash"]                 0B
<missing>      2 weeks ago          /bin/sh -c #(nop) ADD file:396eeb65c8d737180…   77.8MB
```

Container developers may choose to use the `docker-squash` command to collapse some of the intermediate container image layers. Squashing layers has pros and cons. On the one hand, a single large container image may load faster. However, on the other hand, it reduces the opportunity to re-use layers between images making container images less efficient to store.

We can also explore the new container interactively and verify that the python and salmon packages are installed and working:

```
ubuntu@my-host:~ $ docker run -it my-bio-container bash
root@92c93eb5d26e:/# which salmon
/usr/bin/salmon
root@92c93eb5d26e:/# which python3
/usr/bin/python3
root@92c93eb5d26e:/# salmon -v
salmon 1.5.2
```

## Testing a container with application data
To do something meaningful with our new container, it is helpful to have some data to work with. You can obtain sample data by cloning a public data set used for training purposes by Seqera Labs: 

```
ubuntu@my-host:~$ cd ~
ubuntu@my-host:~$ git clone https://github.com/seqeralabs/nf-training-public
Cloning into 'nf-training-public'...
remote: Enumerating objects: 3362, done.
remote: Counting objects: 100% (791/791), done.
remote: Compressing objects: 100% (327/327), done.
remote: Total 3362 (delta 562), reused 664 (delta 464), pack-reused 2571
Receiving objects: 100% (3362/3362), 44.81 MiB | 27.84 MiB/s, done.
Resolving deltas: 100% (2024/2024), done.
```

In the example above, we clone the _nf-training-public_ repository in our home directory. Salmon executes by reading a FASTA format file containing reference transcripts and a set of reads. We can verify that there is a _transcriptome.fa_ file suitable for use with salmon:

```
ubuntu@my-host:~$ ls -al nf-training-public/nf-training/data/ggal/*.fa
-rwxrwxr-x 1 ubuntu ubuntu 173911 Aug 22 17:56 transcriptome.fa
```

The first phase of Salmon’s execution involves creating an index. We can test salmon by running it within the container against the _transcriptome.fa_ file downloaded in the previous step to generate an index:

```
docker run --volume $PWD:$PWD --workdir $PWD my-bio-container:latest \
  salmon index -t $PWD/nf-training-public/nf-training/data/ggal/transcriptome.fa -i transcript-index
```

The command above requires some explanation. First, there are two parameters passed to the docker run command.

* The `--volume $PWD:$PWD` switch causes Docker to mount `/home/ubuntu` on the host to `/home/ubuntu` in the container. The first instance of `$PWD` refers to the path on the host, and the second instance following the colon refers to the path inside the container.
* The switch `--workdir $PWD` makes our home directory the working directory. This means that _transcript-index_ will be written in the container to `/home/ubuntu` which is mapped to the same directory path on the Docker host. We could have specified these details in our Dockerfile using the WORKDIR and VOLUME instructions. 
* The remainder of the command instructs salmon to create an index providing the path to the source transcripts file and the destination index file.

After executing the container, we see the following output:

```ubuntu@my-host:~$ docker run --volume $PWD:$PWD --workdir $PWD my-bio-container:latest \
>   salmon index -t $PWD/nf-training-public/nf-training/data/ggal/transcriptome.fa -i transcript-index
…
[2022-08-22 18:11:28.617] [jLog] [info] building index
out : transcript-index
[2022-08-22 18:11:28.617] [puff::index::jointLog] [info] Running fixFasta
…
[2022-08-22 18:11:28.820] [puff::index::jointLog] [info] finished populating pos vector
[2022-08-22 18:11:28.820] [puff::index::jointLog] [info] writing index components
[2022-08-22 18:11:28.820] [puff::index::jointLog] [info] finished writing dense pufferfish index
[2022-08-22 18:11:28.821] [jLog] [info] done building index
for info, total work write each  : 2.331    total work inram from level 3 : 4.322  total work raw : 25.000
Bitarray          901696  bits (100.00 %)   (array + ranks )
final hash             0  bits (0.00 %) (nb in final hash 0)
ubuntu@my-host:~$
```

The output above is abbreviated. Although salmon ran inside the _my-bio-container_ container, the directory _transcript-index_ is available in our current working directory on the host. We can inspect the _transcript-index_ generated by salmon as follows:

```
ubuntu@my-host:~$ ls -al transcript-index/
total 648
drwxr-xr-x  2 root   root     4096 Aug 22 18:11 .
drwxr-x--- 11 ubuntu ubuntu   4096 Aug 22 18:11 ..
-rw-r--r--  1 root   root       12 Aug 22 18:11 complete_ref_lens.bin
-rw-r--r--  1 root   root      173 Aug 22 18:11 ctable.bin
-rw-r--r--  1 root   root       48 Aug 22 18:11 ctg_offsets.bin
-rw-r--r--  1 root   root       25 Aug 22 18:11 duplicate_clusters.tsv
-rw-r--r--  1 root   root      997 Aug 22 18:11 info.json
-rw-r--r--  1 root   root   113148 Aug 22 18:11 mphf.bin
-rw-r--r--  1 root   root   384656 Aug 22 18:11 pos.bin
-rw-r--r--  1 root   root      496 Aug 22 18:11 pre_indexing.log
-rw-r--r--  1 root   root    21456 Aug 22 18:11 rank.bin
-rw-r--r--  1 root   root       16 Aug 22 18:11 refAccumLengths.bin
-rw-r--r--  1 root   root     2598 Aug 22 18:11 ref_indexing.log
-rw-r--r--  1 root   root       12 Aug 22 18:11 reflengths.bin
-rw-r--r--  1 root   root    42784 Aug 22 18:11 refseq.bin
-rw-r--r--  1 root   root    42872 Aug 22 18:11 seq.bin
-rw-r--r--  1 root   root      126 Aug 22 18:11 versionInfo.json
```

## Publishing a container to a registry
Now that we know our container works, we can publish it to a registry. This will allow others to use the container in their own pipelines.

Several public container registries are available, including Docker Hub, quay.io, Amazon Elastic Container Registry, and others. In our example, we will publish our container to Docker Hub.

If you have not used Docker Hub before, you will need to create a profile on Docker Hub. You will also need to create a repository to hold your containers. In this example, we have created a repository called _<my-org>/nextflow_.

Run the following command from the shell to log in to your Docker Hub account:

```
ubuntu@my-host:~$ docker login
Login with your Docker ID to push and pull images from Docker Hub. If you don't have a Docker ID, head over to https://hub.docker.com to create one.
Username: <my-org>
Password:
Login Succeeded
ubuntu@my-host:~$
```

Tag your image with your Docker Hub username/organization name:

```
ubuntu@my-host:~$ docker tag my-bio-container <my-org>/my-bio-container
```

Push your image to Docker Hub:

```
ubuntu@my-host:~$ docker push <my-org>/my-bio-container
Using default tag: latest
The push refers to repository [docker.io/<my-org>/my-bio-container]
337c27f8f75a: Pushed
605601be3719: Pushed
865b378ebca2: Pushed
629d9dbab5ed: Mounted from library/ubuntu
latest: digest: sha256:5f4cc9e4475ec3b56ab0ca0d191b95d2210e8cd4b04871a1fe88c7bf5cc62e84 size: 1166
ubuntu@my-host:~$
```

The output of the `docker push` command shows the benefit of storing Docker images in layers. Note that Docker Hub recognized the base image. Rather than upload a new copy of it, it simply linked to the existing ubuntu base image.

## Using a container in a Nextflow pipeline
Now that our container has been published to a registry, we can use the container in a Nextflow pipeline.

We can create a simple one-step workflow called _my-pipeline.nf_ that generates an index file based on the transcriptome file in the reference data set. This performs the same steps we used to test _salmon_ running in our container, except we have expressed the logic in Nextflow's DSL.

```
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.outdir = "results"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

workflow {
    index_ch = INDEX(params.transcriptome_file)
}
```

The Nextflow pipeline above does not refer to where _salmon_ can be found or what container technology to use. The pipeline code simply assumes that salmon is available when the process step executes in the target compute environment.

To provide details about the execution environment, we need a companion _nextflow.config_ file that provides these additional details. Create a _nextflow.config_ file in the same directory as your pipeline that contains the following:

```
process.executor = 'local'
docker.enabled = 'true'
process.container = '<my-org>/my-bio-container'
docker.runOptions = '-u $(id -u):$(id -g)'
```

The lines in _nextflow.config_ tell Nextflow the following:

* Task execution will take place on the local host as opposed to a cluster or cloud service.
* Tasks will execute inside a Docker container.
* Process steps will run in the <my-org>my-bio-container image available from Docker Hub.
* Nextflow provides a docker.runOptions setting to pass arguments to the docker run command as specified in the [Docker documentation](https://docs.docker.com/engine/reference/run/). We override the default execution user in the container (normally root) and substitute this with the numeric user-id and group-id of the host operating system user to ensure that files created inside the container are readable when written to the mapped volume.

Nextflow assembles its configuration information from multiple files, and precedence rules apply. It is a good idea to run `nextflow config` from the directory containing your pipeline and _nextflow.config_ file to make sure that your configuration is correct:

```
ubuntu@my-host:~ $ nextflow config
process {
   executor = 'local'
   container = '<my-org>/my-bio-container'
}

docker {
   enabled = 'true'
   runOptions = '-u $(id -u):$(id -g)'
}
```

Now run the pipeline:

```
ubuntu@my-host:~$ nextflow run my-pipeline.nf
N E X T F L O W  ~  version 22.04.5
Launching `my-pipeline.nf` [pensive_nightingale] DSL2 - revision: 114efa731c
R N A S E Q - N F   P I P E L I N E
===================================
transcriptome: /home/ubuntu/nf-training-public/nf-training/data/ggal/transcriptome.fa
outdir       : results

executor >  local (1)
[23/df492a] process > INDEX [100%] 1 of 1 ✔
```

Congratulations! If you have gotten this far, you have successfully executed a pipeline using your own custom-built scientific container.

## Some helpful examples
The _nf-core_ project provides excellent examples of how production containers are built "in the wild." Pipeline developers frequently include the Dockerfiles used to construct containers in their pipeline repos on GitHub. These examples are instructive because _nf-core_ pipelines and containers are developed based on best practices by experienced bioinformaticians. These pipelines and containers tend to be extensively peer-reviewed.

The [nf-core rnafusion](https://github.com/nf-core/rnafusion) pipeline is an excellent example. This project repo contains Dockerfiles for five containers referenced by the pipeline in a [containers](https://github.com/nf-core/rnafusion/tree/master/containers) directory. Bioinformaticians often use _conda_ as a tool to manage software installations. These examples show Dockerfiles and accompanying YAML files that define each container's _conda_ environment.

## Life sciences container registries
In this guide, we built our own container, but thousands of curated containers are available from multiple registries. Before building your own container, there is nothing wrong with leveraging publicly available bioinformatics containers from trusted sources and using them in your workflow. Some registries of curated bioinformatics containers are:

* https://hub.docker.com/u/nfcore  - a list of curated containers used in nf-core pipelines
* https://biocontainers.pro/registry - a community-driven project to create and manage bioinformatics software containers
* https://quay.io – A searchable repository supporting multiple container formats operated by Red Hat
* https://depot.galaxyproject.org/singularity  - an extensive set of life sciences Singularity containers

## More about containers in Nextflow
Containers are treated as "first-class citizens" in Nextflow. While running workflow steps outside of containers is possible, containerized deployments are the norm. Pipelines are written to be independent of underlying compute and container environments. Details about containers and compute environments are maintained in a separate _nextflow.config_ file deployed along with each pipeline. 

As we saw in our example above, the _process.container_ directive in the _nextflow.config_ file tells Nextflow where to find the container(s) in which process steps execute. For example, for the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline hosted on GitHub, a single _process.container_ is specified in the [nextflow.config](https://github.com/nextflow-io/rnaseq-nf/blob/master/nextflow.config#L35) file:

```
process.container = 'quay.io/nextflow/rnaseq-nf:v1.1'
```

In the rnaseq-nf pipeline, all process steps execute in the same _nextflow/rnaseq-nf_ container. If the container is not on the execution host at runtime, Nextflow automatically retrieves it from Red Hat's quay.io repository. 

With complex pipelines, pipeline authors may choose to run each process step in a different container for maintainability. Several of the curated nf-core pipelines take this approach. Using a separate container for each workflow step can make pipelines easier to maintain since process steps can be tested and debugged individually. Developers can be confident that a change to a container used by one step will not affect other steps in the pipeline.

Nextflow provides configuration flexibility. In the example below, all process steps run by default in the _nextflow/rnaseq-nf_ image available from Docker Hub. The FASTQ pipeline step, however, will execute in a separate container listed on Docker Hub provided by _biocontainers_. Pipeline authors can optionally specify different containers for each process step using the `withName:` syntax illustrated below. 

```
process.container = 'nextflow/rnaseq-nf'

process {
    withName:FASTQC {
        container = 'biocontainers/fastqc'
    }
}
```

## Wrapping up
In this article, we have provided an overview of scientific containers and explained how to use them. We have also provided a step-by-step tutorial explaining how to:

* Build your own container
* Publish it to a contain registry
* Use a container in a Nextflow pipeline

Once you get the hang of building and maintaining containers, they are surprisingly easy to use. After a time, they will become so familiar you will wonder how you ever worked without them.  

To learn more, you can inspect the pipelines available at [nf-core](https://nf-co.re/) for examples of production-quality containerized pipelines. These pipelines use prebuilt [BioContainers](https://biocontainers.pro/) which offer a great way to get started with commonally used software packages and tools.

Nevertheless, challenges still remains when managing containers at scale. Modern pipelines can require the use of dozens of containers that can be difficult to manage and kept updated. Stay tuned for a follow up article on how to tackle these challenges.
