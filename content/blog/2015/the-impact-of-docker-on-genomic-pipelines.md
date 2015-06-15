title=The impact of Docker containers on the performance of genomic pipelines
date=2015-06-15
type=post
tags=docker,reproducibility,pipelines,nextflow,genomic
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~

In a recent publication we assessed the impact of Docker containers technology 
on the performance of bioinformatic tools and data analysis workflows. 

We benchmarked three different data analyses: a RNA sequence pipeline for gene expression, 
a consensus assembly and variant calling pipeline, and finally a pipeline for the detection 
and mapping of long non-coding RNAs.  

We found that Docker containers have only a minor impact on the performance 
of common genomic data analysis, which is negligible when the executed tasks are demanding 
in terms of computational time.

*[This publication is available as PeerJ preprint at this link](https://peerj.com/preprints/1171/).*
