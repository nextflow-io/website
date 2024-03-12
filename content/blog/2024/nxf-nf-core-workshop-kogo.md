title=Nextflow workshop at the 20th KOGO Winter Symposium
date=2024-03-15
type=post
description=This is the first time we have a Nextflow workshop in Korea, and the feedback was amazing!
image=img/blog-2024-03-15--share.jpg
tags=nextflow,nf-core,workshop
status=published
author=Yuk Kei Wan
icon=yukkei.jpg
~~~~~~

Through a partnership between AWS Asia Pacific and Japan, and Seqera, Nextflow touched ground in South Korea for the first time with a training session at the Korea Genome Organization (KOGO) Winter Symposium. The objective was to introduce participants to Nextflow, empowering them to craft their own pipelines. Recognizing the interest among bioinformaticians, MinSung Cho from AWS Korea’s Healthcare & Research Team decided to sponsor this 90-minute workshop session. This initiative covered my travel expenses and accommodations.

<!-- end-archive-description -->

<div style="margin-top: 2rem; margin-bottom: 2rem;">
    <img src="/img/blog-2024-03-15-kogo-img1a.jpg" alt="Nextflow workshop at KOGO Winter Symposium 2024" />
</div>

The training commenced with an overview of Nextflow pipelines, exemplified by the [nf-core/nanoseq](https://nf-co.re/nanoseq/3.1.0) Nextflow pipeline, highlighting the subworkflows and modules. nfcore/nanoseq is a bioinformatics analysis pipeline for Nanopore DNA/RNA sequencing data that can be used to perform base-calling, demultiplexing, QC, alignment, and downstream analysis. Following this, participants engaged in a hands-on workshop using the AWS Cloud9 environment. In 70 minutes, they constructed a basic pipeline for analyzing nanopore sequencing data, incorporating workflow templates, modules, and subworkflows from [nf-core/tools](https://github.com/nf-core/tools). If you're interested in learning more about the nf-core/nanoseq Nextflow pipeline, I recorded a video talking about it in the nf-core bytesize meeting. You can watch it [here](https://www.youtube.com/watch?v=KM1A0_GD2vQ).

<div style="margin-top: 2rem; margin-bottom: 2rem;">
    <img src="/img/blog-2024-03-15-kogo-img1b.png" alt="Slide from Nextflow workshop at KOGO Winter Symposium 2024" />
</div>

You can find the workshop slides [here](https://docs.google.com/presentation/d/1OC4ccgbrNet4e499ShIT7S6Gm6S0xr38_OauKPa4G88/edit?usp=sharing) and the GitHub repository with source code [here](https://github.com/yuukiiwa/nf-core-koreaworkshop).

The workshop received positive feedback, with participants expressing interest in further sessions to deepen their Nextflow proficiency. Due to this feedback, AWS and the nf-core outreach team are considering organizing small-group local or Zoom training sessions in response to these requests.

It is imperative to acknowledge the invaluable contributions and support from AWS Korea’s Health Care & Research Team, including MinSung Cho, HyunMin Kim, YoungUng Kim, SeungChang Kang, and Jiyoon Hwang, without whom this workshop would not have been possible. Gratitude is also extended to Charlie Lee for fostering collaboration with the nf-core/outreach team.

