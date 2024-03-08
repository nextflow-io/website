---
title: Nextflow and nf-core mentorship, Round 1
date: 2022-09-18
type: post
description: With the first round of the CZI Nextflow / nf-core mentorship now complete, we look back at the success of the program so far.
image: img/mentorships-social-card.jpg
tags: nextflow,nf-core,czi,mentorship,training
author: Chris Hakkaart
icon: chris_h.png
---

## Introduction

<div class="pull-right" style=" max-width: 350px; border: 1px solid #ededed; margin-left: 2rem; padding: 2rem;">
    <img src="/img/mentorships-round1-wordcloud.png" alt="Word cloud">
    <p><em>Word cloud of scientific interest keywords, averaged across all applications.</em></p>
</div>

Our recent [The State of the Workflow 2022: Community Survey Results](https://seqera.io/blog/state-of-the-workflow-2022-results/) showed that Nextflow and nf-core have a strong global community with a high level of engagement in several countries. As the community continues to grow, we aim to prioritize inclusivity for everyone through active outreach to groups with low representation.

Thanks to funding from our Chan Zuckerberg Initiative Diversity and Inclusion grant we established an international Nextflow and nf-core mentoring program with the aim of empowering those from underrepresented groups. With the first round of the mentorship now complete, we look back at the success of the program so far.

From almost 200 applications, five pairs of mentors and mentees were selected for the first round of the program. Over the following four months they met weekly to work on Nextflow based projects. We attempted to pair mentors and mentees based on their time zones and scientific interests. Project tasks were left up to the individuals and so tailored to the mentee's scientific interests and schedules.

People worked on things ranging from setting up Nextflow and nf-core on their institutional clusters to developing and implementing Nextflow and nf-core pipelines for next-generation sequencing data. Impressively, after starting the program knowing very little about Nextflow and nf-core, mentees finished the program being able to confidently develop and implement scalable and reproducible scientific workflows.

![Map of mentor / mentee pairs](/img/mentorships-round1-map.png)<br>
_The mentorship program was worldwide._

## Ndeye Marième Top (mentee) & John Juma (mentor)

For the mentorship, Marième wanted to set up Nextflow and nf-core on the servers at the Institut Pasteur de Dakar in Senegal and learn how to develop / contribute to a pipeline. Her mentor was John Juma, from the ILRI/SANBI in Kenya.

Together, Marème overcame issues with containers and server privileges and developed her local config, learning about how to troubleshoot and where to find help along the way. By the end of the mentorship she was able to set up the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline for the genomic surveillance analysis of SARS-Cov2 sequencing data from Senegal as well as 17 other countries in West Africa, ready for submission to [GISAID](https://gisaid.org/). She also got up to speed with the [nf-core/mag](https://nf-co.re/mag) pipeline for metagenomic analysis.

<blockquote>
    <em>"Having someone experienced who can guide you in my learning process. My mentor really helped me understand and focus on the practical aspects since my main concern was having the pipelines correctly running in my institution."</em> - Marième Top (mentee)
</blockquote>
<blockquote>
    <em>"The program was awesome. I had a chance to impart nextflow principles to someone I have never met before. Fully virtual, the program instilled some sense of discipline in terms of setting and meeting objectives."</em> - John Juma (mentor)
</blockquote>

## Philip Ashton (mentee) & Robert Petit (mentor)

Philip wanted to move up the Nextflow learning curve and set up nf-core workflows at Kamuzu University of Health Sciences in Malawi. His mentor was Robert Petit from the Wyoming Public Health Laboratory in the USA. Robert has developed the [Bactopia](https://bactopia.github.io/) pipeline for the analysis of bacterial pipeline and it was Philip’s aim to get this running for his group in Malawi.

Robert helped Philip learn Nextflow, enabling him to independently deploy DSL2 pipelines and process genomes using Nextflow Tower. Philip is already using his new found skills to answer important public health questions in Malawi and is now passing his knowledge to other staff and students at his institute. Even though the mentorship program has finished, Philip and Rob will continue a collaboration and have plans to deploy pipelines that will benefit public health in the future.

<blockquote>
    <em>"I tried to learn nextflow independently some time ago, but abandoned it for the more familiar snakemake. Thanks to Robert’s mentorship I’m now over the learning curve and able to deploy nf-core pipelines and use cloud resources more efficiently via Nextflow Tower"</em> - Phil Ashton (mentee)
</blockquote>
<blockquote>
    <em>"I found being a mentor to be a rewarding experience and a great opportunity to introduce mentees into the Nextflow/nf-core community. Phil and I were able to accomplish a lot in the span of a few months, and now have many plans to collaborate in the future."</em> - Robert Petit (mentor)
</blockquote>

## Kalayanee Chairat (mentee) & Alison Meynert (mentor)

Kalayanee’s goal for the mentorship program was to set up and run Nextflow and nf-core pipelines at the local infrastructure at the King Mongkut’s University of Technology Thonburi in Thailand. Kalayanee was mentored by Alison Meynert, from the University of Edinburgh in the United Kingdom.

Working with Alison, Kalayanee learned about Nextflow and nf-core and the requirements for working with Slurm and Singularity. Together, they created a configuration profile that Kalayanee and others at her institute can use - they have plans to submit this to [nf-core/configs](https://github.com/nf-core/configs) as an institutional profile. Now she is familiar with these tools, Kalayanee is using [nf-core/sarek](https://nf-co.re/sarek) and [nf-core/rnaseq](https://nf-co.re/rnaseq) to analyze 100s of samples of her own next-generation sequencing data on her local HPC environment.

<blockquote>
    <em>"The mentorship program is a great start to learn to use and develop analysis pipelines built using Nextflow. I gained a lot of knowledge through this program. I am also very lucky to have Dr. Alison Meynert as my mentor. She is very knowledgeable, kind and willing to help in every step."</em> - Kalayanee Chairat (mentee)
</blockquote>
<blockquote>
    <em>"It was a great experience for me to work with my mentee towards her goal. The process solidified some of my own topical knowledge and I learned new things along the way as well."</em> - Alison Meynert (mentor)
</blockquote>

## Edward Lukyamuzi (mentee) & Emilio Garcia-Rios (mentor)

For the mentoring program Edward’s goal was to understand the fundamental components of a Nextflow script and write a Nextflow pipeline for analyzing mosquito genomes. Edward was mentored by Emilio Garcia-Rios, from the EMBL-EBI in the United Kingdom.

Edward learned the fundamental concepts of Nextflow, including channels, processes and operators. Edward works with sequencing data from the mosquito genome - with help from Emilio he wrote a Nextflow pipeline with an accompanying Dockerfile for the alignment of reads and genotyping of SNPs. Edward will continue to develop his pipeline and wants to become more involved with the Nextflow and nf-core community by attending the nf-core hackathons. Edward is also very keen to help others learn Nextflow and expressed an interest in being part of this program again as a mentor.

<blockquote>
    <em>"Learning Nextflow can be a steep curve. Having a partner to give you a little push might be what facilitates adoption of Nextflow into your daily routine."</em> - Edward Lukyamuzi (mentee)
</blockquote>
<blockquote>
    <em>"I would like more people to discover and learn the benefits using Nextflow has. Being a mentor in this program can help me collaborate with other colleagues and be a mentor in my institute as well."</em> - Emilio Garcia-Rios (mentor)
</blockquote>

## Suchitra Thapa (mentee) & Maxime Borry (mentor)

Suchitra started the program to learn about running Nextflow pipelines but quickly moved on to pipeline development and deployment on the cloud. Suchitra and Maxime encountered some technical challenges during the mentorship, including difficulties with internet connectivity and access to computational platforms for analysis. Despite this, with help from Maxime, Suchitra applied her newly acquired skills and made substantial progress converting the [metaphlankrona](https://github.com/suchitrathapa/metaphlankrona) pipeline for metagenomic analysis of microbial communities from Nextflow DSL1 to DSL2 syntax.

Suchitra will be sharing her work and progress on the pipeline as a poster at the [Nextflow Summit 2022](https://summit.nextflow.io/speakers/suchitra-thapa/).

<blockquote>
    <em>"This mentorship was one of the best organized online learning opportunities that I have attended so far. With time flexibility and no deadline burden, you can easily fit this mentorship into your busy schedule. I would suggest everyone interested to definitely go for it."</em> - Suchitra Thapa (mentee)
</blockquote>
<blockquote>
    <em>"This mentorship program was a very fruitful and positive experience, and the satisfaction to see someone learning and growing their bioinformatics skills is very rewarding."</em> - Maxime Borry (mentor)
</blockquote>

## Conclusion

Feedback from the first round of the mentorship program was overwhelmingly positive. Both mentors and mentees found the experience to be a rewarding opportunity and were grateful for taking part. Everyone who participated in the program said that they would encourage others to be a part of it in the future.

<blockquote>
    "This is an exciting program that can help us make use of curated pipelines to advance open science. I don't mind repeating the program!"</em> - John Juma (mentor)
</blockquote>

![Screenshot of final zoom meetup](/img/mentorships-round1-zoom.png)

As the Nextflow and nf-core communities continue to grow, the mentorship program will have long-term benefits beyond those that are immediately measurable. Mentees from the program are already acting as positive role models and contributing new perspectives to the wider community. Additionally, some mentees are interested in being mentors in the future and will undoubtedly support others as our communities continue to grow.

We were delighted with the high quality of this year’s mentors and mentees. Stay tuned for information about the next round of the Nextflow and nf-core mentorship program. Applications for round 2 will open on October 1, 2022. See [https://nf-co.re/mentorships](https://nf-co.re/mentorships) for details.

<p class="text-center" style="margin: 30px 0;">
    <a href="https://nf-co.re/mentorships" target="_blank" class="btn btn-color btn-xxl">Mentorship Round 2 - Details</a>
</p>
