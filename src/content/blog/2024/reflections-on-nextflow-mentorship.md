---
title: One-Year Reflections on Nextflow Mentorship
date: 2024-04-10
type: post
description: From a mentee in the Nextflow and nf-core mentorship program to a Senior Bioinformatician at ZS. Come and read this post as I share a bit of my journay in the Nextflow community.
image: /img/blog-2024-04-10--share.jpg
tags: nextflow,nf-core,workshop
status: published
author: Anabella Trigila
icon: anabella.jpg
---

<div class="alert alert-info">
This post has been written by our valued community members 🤝
</div>

During December 2022 to March 2023, I was part of the second cohort of the Nextflow and nf-core mentorship program, which spanned four months and attracted participants globally. I could not have anticipated the extent to which my participation in this program and the associated learning experiences would positively change my professional growth.
The program aims to foster collaboration, knowledge exchange, flexible learning, collaborative coding, and contributions to the nf-core community. It is funded by the Chan Zuckerberg Initiative and is guided by experienced mentors in the community.
In the upcoming paragraphs, I'll be sharing more details about the program—its structure, the valuable learning experiences it brought, and the exciting opportunities it opened up for me.

<!-- end-archive-description -->

# Meeting my mentor

One of the most interesting aspects of the mentorship is that the program emphasizes that mentor-mentee pairs share research interests and that the mentor has significant experience in the areas where the mentee wants to develop. I found this extremely valuable, as it makes the program very flexible while also considering individual goals and interests. My goal as a mentee was to transition from a Nextflow user to a Nextflow developer.

I was lucky enough to have Matthias De Smet as a mentor. He is a member of the Center for Medical Genetics in Ghent and has extensive experience working with open-source projects such as nf-core and Bioconda. His experience working in clinical genomics was a common ground for us to communicate, share experiences and build effective collaboration.

During my first days, he guided me to the most useful Nextflow resources available online, tailored to my goals. Then, I drafted a pipeline that I wanted to build and attempted to write my first lines of code in Nextflow. We communicated via Slack and Matthias reviewed and corrected my code via GitHub. He introduced me to the supportive nf-core community, to ask for help when needed, and to acknowledge every success along the way.

<div style="margin-top: 2rem; margin-bottom: 2rem;">
    <img src="/img/blog-2024-04-10-img1a.png" alt="Mentor compliment about new module added" />
</div>

# Highlights of the program

We decided to start small, setting step-by-step goals. Matthias suggested that a doable goal would be to create my first Nextflow module in the context of a broader pipeline I wanted to develop. A module is a building block that encapsulates a specific functionality or task within a workflow. We realized that the tool I wanted to modularize was not available as part of nf-core. The nf-core GitHub has a community-driven collection of Nextflow modules, subworkflows and pipelines for bioinformatics, providing standardized and well-documented modules. The goal, therefore, was to create a module for this missing tool and then submit it as a contribution to nf-core.

For those unfamiliar, contributing to nf-core requires another member of the community, usually a maintainer, to review your code. As a newcomer, I was obviously curious about how the process would be. In academia, where anonymity often prevails, feedback can occasionally be a bit stringent. Conversely, during my submission to the nf-core project, I was pleasantly surprised that reviewers look for collective improvement, providing quick, constructive and amicable reviews, leading to a positive environment.

<div style="margin-top: 2rem; margin-bottom: 2rem;">
    <img src="/img/blog-2024-04-10-img1b.png" alt="Review comment in GitHub" />
</div>

For my final project in the mentorship program, I successfully ported a complete pipeline from Bash to Nextflow. This was a learning experience that allowed me to explore a diverse range of skills, such as modularizing content, understanding how crucial the meta map is, and creating Docker container images for software. This process not only enhanced my proficiency in Nextflow but also allowed me to interact with and contribute to related projects like Bioconda and BioContainers.

# Life after the mentorship

With the skills I acquired during the mentorship as a mentee, I proposed and successfully implemented a custom solution in Nextflow for a precision medicine start-up I worked at the time that could sequentially do several diagnostics and consumer-genetics applications in the cloud, resulting in substantial cost savings and increasing flexibility for the company.
Beyond my immediate projects, I joined a group actively developing an open-source Nextflow pipeline for genetic imputation. This project allowed me to be in close contact with members of the nf-core community working on similar projects, adding new tools to this pipeline, giving and receiving feedback, and continuing to improve my overall Nextflow skills while also contributing to the broader bioinformatics community. You can learn more about this project with the fantastic talk by Louis Le Nézet at Nextflow Summit 2023 Nextflow Summit 2023 - Louis Le Nézet
Finally, I was honored to become a Nextflow ambassador. The program’s goal is to extend the awareness of Nextflow around the world while also building a supportive community. In particular, the South American community is underrepresented, so I serve as a point of contact for any institution or newcomer who wants to implement pipelines with Nextflow.
As part of this program, I was invited to speak at the second Chilean Congress of Bioinformatics, where I gave a talk about how Nextflow and nf-core can support scaling bioinformatics projects in the cloud. It was incredibly rewarding to introduce Nextflow to a community for the first time and witness the genuine enthusiasm it sparks among students and attendees for the potential in their research projects.

<div style="margin-top: 2rem; margin-bottom: 2rem;">
    <img src="/img/blog-2024-04-10-img1c.png" alt="Second Chilean Congress of Bioinformatics" />
</div>

# What’s next?

The comprehensive skill set acquired in my journey resonated positively and allowed me to join the ZS Discovery Team. This organization accelerates transformation in research and early development with direct contribution to impactful bioinformatics projects with a globally distributed, multidisciplinary talented team.

In addition, we organized the first Nextflow Hackathon in Argentina, fostering a space to advance our skills in workflow management collectively. It was a pleasure to see how beginners got their first PRs approved and how they interacted with the nf-core community for the first time.

<div style="margin-top: 2rem; margin-bottom: 2rem;">
    <img src="/img/blog-2024-04-10-img1d.png" alt="nf-core March 2024 Hackathon site in Argentina" />
</div>

My current (and probably future!) day-to-day work involves working and developing pipelines with Nextflow, while also mentoring younger bioinformaticians into this language. The commitment to open-source projects remains a cornerstone of my journey and I am thankful that it has provided me the opportunity to collaborate with individuals from diverse backgrounds all over the world.

Whether you're interested in the mentorship program, curious about the hackathon, or simply wish to connect, feel free to reach out at the nf-core Slack!

<br>
<div class="footer-wrapper">
<img src="/img/nextflow_ambassador_logo.svg" height="70" class="pull-right" style="padding: 1rem 1rem;">

> This post was contributed by a Nextflow Ambassador, passionate individuals who support the Nextflow community. Interested in becoming an ambassador? Read more about it <a href="https://www.nextflow.io/ambassadors.html">here</a>.

</div>
