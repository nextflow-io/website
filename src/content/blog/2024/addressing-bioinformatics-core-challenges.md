---
title: "Addressing Bioinformatics Core Challenges with Nextflow and nf-core"
date: 2024-09-11
type: post
description: From managing complex pipelines to optimizing resource utilization, Nextflow offers a range of benefits that can streamline workflows and improve productivity for bioinformatics core groups.
image: /img/blog-2024-09-11--share.png
tags: nextflow,ambassador_post
status: published
author: Francesco Lescai
icon: francesco_nxf.png
community_post: true
ambassador_post: true
---

I was honored to be invited to the ISMB 2024 congress to speak at the session organised by the COSI (Community of Special Interest) of Bioinformatics Cores. This session brought together bioinformatics professionals from around the world who manage bioinformatics facilities in different institutions to share experiences, discuss challenges, and explore solutions for managing and analyzing large-scale biological data. In this session, I had the opportunity to introduce Nextflow, and discuss how its adoption can help bioinformatics cores to address some of the most common challenges they face. From managing complex pipelines to optimizing resource utilization, Nextflow offers a range of benefits that can streamline workflows and improve productivity. <!-- end-archive-description --> In this blog, I'll summarize my talk and share insights on how Nextflow can help overcome some of those challenges, including meeting the needs of a wide range of users or customers, automate reporting, customising pipelines and training.

### Challenge 1: running multiple services

_Challenge description: “I have a wide range of stakeholders, and my pipelines need to address different needs in multiple scientific domains”_

One of the biggest challenges faced by bioinformatics cores is catering to a diverse range of users with varying applications. On one hand, one might need to run analyses for researchers focused on cancer or human genetics. On the other hand, one may also need to support scientists working with mass spectrometry or metagenomics. Fortunately, the nf-core community has made it relatively easy to tackle these diverse needs with their curated pipelines. These pipelines are ready to use, covering a broad spectrum of applications, from genomics and metagenomics to immunology and mass spectrometry. In one of my slides I showed a non-exhaustive list, which spans genomics, metagenomics, immunology, mass spec, and more: one can find best-practice pipelines for almost any bioinformatics application imaginable, including emerging areas like imaging and spatial-omics. By leveraging this framework, one can not only tap into the expertise of the pipeline developers but also engage with them to discuss specific needs and requirements. This collaborative approach can significantly ease the deployment of a workflow, allowing the user to focus on high-priority tasks while ensuring that the analyses are always up to date and aligned with current best practices.

### Challenge 2: customising applications

_Challenge description: “We often need to customise our applications and pipeline, to meet specific in-house needs of our users”_

While ready-to-use applications are a huge advantage, there are times when customisation is necessary. Perhaps the standard pipeline that works for most users doesn't quite meet the specific needs of a facilities user or customer. Fortunately, the nf-core community has got these cases covered. With over 1,300 modules at everyone’s disposal, one can easily compose their own pipeline using the nf-core components and tooling. Should that not be enough though, one can even create a pipeline from scratch using nf-core tools. For instance, one can run a simple command like “nf-core create” followed by the name of the pipeline, and voilà! The software package will create a complete skeleton for the pipeline, filled with pre-compiled code and placeholders to ease customisation. This process is incredibly quick, as I demonstrated in a video clip during the talk, where a pipeline skeleton was created in just a few moments.

Of course, customisation isn't limited to pipelines. It also applies to containers, which are a crucial enabler of portability. When it comes to containers, Nextflow users have two options: an easy way and a more advanced approach. The easy way involves using Seqera Containers, a platform that allows anyone to compose a container using tools from bioconda, pypi, and conda-forge. No need for logging in, just select the tools, and the URL of your container will be made available in no time. One can build containers for either Docker or Singularity, and for different platforms (amd64 or arm64).

If one is looking for more control, they can use Wave as a command line. This is a powerful tool that can act as an intermediary between the user and a container registry. Wave builds containers on the fly, allowing anyone to pass a wave build command as an evaluation inside a docker run command. It's incredibly fast, and builds containers from conda packages in a matter of seconds. Wave, which is also the engine behind Seqera Containers, can be extremely handy to allow other operations like container augmentation. This feature enables a user to add new layers to existing containers without having to rebuild them, thanks to Docker's layer-based architecture. One can simply create a folder where configuration files or executable scripts are located, pass the folder to Wave which will add the folder with a new layer, and get the URL of the augmented container on the fly.

### Challenge 3: Reporting

_Challenge description: “I need to deliver a clear report of the analysis results, in a format that is accessible and can be used for publication purposes by my users”_

Reporting is a crucial aspect of any bioinformatics pipeline, and as for customisation Nextflow offers different ways to approach it. suitable for different levels of expertise. The most straightforward solution involves running MultiQC, a tool that collects the output and logs of a wide range of software in a pipeline and generates a nicely formatted HTML report. This is a great option if one wants a quick and easy way to get a summary of their pipeline's results. MultiQC is a widely used tool that supports a huge list (and growing) of bioinformatics tools and file formats, making it a great choice for many use cases.

However, if the developer needs more control over the reporting process or wants to create a custom report that meets some specific needs, it is entirely possible to engineer the reports from scratch. This involves collecting the outputs from various processes in the pipeline and passing them as an input to a process that runs an R Markdown or Quarto script. R Markdown and Quarto are popular tools for creating dynamic documents that can be parameterised, allowing anyone to customize the content and the layout of a report dynamically.
By using this approach, one can create a report that is tailored to your specific needs, including the types of plots and visualizations they want to include, the formatting and layouting, branding, and anything specific one might want to highlight.

To follow this approach, the user can either create their own customised module, or re-use one of the available notebooks modules in the nf-core repository (quarto [here](https://github.com/nf-core/modules/tree/master/modules/nf-core/quartonotebook), or jupyter [here](https://github.com/nf-core/modules/tree/master/modules/nf-core/jupyternotebook)).

### Challenge 4: Monitoring

_Challenge description: “I need to be able to estimate and optimise runtimes as well as costs of my pipelines, fitting our cost model”_

Monitoring is a critical aspect of pipeline management, and Nextflow provides a robust set of tools to help you track and optimise a pipeline's performance. At its core, monitoring involves tracking the execution of the pipeline to ensure that it's running efficiently and effectively. But it's not just about knowing how long a pipeline takes to run or how much it costs - it's also about making sure each process in the pipeline is using the requested resources efficiently.
With Nextflow, the user can track the resources used by each process in your pipeline, including CPU, memory, and disk usage and compare them visually with the resources requested in the pipeline configuration and reserved by each job. This information allows the user to identify bottlenecks and areas for optimisation, so one can fine-tune their pipeline for a better resource consumption. For example, if the user notices that one process is using a disproportionate amount of memory, they can adjust the configuration to better match the actual usage.

But monitoring isn't just about optimising a pipeline's performance - it's also about reducing the environmental impact where possible. A recently developed Nextflow plugin allows to track the carbon footprint of a pipeline, including the energy consumption and greenhouse gas emissions associated with running that pipeline. This information allows one to make informed decisions about their environmental impact, and gaining better awareness or even adopting greener strategies to computing.

One of the key benefits of Nextflow’s monitoring system is its flexibility. The user can either use the built-in html reports for trace and pipeline execution, or could monitor a run live by connecting to Seqera Platform and visualising its progress on a graphical interface in real time. More expert or creative users could use the trace file produced by a Nextflow execution, to create their own metrics and visualisations.

### Challenge 5: User accessibility

_Challenge description: “I could balance workloads better, by giving users a certain level of autonomy in running some of my pipelines”_

User accessibility is a crucial aspect of pipeline development, as it enables users with varying levels of bioinformatics experience to run complex pipelines with ease. One of the advantages of Nextflow, is that a developer can create pipelines that are not only robust and efficient but also user-friendly. Allowing your users to run them with a certain level of autonomy might be a good strategy in a bioinformatics core to decentralise straightforward analyses and invest human resources on more complex projects. Empowering a facility’s users to run specific pipelines independently could be a solution to reduce certain workloads.

The nf-core template includes a parameters schema, which is captured by the nf-core website to create a graphical interface for parameters configuration of the pipelines hosted under the nf-core organisation on GitHub. This interface allows users to fill in the necessary fields for parameters needed to run a pipeline, and allows even users with minimal experience with bioinformatics or command-line interfaces to quickly set up a run. The user can then simply copy and paste the command generated by the webpage into a terminal, and the pipeline will launch as configured. This approach is ideal for users who are familiar with basic computer tasks, and have a very minimal familiarity with a terminal.

However, for users with even less bioinformatics experience, Nextflow and the nf-core template together offer an even more intuitive solution. The pipeline can be added to the launcher of the Seqera Platform, and one can provide users with a comprehensive and user-friendly interface that allows them to launch pipelines with ease. This platform offers a range of features, including access to datasets created from sample sheets, the ability to launch pipelines on a wide range of cloud environments as well as on HPC on-premise. A simple graphical interface simplifies the entire process.The Seqera Platform provides in this way a seamless and intuitive experience for users, allowing them to run pipelines without requiring extensive bioinformatics knowledge.

### Challenge 6: Training

_Challenge description: “Training my team and especially onboarding new team members is always challenging and requires documentation and good materials”_

The final challenge we often face in bioinformatics facilities is training. We all know that training is an ongoing issue, not just because of staff turnover and the need to onboard new recruits, but also because the field is constantly evolving. With new tools, techniques, and technologies emerging all the time, it can be difficult to keep up with the latest developments. However, training is crucial for ensuring that pipelines are robust, efficient, and accurate.

Fortunately, there are now many resources available to help with training. The Nextflow training website, for example, has been completely rebuilt recently and now offers a wealth of material suitable for everyone, from beginners to experts. Whether you're just starting out with Nextflow or are already an experienced user, you'll find plenty of resources to help you improve your skills. From introductory tutorials to advanced guides, the training website has everything you need to get the most out of this workflow manager.

Everyone can access the material at their own pace, but regular training events have been scheduled during the year. Additionally, there is now a network of Nextflow Ambassadors who often organise local training events across the world. Without making comparisons with other solutions, I can easily say that the steep learning curve to get going with Nextflow is just a myth nowadays. The quality of the training material, the examples available, the frequency of events in person or online you can attend to, and more importantly a welcoming community of users, make learning Nextflow quite easy.

In my laboratory, usually in a couple of months bachelor students are reasonably confident with the code and with running pipelines and debugging common issues.

### Conclusions

In conclusion, the presentation at ISMB has gathered quite some interest because I believe it has shown how Nextflow is a powerful and versatile tool that can help bioinformatics cores address those common challenges everyone has experienced. With its comprehensive tooling, extensive training materials, and active community of users, Nextflow offers a complete package that can help people streamline their workflows and improve their productivity.
Although I might be biased on this, I also believe that by adopting Nextflow one also becomes part of a community of researchers and developers who are passionate about bioinformatics and committed to sharing their knowledge and expertise. Beginners not only will have access to a wealth of resources and tutorials, but more importantly to a supportive network of peers who can offer advice and guidance, and which is really fun to be part of.
