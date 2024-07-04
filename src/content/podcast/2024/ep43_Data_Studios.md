---
title: Data Studios
episode: 43
description: Iterative Nextflow development in the cloud & tertiary analysis.
date: 2024-08-06
type: podcast
subtype: Technical discussion
youtubeid: z4oTPtQNsjw
image: /img/podcast_ep43.png
tags: Nextflow, Seqera
author: Developer advocates
icon: logo_podcast_channels.jpg
---

In this episode [Phil Ewels](https://x.com/tallphil) hosts [Rob Newman](https://www.linkedin.com/in/roblnewman/) and [Rob Syme](https://x.com/robsyme) to explore Data Explorer and Data Studios.
These can revolutionise bioinformatics development, reducing friction with cloud computing and speeding up both pipeline scripting and tertiary analysis.
We talk about:

- How to seamlessly develop and troubleshoot Nextflow pipelines directly in the cloud.
- How to use the Nextflow CLI with super simple AWS Batch job submission _inside Data Studios_, for rapid iterative development.
- The power of Data Explorer for managing and browsing cloud storage across AWS S3, Azure Blob Storage, and Google Cloud Storage.
- Spinning up Studios with Jupyter Notebooks, RStudio, and Visual Studio Code for interactive, reproducible, and efficient tertiary analysis.
- How session persistence and real-time collaboration can transform team workflows.
- A sneak peek into upcoming features like custom containers and Data Studios session lifetime management.

Key links:

- [Data Explorer documentation](https://docs.seqera.io/platform/latest/data/data-explorer)
- [Data Studios documentation](https://docs.seqera.io/platform/latest/data/data-studios)
- [Rob Syme's Data Studios demo](https://summit.nextflow.io/2024/boston/agenda/05-24--practical-data-studios/) at the Nextflow Summit Boston 2024
- [Seqera feedback](https://feedback.seqera.io/) for ideas and feature requests
- [Seqera community forum](https://community.seqera.io/) for help and discussion
- [Barcelona Summit 2024](https://summit.nextflow.io/2024/barcelona/) - registration now open!

# Podcast overview

This episode focuses on how Data Studios changes the landscape for anyone working with Nextflow - providing an interactive environment for data analysis right within Seqera Platform.

## What is Data Studios?

Data Studios streamlines the often cumbersome tertiary analysis phase, enabling users to spin up ephemeral or timeboxed containers with various environments such as Jupyter Notebooks, RStudio, and Visual Studio Code. This feature allows real-time analysis and troubleshooting without needing extensive data transfers.

This can help span the "final mile" of analysis after the end of Nextflow pipeline runs, which require manual intervention and data transfers. Data Studios bridges this gap by providing a seamless, reproducible, and interactive environment for completing the analysis process.

## Key Features

### Reproducibility and Standardization

One of the core principles of Nextflow and Seqera is reproducibility. Data Studios standardizes environments to ensure consistent results, eliminating the unpredictable "wild west" scenario often encountered in final analysis stages.

### Moving Compute to the Data

A significant advantage of Data Studios is its capability to operate where the data resides; whether in AWS S3, Azure Blob Storage, or Google Cloud Storage. This reduces data transfer issues, making the process more efficient and cost-effective.

Data Explorer integrates seamlessly with Data Studios, allowing users to browse and interact with remote data repositories, select input data for pipelines, and preview or download data as needed. This ease of exploring and mounting data from various cloud providers directly into Data Studios significantly enhances workflows.

### Real-time Collaboration

A standout feature of Data Studios is its real-time collaboration capability. Distributed teams can easily work together within the same Data Studio session via a simple shareable link, enhancing productivity and communication.

### Session Persistence

Session persistence creates a complete snapshot of each Data Studios session when stopped. This ensures that no data is lost, providing continuity and reproducibility by allowing users to resume work from any checkpoint.

## Practical demos

### RStudio

The classic use case for Data Studios is to launch a notebook environment for interactive analysis. Rob Syme shows us how to load up RStudio in the browser, running in Studios on AWS. He mounts the results output from a [nf-core/differentialabundance](https://nf-co.re/differentialabundance/) pipeline run and launches the Shiny app that is bundled with the pipeline. In a new browser tab he then uses this customised analysis app to visualise the data analysis results interactively.

### Nextflow + VSCode ✨

Tertiary analysis is great, but what about developing Nextflow pipelines? Most bioinformaticians working with Nextflow typically work locally and push to the cloud once the code reaches a certain state of readiness. But sometimes it can be useful or necessary to run jobs on the cloud - especially if you need access to specific data or compute resources.

In this demo, Rob shows how to launch a Data Studios session running VSCode in the browser. Within that IDE he writes a minimal Nextflow pipeline and config and launches it via the VSCode terminal. The real magic is in how he takes the AWS Batch queues names from the Compute Environment used by the Data Studios session. Because the Studios session is running in the same environment all permissions are already set up automatically, making submission and infrastructure setup incredibly fast and simple.

He shows what's needed to set up the proper configuration files and how to use Data Explorer mounted data links in blob storage, with Fusion and Wave for seamless integration. The process emphasizes the benefit of handling complex tasks interactively without worrying about infrastructure, credentials, or data egress, thus simplifying the development and debugging of Nextflow pipelines.

### Custom Environments

Whilst he didn't actually demo it this time, Rob explains the demo he did for his [Nextflow Summit Boston 2024 talk](https://summit.nextflow.io/2024/boston/agenda/05-24--practical-data-studios/), using a custom environment running a graphical Ubuntu desktop to do analysis of microscopy data using Python GUI tools and applications like IGV.

## Looking forward

### Custom Containers

The big new feature that is coming soon is the ability to use Custom Containers, supplied by the end user. This is what was demonistrated at Boston (see above) and should soon be released for general usage.

This functionality will have a huge impact - allowing anyone to package any software for interactive usage on the cloud. We talk about how it will work on the back end and the kinds of things that we think people might use this new power-feature for.

### Managing Session Lifespan

Another important upcoming feature is the lifetime management of sessions. This allows for setting session limits to prevent runaway compute costs, ensuring that sessions are only active when needed and automatically stopping idle sessions. This feature is crucial for cost management and resource optimization.

## Find out more

The Nextflow Barcelona Summit will see Data Studios move to general availability, showcasing more advanced features and integrations. This event promises to be an unmissable opportunity for bioinformatics enthusiasts! Join us at [summit.nextflow.io](https://summit.nextflow.io/2024/barcelona/) to stay ahead of the curve in data analysis.

We'd love to hear what you think of these new feateures! User feedback is critical for us to refine and enhance Data Studios. Please do have a play, try to break it and tell us all your wildest ideas at [feedback.seqera.io](https://feedback.seqera.io/).

# Full transcript

:::note{title=Transcript}

## Introductions

**Phil Ewels:** Hello and welcome to Channels: the Nextflow podcast. You're joining us for Episode 43, going out on August the 6th 2024. Today we've got a really good lineup for you: We've got two Robs from Seqera joining me today. Two of many Robs actually that we have at Seqera. Rob Newman and Rob Syme. We're gonna talk today all about Data Studios, this really exciting new feature coming out of Seqera.

It's really relevant for anyone working with Nextflow. We're going to tell you why and how it works.

Rob Newman, can you tell us who you are and what you do?

**Rob Newman:** I'm a product manager for the platform. Specifically with the focus on anything related to data. So that's Data Explorer, Data Studios. I work with the engineering team to build out the functionality of Data Studios. I'm super excited to talk about it today.

**Phil Ewels:** Rob Syme

**Rob Syme:** Hi all, I'm Rob Syme, the scientific support lead here at Seqera. I'm part of the scientific development team. So we're the bioinformatician / computational biology dogsbodies. We end up doing a bunch of stuff across a bunch of different platforms. Previously to this I worked in academia and research environments.

**Rob Newman:** I would say, Rob, you keep us honest. When we're building products, you make sure that we're actually delivering what people need.

**Rob Syme:** I'd be using Seqera if I didn't work for Seqera.

**Phil Ewels:** It's a chunk of the company, have come from being users to now being developers or everything else. I'm definitely in that same category.

Anyone joining the podcast who hasn't heard my voice before, my name is Phil Ewels. I'm also a product manager at Seqera, working with open source.

## Data Studios

**Phil Ewels:** Okay, Rob, you mentioned these words, Data Studios, Data Explorer, what is this stuff?

**Rob Newman:** What we were seeing with our users was that they were running pipelines within our platform, and they had a lot of pain when it actually came to doing their tertiary analysis or troubleshooting workflows. We realized that it would be really great to be able to provide an environment within the platform that allows them to interactively troubleshoot their pipelines or review the data that was output from the pipelines. And also do their own tertiary analysis, within the data provided by their pipeline runs.

Data Studios was created to simplify that process. It allows you to spin up ephemeral or timeboxed containers with several different environments that we provide out of the box, whether that's Jupyter Notebooks or RStudio Notebooks or Visual Studio Code. You can mount data within your cloud storage buckets that are attached to the platform through Data Explorer, explore that data directly within your notebook environments.

So if you have a data science team, or a machine learning team, or analysts, or even just bioinformaticians who want to work through the data that is created by their pipelines, you can do that interactively in the platform using your existing credentials, your existing compute environments.

And it just makes the whole process that much more seamless.

**Phil Ewels:** That's a use case I'm very familiar with, coming from nf-core side. We have these pipelines and we always say: we know that there's a limitation in how far we can push things. We know that the nf-core/rnaseq pipeline is not going to run, and take you all the way up to your publication grade figures at the end. There's always going to be some final stage of the analysis.

We can standardize up to a point, but at the end of the day, every experiment is different. Every study is different. There's always going to be that final deep dive on the specifics, generation of figures and validation of hypotheses.

Until now, that's always been a case of: you run your pipeline, and then you download your data. And you go off and do whatever else. It's been a bit of a black hole in the Nextflow ecosystem.

**Rob Newman:** One of the founding principles of Nextflow and Seqera as a whole is this idea of reproducibility. Making sure that everything is standardized. With Data Studios, one of the things that we've seen with data scientists is that: you can spin up your own bespoke environment, with different libraries that are installed, whether it's the version of Python, whether it's the version of whatever library that you use.

With Data Studios, that's standardized. And it's pinned. So you know that if you're doing analysis within that Data Studio session, you're all using the same libraries. You will all get the same results from the analyses that you do.

**Phil Ewels:** Yeah, there's no point in being super rigorous about reproducibility for the first three quarters of your analysis. And at the end, it's just a wild west.

## Move the compute to the data

**Rob Syme:** I think one of the other nice things about it is: that if you're running on small scale data, or if you're doing preliminary experiments, you might think _&quot;Oh, I could just move my data around&quot;_. But as your team grows, and as the data volumes grow, particularly if you're using blob storage, suddenly moving that data around becomes a real problem. It might take time. It's a bit flaky.

The other nice thing that I really love about Data Studios is that you're moving the compute to the data. So you're operating on that _in-situ,_ where the data lives at rest: in S3 or in blob storage.

**Rob Newman:** One of the key pieces of information that we've heard from our customers is that when they want to do tertiary analysis, they're having to go into their cloud provider and spin up an EC2 container. And then they have to manage all the permissions to pull in the data. And then they have to manage that instance, or maybe it's not even them that are managing it. Maybe they have a system administrator who has to manage it. And that system administrator has to keep track of all the costs associated with it. And it's a huge headache for them.

So having this ability to spin up this ephemeral analysis container in the Platform is a huge time saver and efficiency saver. And is much more transparent.

**Phil Ewels:** We've had this catchphrase about moving compute to the data for years. And again, it's no point doing that for the three quarters of the analysis -if you then still have to download all the data at the end, to do the final bit.

It's nice to be able to actually push this the whole way and really follow those principles through to the end of your analysis.

## Real-time collaboration

**Rob Newman:** Another nice feature of Data Studios sessions is that you can do real time collaboration. Whether it's using Jupyter or VSCode, there's a simple one click where you can just say, _&quot;Hey, share this link with a collaborator&quot;_.

Provided that user has a role that has permissions to use Data Studios in the platform, they can click that link and they can go and do real-time collaboration in the same notebook environment.

**Phil Ewels:** Especially as more and more teams are distributed: the three of us are all sat in different countries, different time zones, the idea that you can just hop in and work together with someone as if they were sat next to you at your desk, is really powerful.

## Data Explorer

**Phil Ewels:** We've mentioned about mounting data there, and making data accessible. Can you tell us a little bit about what Data Explorer is and how this fits into this picture?

**Rob Newman:** Data Explorer allows you to browse and interact with remote data repositories from your organization workspaces in the platform.

It supports AWS S3, Azure Blob Storage and Google Cloud Storage. Based on the credentials that you attach to the workspace, it will automatically create these data links within the Data Explorer tab in your workspace.

This allows you to interrogate the contents of the buckets. So you can click down, making sure that all the data structures are there as you need.

If you have the correct permissions, you can preview and download the data. The previews that are supported are the Nextflow output files: the .command files, the fusion files. You can also preview text files like CSVs and TSVs, same with PDF files and also HTML files, and then images as well. If you have the correct permissions, you can actually download that data.

This is not only integrated with Data Studios, it's actually incorporated into the pipeline launch and the run details page. If you want to select input data for your Nextflow pipeline, you can open up the Data Explorer modal window, and you can drill down into your bucket and choose the data that you want as input data.

Once you have fired off your run and it's running, it's producing data, you can then go to the run details page, and you can click on the Data Explorer link and it will open up the output directory. In the Platform, in Data Explorer. You can start interrogating and looking at the files that are being created by your pipeline.

**Rob Syme:** I love that feature. The task details, and being able to drill down to a specific erroring task is one of my favorite features of the Platform.

Often in computational biology, you have these huge workflows that are spinning up hundreds, thousands, even tens of thousands of individual tasks. And if one of those tasks goes wrong, one 10,000th of a pipeline, being able to drill down to that specific error and say: &quot;what were the input files for this erroring out task?&quot;, &quot;Can I check the logs just for this one little piece of the pipeline in isolation, separated from the rest of the workflow?&quot;

Because, like a lot of computational work, things will go wrong eventually. It's only a matter of time. The work is complicated.

So being able to get that sort of precision debugging experience is really nice.

**Phil Ewels:** This is something I personally really missed when I switched from HPC to cloud. Before Seqera Platform existed, I was running Nextflow on an HPC and it wasn't too bad, because the task hash is there, and you could cd into the work directory pretty quickly.

But then when I started using cloud instead, it becomes much more onerous to actually get into that directory. i'd often end up having a browser tab there with a file explorer in AWS.

## Access to public data

**Rob Newman:** One of the really nice features of Data Explorer is that if you want to get access to public S3 buckets, like the 1000 Genomes project: that has a public S3 bucket. You can mount those public buckets in Data Explorer.

You can mount those as read only. You can explore them and use them an input to a pipeline, or whether it's helping with the analyses that you're doing in Data Studios, post pipeline run.

**Rob Syme:** A lot of public data sets, they say, &quot;Oh, we made the data public&quot;. &quot;It's on an S3 bucket&quot;. And some of those have an index listing, an HTML page that's generated, which shows you, but not all of them. And sometimes those HTML listings or indices fall out of date.

So it's really nice to be able to see exactly: What is the data that's in that public directory? What is the data that has been released?

Nextflow, will do the conversion for you. You can just apply an S3 path and use it as input to a pipeline. But sometimes, you need to know what that path is going to be. You need to actually get in there.

And yes, if you can drop down to the command line, run 'aws s3 ls' and step through that process, you can do that. It's not super fun. It's lovely to be able to just click through the data.

## Data Explorer demo

**Phil Ewels:** We talked about a few things here and visually, it's maybe more obvious as to how this works.

Don't worry if you're listening with audio only, you're not going to miss anything.

I'm just looking at the community showcase here on Seqera Platform. We've got Data Explorer along the top, So Rob you mentioned 1000 Genomes and I can see that's listed at the top here.

So if I click that, there you go. All these files and then click through different directories.

If I look at some pipeline run here, we can dive into specific tasks. So if I drop down and I find FastQC here. And then this is the new tab you were talking about, the Data Explorer tab. And I'm straight into the work directory.

That's so cool. And I can see all the HTML files. And they load. So it's literally two, three clicks from a pipeline run page.

**Rob Newman:** And also you can see the little link icon: you can copy that to the clipboard, and then you can paste that to a colleague, and just say, &quot;Hey, check out these results.&quot; So you don't have to be sending around files themselves. You can just say, &quot;Hey, just grab this link and go see it.&quot;

**Phil Ewels:** Also here, I was clicking around randomly in Data Explorer, it's the iGenomes here. This is on Azure. We've got Azure and AWS here, and Google Cloud as well, and all these different cloud providers sat together in one place. And I've got the same interface and I don't have to think really about how I'm accessing them.

**Rob Newman:** You can see in the top right of the bucket, it says download current directory. So you can choose to download everything in that current view. Or if you click the code button, this provides you with the cloud provider syntax to be able to download that whole view that you have within Data Explorer.

So this is hugely beneficial. Maybe you do want to download a bunch of files, and have those on your local host to do some analysis. Provided you have your credentials set up, whether that's AWS, Azure or GCP, you can very easily just copy and paste this into your terminal prompt, and you can download all the files locally.

**Phil Ewels:** Yeah, so I don't want to do this through the browser because it's about five terabytes. So it would take a while. it's nice. I don't have to go to Stack Overflow and try and remind myself the syntax of how to fetch things off the command line.

**Rob Syme:** It's nice that we give you both the in-browser experience, but also the ability and the help drop down.

**Rob Newman:** Downloading within the browser itself is limited, so we do put a cap on that of a thousand files, but within the command line interface, it's unlimited. Because you're just doing it locally on your machine. It's not going through the Platform.

**Phil Ewels:** Yeah, it makes a lot of sense. It's basically not using Data Explorer at that point. It's just using the native underlying systems.

I think that's a nice thread that I feel like we try and do quite a lot. We build a graphical interface to make stuff easier, but we don't hide the complexity too much. You can still drill down to the fundamentals pretty quickly if you need to. I think that's quite important.

So if I hop into the launch pipeline for Sarek, and then sure enough for the input and then the output directory, I've got browse buttons. Super cool.

And so now I can jump through all these different buckets I've got set up, select a bucket there in one click and I don't have to figure out the S3 path.

Happy days!

Sorry, a bit of a tangent there. I wasn't really intending to talk about Data Explorer quite as much, but it's quite a nice feature.

It's it's one of these things that I feel like I've endured quite a lot of frustration and pain myself personally. So seeing this stuff all slot together is really gratifying.

## Data Studios setup

**Phil Ewels:** Back to Data Studios. We've got our data links in the Seqera platform. We talked about how we can have all the different cloud providers side by side.

I guess that's all managed through centralized credential management within platform. So I don't have to think too much about that. It could even be set up by someone else. I don't have to know what the credentials are.

So now, what do I need to get Data Studios running?

**Rob Newman:** Typically, all you need is a compute environment that has the right credentials associated with it.

You could use an existing compute environment that had the right permissions associated with it. But you can go and create a new compute environment. One of the key pieces is that you need to define the credential associated with the compute environment, because that will dictate what you can and can't mount in terms of data.

If you want to read and write, within the compute environment creation process, in the &quot;allowed buckets&quot; form field, you have to list the of buckets that you want to write back to.

That's by design. We don't want people to be able to write data from their analyses to any bucket. We want to be very stringent about that. You don't want to pollute the data in your buckets. So you have to explicitly define which S3 buckets you're going to write your data back to when you create the compute environment.

Then when you create your Data Studio session, you can choose to mount data from one or more buckets. And those will show up as read and writable. They will also show up if they're in the same region as your compute environment, because the data ingress and egress costs aren't necessarily minimal, especially if you're pulling in really big files. So it's transparent to our users about costs that might be incurred.

**Phil Ewels:** That's another one of those things, which is really nice about working in the cloud is if you find a new public dataset and it's in the US instead of being in Sweden, I can just create a new compute environment in the US and then launch a Data Studio there. And suddenly all my data egress fees just go away.

**Rob Syme:** Creating a compute environment is made a little bit easier by the Platform with this forge feature. It's a wizard for creating compute environments. It'll do things like setting up the IAM permissions for you.

It'll set up the auto scaling groups, the AWS batch queues and compute environments, all those little pieces. Which are a bit of a pain, particularly if you're a biologist who doesn't spend a lot of time on cloud, it's lovely to be able to have the Platform stand that up on your behalf, inside of your account, ready to spin up either AWS Data Studio jobs or Nextflow jobs.

**Phil Ewels:** If you don't know what all those acronyms mean, like IAM and stuff, that's the point: you don't need to. Fantastic.

**Rob Syme:** You're living a better life.

**Rob Newman:** Also one of the nice things about the compute environment is you can specify the instance types. So let's say you're going to build a machine learning model in your Jupyter notebook. You would create a compute environment with GPU enabled instances. So you have the horsepower that you need.

**Phil Ewels:** We've set up our compute environment. We've mounted our data. Those things, I do as a one off, right? I don't need to create a new compute environment every time.

**Rob Newman:** It depends on the number of Data Studios that you create and have running at the same time. If they're all using the same compute environment, you have to be aware of those restrictions because you do define the number of CPUs for the compute environment. And if you create 20 Data Studios all with 10 CPUs, then you will potentially consume everything that the compute environment can provide.

Under the hood, it's all working within the batch, and currently only AWS batch is supported. So yeah, you've got to be cognizant of that.

But let's say you've got your compute environment created. So you can create your Data Studios. sessions by selecting the Data Studios tab in your on your workspace, and then there's a really nice, simple interface where you can just create or add a Data Studio.

Within that, you're defining the container template image that you want to use. Right now out of the box, we offer Jupyter, RStudio and VSCode. Jupyter and RStudio that's a notebook environment. When you are doing some tertiary analysis. You can also create a VSCode container. And that is super helpful for troubleshooting pipelines, or creating or updating pipelines.

So you select your container image template. You give it a name. You can customize that name. We do auto generate those for you and then you can give it an optional description as well. So maybe you're doing some very specific analysis and you want to add that metadata.

Exactly. Super imaginative. &quot;Hello, world.&quot; It's great.

And then you choose the compute environment. And again, the computer environment is super important because it provides the amount of resources that you can use for this Data Studio session.

Then you mount the data. So you click the mount data button, and that takes you to this Data Explorer modal, and it shows data at the bucket level.

**Phil Ewels:** What I like about this is it's showing a big red exclamation mark to warn me that this is not going to be free, because these are in different regions.

But even know this is Azure and my Data Studios is going to be running on AWS. That's fine, right? It's still going to mount the Azure data, so I can still access it within my AWS environment. Which still blows my mind every time I see this.

**Rob Syme:** We can mount the root of the bucket, but sometimes you want to be a little bit more specific. Maybe you have an enormous bucket that is currently hosting projects from lots of people or lots of different times. You can create data links, which points to a particular key or prefix you can think of as a directory inside of a bucket.

And so you're only mounting just a little piece of that bucket that you need for your tertiary analysis. So you can be as discriminatory or specific or as general as you want.

**Rob Newman:** And then you hit the add button. What it's doing under the hood is it's taking that container template image, which is just a docker file, and it's creating a new container within AWS batch.

It's creating a new job, essentially. Once you start it up, it's continuously running. And you can then go in within the Data Studios interface, hit the connect button and then it connects you to that running batch job under the hood, and you have your notebook environment or your VSCode IDE.

It's mounting data using our Fusion file system. So it's super fast, super responsive. And you can then go in and open up data.

So when your Data Studio session is running, it's incurring costs. It's a batch job that is in a running state. What you often want to do is start and stop that analysis container at will. Maybe you're doing some analysis for three or four hours and then you want to stop.

## Session persistance

**Rob Newman:** What we've recently rolled out is this idea of session persistence. Every time a running Data Studio session is stopped, a complete snapshot is taken and saved in the working directory of the associated compute environment.

It's created as a named checkpoint, and that checkpoint includes all the intermediary work and data not explicitly written back to the cloud storage bucket, and all the notebook configuration. So the environment set up the settings, the configuration and any dependencies.

And then you can, at any point, you can restart that Data Studio session. Or indeed, you can also start a new Data Studio session from that existing snapshot. And no data is lost. The state from the previous session is preserved.

**Phil Ewels:** It's almost like git history where you finish your analysis session and it creates a new commit and in the future you can branch off that, or you can continue working from that state.

**Rob Newman:** So this results in like much lower costs. You only turn on and off your Data Studio session or start and stop your Data Studio session on demand.

It's fault tolerant. It provides improved reproducibility and also portability of your analyses. Whether it's just for your individual user, or your collaborative teams working.

It also means that if you do accidentally leave your Data Studio session running, and you go off and get a coffee or grab lunch, or maybe you're done with work for the day. An administrator for the workspace can come in and they can safely stop that Data Studio session.

So costs don't needlessly escalate. It's just a one click. And if you compare that to what some of our customers were previously using, having to manage EC2 instances, that's a much more complicated process.

So administrators can cleanly stop a running Data Studio session, and the environment state won't be impacted. So it's a huge benefit.

Also, one of the really nice things about this session persistence is these checkpoints are created every time you stop a running Data Studio, but you can use those checkpoints as the basis for a new Data Studio session.

So let's say you want to do some exotic analysis. Maybe there's some library that you want to install, but you don't want to impact the parent Data Studio session.

You can go to the checkpoints tab within Data Studios, and you can say, Hey, I want to start a new Data Studio session from this pre-existing checkpoint, and it will create a new Data Studio session. It's essentially a fork of that environment, and you can choose to keep that environment and keep doing that exotic analysis, or you can choose to, roll that those changes back into the parent Data Studio session as needed.

## Data Studios RStudio demo

**Phil Ewels:** Rob, can you show us some of the more typical workflows that people might use this for?

**Rob Syme:** Yeah, sure. I think the thing that comes to mind for most people when they hear about these interactive analysis, is taking the outputs from a run and to do this work that we talked about at the beginning of the podcast: generating HTML reports or doing some sort of interactive analysis on these results.

So let's look.

I have a workspace here. It includes Data Explorer. And I have some data links here. So these links into a bucket, into a particular directory in that bucket. So this is the Data Explorer that we saw before.

So I could look at one of these directories, and I can see there are some files in here, including a shiny NGS app. These are different nf-core differential abundance outputs.

So I can see that data, but let's say I want to interact with the data. I can go to the Data Studios tab. I have a running Data Studio instance here.

Before I enter it, I'm just going to show you that I have mounted data in this Data Studio instance. So when I spun this up, I selected these two data links. So these two directories inside of blob storage to be mounted in my instance. And I've given them these names, differential abundance tests, pool results, and RNA seq test.

When I say I, but it was actually my colleague, John. He was actually one of the developers of the nf-core differential abundance pipeline.

And before we jump in I think it's also interesting to have a look at these checkpoints. These are exactly what Rob Newman was describing before. So I can see that this particular Data Studio has been stopped twice, and every time it's stopped, we have this checkpoint and I can give it a new name or if it's I can also start a new Data Studio based on any one of those checkpoints.

So maybe I've installed a bunch of software and I like the software versions that I have here. I'm going to save that and have that as a checkpoint. And maybe I start all my new analysis or all my new projects from that checkpoint.

Let's have a look. This is the Data Studio running. It's just a little RStudio notebook, so I can do all the normal RStudio things.

**Phil Ewels:** You hit connect there, and then we're still in the browser, right? This looks like RStudio, but we're still running in the web browser and the URL at the top, that's still seqera.io.

**Rob Syme:** Exactly. And in this URL, I could even save and pass to my colleagues. If they have access to that Seqera Platform workspace, all they need is that URL to connect and see these results and see this analysis.

So I can run some code here and this particular code takes the data that I've launched. So here on the root workspace is where you generally live inside of Data Studio under root workspace. And under workspace/data, you have these two data links that I mounted.

So if I mount more data links, they will appear as directories, which look like normal POSIX-y file systems to my RStudio code, my Python code.

It's a really nice thing not to have to think about blob storage. We've all had the luxury of dealing with Nextflow for a long time, and not had to really think about the differences between blob storage and local file systems, because Nextflow has this abstraction layer, which takes care of all that for us. I can just point it to an HTTP file and it'll work. I point it to some S3 blob storage and it'll work.

But in my custom analysis, often my R scripts will expect local files. I don't want to have to write code that deals with both blob or remote storage, like local files. But these emulate a POSIX-y file system. It's a FUSE file system under the hood. This is the fusion technology that Rob Newman was talking about. And I can browse these and treat them just like local file systems, even though they're files on remote storage.

**Phil Ewels:** And these could be public buckets we talked about before. You can mount anything here, no matter how big or who owns it, as long as there's access permissions.

**Rob Syme:** Exactly. So there's like that complexity that goes away. And the other nice thing is that I can read and write enormous pieces of data, to and from those blob storage.

A lot of the sort of work that is happening at Seqera is ensuring that computational biologists or bioinformaticians or researchers in general, don't have to think about infrastructure.

Writing directly to blob storage is a huge boon, because it means that I don't have to think about how much disk before I start my analysis, will I need to provision to accommodate all the things I want to do? Often if I'm doing interactive work, we don't necessarily know where it's going to lead. We don't necessarily know how complicated this job is going to get and how much disk we're going to need.

So it's very frustrating to start, spin up an EC2 instance begin your work and then realizing that you've only provisioned 200 gigs of data and you'll really need about 400 gigs. To have to tear that down, maybe attach some new block storage.

All of that complexity, it just gets in the way and it's nice to have that sort of fall away.

**Phil Ewels:** And inevitably you end up being super cautious and provisioning way more than you need, which then gets super expensive.

**Rob Syme:** Exactly. So you can treat S3 as this basically infinite disk, which is lovely.

So here just run some R code and this is a nice little interactive interactive exploration of my differential abundance results. So I

**Phil Ewels:** so this is a shiny app that's running out of the RStudio session and it's opened up a new browser window here, but even though it's the second browser window, it's still running in the same Data Studios session. Running off that background, R server that's in the RStudio session.

**Rob Syme:** Yeah, it's complicated deep stack here. So we have RShiny, running in an AWS batch job, running ECS agent, running on a EC2 VM.

But I don't have to think about any of that. It's just RStudio with my data.

**Phil Ewels:** Magic.

**Rob Syme:** It's really nice not to have to think about that stack.

## Nextflow development in Data Studios

**Rob Syme:** This is a sort of like the classic analysis, but there's something else that I want to show and that is interactive development work.

Often you will do a lot of Nextflow development on your local machine, using a small test data set. I strongly recommend you use the nf-core standards of having a test profile, which has the little test data set, that you can run in CI/CD or while you're interactively developing.

But there are often times where you are doing some debugging and you need a little bit more horsepower. You're dealing with data that is on blob storage.

I have a second Data Studios instance running here that I've called dev environment. And it is running a VSCode instance. So an image running VSCode.

I can connect to that. And here I'm dropped down to what looks like just VSCode, but it's running in the browser. I'm running in a Seqera Platform.

And I've got a little Nextflow project. This is some very advanced computational biology work going on here. You'll see, this is a workflow with a single process that creates a text file that says hello.

**Phil Ewels:** You should use that Data Studio I set up, it had a description, hello world.

**Rob Syme:** Yeah, exactly. Would be very appropriate.

But one of the nice things about this is that we are running this Data Studio in an AWS batch job, that has permissions, all the permissions that those AWS batch jobs run. When you submit Nextflow in Seqera Platform, Nextflow runs as an AWS batch job, and then submits all the component tasks.

So I know already ahead of time that this EC2 instance has permissions to submit jobs to AWS batch.

I can see this is running in this compute environment with this name &quot;on_demand_data_studios&quot;. If I look in the compute environments tab, I can inspect the details of this compute environment and it actually lists the head and the compute queue names.

So these are the names that correspond to AWS batch queues.

**Phil Ewels:** So Just to reiterate this is the compute environment that was spun up for the Data Studios. And when you do that with Tower Forge, as it's called, it, it creates all the AWS infrastructure for you. And AWS Batch works a little bit like an HPC job scheduler, I guess it's got these different queues.

And here, this is the native AWS Batch identifier for that compute queue.

**Rob Syme:** Right. And so this will work if you create the compute environment with Tower Forge. But even if you created it manually, if you're in a larger organization that has very specific policies about how these AWS environments will get set up, and you've created head and compute queues manually. Those compute queues will still be listed here in the details page.

And I can still submit jobs to them because I know that like these compute environments will need to have the permissions necessary to submit jobs.

So I can just copy that and I've added just a small amount of configuration. I've just dropped this in a nextflow.config in the root directory of the workflow. So I know it's going to get picked up. But I could just as easily pop it in my home directory. So in here, ~/.nextflow/config. And maybe I could put it under a profile, maybe AWS or remote or AWS batch.

Then I could switch between running tasks locally and , submitting tasks to AWS batch. This might be particularly handy if you have specific processes that are a bit larger and you would like those to be submitted to the remote AWS batch and have the rest run directly on your Data Studios instance.

**Phil Ewels:** Okay, so I could run a workflow in a Data Studio and most of it could run on the local executor, which would be the head job running on AWS batch. But if I have anything particularly big, a single process in the pipeline, which is heavy, I could tell it to use a different queue for that one process. And it would submit those tasks off to a different AWS batch queue, which has more heavy duty instance types. Is that right?

**Rob Syme:** Spot on. The caveat that you'll need to run Fusion and Wave to provision Fusion, probably inside of singularity containers for those local tasks.

I suspect that most people will probably just want to submit everything to AWS Batch. It's a lot cleaner and a lot simpler, but you have the option of doing that division if you want to get really fancy.

And for this podcast, like listeners to this podcast are the sort of people I imagine who probably do want to do the slightly more complicated and fancy things. So it's nice to know that's available. But for most people, these sort of five or six lines of configuration should be enough.

**Phil Ewels:** So we've got a config file. We've got process executor, AWS batch, and then you're just setting queue and that's where you've pasted in what you took from the other page.

**Rob Syme:** Exactly. This is the queue name that I've taken from the compute environment details page. The other important detail is I'm setting the work directory to somewhere in blob storage. A particular directory in a bucket that I know this compute environment has access to and it's at the region.

I've turned on fusion and wave. So this enables me to use this S3 blob storage as a work directory.

**Phil Ewels:** And that's exactly the same as if you were running a Nextflow manually, right? If you were submitting to AWS manually, this is very much how it would look.

**Rob Syme:** Yeah, exactly. So you can even do this on a local machine. If you want, you can use this sort of environment.

**Phil Ewels:** And if you're running locally, you'd have to set up all the credentials and all that kind of stuff. Whereas here you just magically have those credentials available!

**Rob Syme:** Exactly. It means that one person can worry about all those credential pieces. You might have one person who sets up the compute environment and needs to think about that. And then everyone else, the credentials are just done. You never need to worry about it.

So let's just run it actually and see if it does work. I've talked a big game here.

So you can see here, I'm in this directory main.nf, nextflow.config and I've actually symlinked here, the work directory. So the work directory that I've specified here in configuration was actually one of the data links that I mounted when I spun up this Data Studio. So even though I'm using the S3 as a work directory, I can cd into the work directory.

And see it as if it was a local local directory. So if I want to do some debugging, inspect the command.run.

**Phil Ewels:** And I can even see it in your sidebar there, in the explorer. It's just literally as if it you were running locally here. You're like hacking the system. Sort of making the AWS Batch part invisible.

**Rob Syme:** So if I run &quot;nextflow run main.nf&quot;. It's going to take a second because it's going to submit a job to AWS batch.

But this is a nice, much more interactive way of doing this sort of development when you need a little bit more oomph for your Nextflow pipelines.

I would still recommend if you're doing real work, that you want to make available to your colleagues, and you want to have all the nice provenance tracking, that you should submit those jobs via the Seqera Platform launchpad. You can of course run &quot;-with-tower&quot; here. But I think this sort of work, the benefit is only for development.

So you can, I can see here that it's run on AWS Batch.

**Phil Ewels:** Fantastic. This is very cool. This is something I theoretically knew should be possible with Data Studios, but I don't think I've ever actually seen someone do it.

So to see it, actually run like this. I love this. Personally, I'm a bit I'm a bit old school. I've been using HPC a lot longer than I've been using cloud. I think the promise of cloud is awesome. But I find it gets in the way quite a lot if I'm actually trying to do anything myself.

So this is perfect for me, because I can pretend that I'm not using the cloud, and I'm developing locally, or developing on an HPC. And it all works the same way that I'm used to, but under the hood, I get that infinite storage and I get all that kind of infinite scalability and all the benefits that come with running on the cloud.

And then we've got git installed in this compute environment. So if I'm happy with how this pipeline is running now, I've done a test pipeline run, I can then kind of Git commit, Git push, and then launch the same workflow through the Seqera Platform UI, as I would normally.

**Rob Syme:** Yep. You could even use the &quot;tw&quot; command line to launch it, directly from this Data Studios. So you wouldn't even necessarily need to leave Data Studios to launch the run, and have the full provenance tracking and have it available to your colleagues in the runs page here.

**Phil Ewels:** Very cool. Full circle.

I feel like there's going to be an avalanche of people trying this after listening to this podcast , we're going to see if the number of support requests suddenly jacks up as everyone tries it.

## Future development

**Phil Ewels:** So, let's start to bring this to a close then. Rob, where are we in terms of what's enabled and what's available to people listening to this today?

**Rob Newman:** It's in public preview. And it's also available for enterprise customers as well. If you're an enterprise customer, it does require a little bit of configuration for your Seqera Platform instance, but it is available it's shipped in the latest release of 24.1. For cloud users, it should be enabled for all users at this point, in all workspaces.

The plan is for it to go to general availability at the Nextflow Barcelona Summit at the end of October of this year. There are some new features that we're going to be rolling out. It is under heavy development. It's one of our core offerings now.

## Custom containers

**Rob Newman:** One of the key things that we're going to be providing in the next few upcoming releases are going to be this idea of bringing your own custom environments. We talked about how you have Jupyter, RStudio, and VSCode. Docker files under the hood, creating those containers for you in AWS Batch.

But for those of you that saw the Boston Summit presentation by Rob Syme and Florian: if you have the ability to define your own custom containers, the sky's the limit. You can have your own Linux machine. A remote desktop even, and you can install any of your own graphical user interface applications.

Which Rob did do in that demo. So this idea of actually being able to bring your own, or define your own custom containers is part of our roadmap and is the next big thing that we're working towards in the next couple of months.

**Phil Ewels:** I think out of all the things I've seen on the roadmap, it's probably the custom containers I'm most excited about. And I think the sentiment will probably be shared by many people listening and many people working in this space.

Being able to launch a Jupyter notebook is super cool. But we're not unique in being able to do that.

The idea of being able to launch whatever crazy kind of custom image I want to, in that environment space, that flexibility and that power is pretty unrivaled.

**Rob Newman:** There's a little bit of extra configuration, because it is using Fusion under the hood and the system, we've abstracted away. It looks really simple: that there's an RStudio Notebook or a Jupyter Notebook running in the browser. But we have our own web server called Connect, that is part of that. If you're going to build custom containers, you do need to both have these connect Components and also the Fusion file system installed as well.

There are some pieces that we are going to be developing to make that work. But at what it boils down to is: we're going to be using Wave; this super powerful piece of software that Seqera has created for doing this container augmentation on the fly.

And we're going to be leveraging that extensively in providing these custom environments for Data Studios.

**Phil Ewels:** We've talked about Wave in earlier podcast episodes. If people are interested in just how that works. The concept is that I can take any public container image, or build any image I want, throw it at Data Studios and Wave will basically add in all the components it needs en-route, without me having to do that manually.

**Rob Newman:** The releaseo of the Seqera Containers project, that also happened at Boston. You could define your own custom container there. The goal is to make those Data Studios compatible if you want to as well. So you could take that container that was created and is stored in Seqera Containers and just drop that into Data Studios as well. It even gets a security scan, and we're going to be leveraging all of those functionalities and capabilities within this upcoming release of Data Studios when it goes to general availability.

**Phil Ewels:** I love this. It's like we're starting to zoom out and all these different little components we've built, that have been isolated until now, they're all starting to link up. And you're like, ah, I get it now.

**Rob Syme:** All weaving together. It's really satisfying to see.

## Boston Summit demo

**Phil Ewels:** Rob, your talk at Boston got mentioned there. It went down pretty fantastically. Not many talks get a standing ovation, and it went semi viral on social media afterwards. There's a YouTube link for anyone who wants to watch it. I'll put it in the show notes.

Can you tell us a little bit about what you did? The Boston demo is quite different from what you showed just now, right?

**Rob Syme:** Yeah. Fundamentally different sort of workflow. We were demoing this Custom Containers feature. My colleague Florian had set up a container that ran an Ubuntu install, but also had a bunch of imaging applications.

Florian's background is in image analysis, which is an emerging field in computational biology. But a lot of the tools are built to run on a laptop. Like a GUI over a Python script. The tools are certainly not developed to run in a web server, or in that sort of way. So we needed a way of getting these tools and marrying those up with our cloud storage.

This custom container allows us to run those custom Python GUIs on a Data Studios instance, to read in data that was produced by an Nextflow pipeline.

We had an image analysis pipeline that produced some data, some test data. We took that test data interactively, manipulated it and generated a model, a segmentation model. Once we've done that interactive work, we took that segmentation model, which was written out to blob storage, and then used it as an input. into the Nextflow pipeline again. And ran our full analysis set.

So it was a way of joining together a Nextflow analysis pipeline, married together with some interactive analysis, and then straight back into another Nextflow pipeline, all without leaving the Platform.

**Phil Ewels:** Really putting proof to that idea of taking the compute to the data. That's the kind of thing that'd be uploading and downloading all this data, every step of the way. But being able to just go boom, and stitch it all together. Pretty fantastic.

**Rob Syme:** This is something that I think trips everyone up when they're first beginning computational work. Every time you move data around, It's so easy to get lost. And you end up with four copies of the data. &quot;Which was the data that you used?&quot; It's very easy to tie yourself in knots.

To have a single source of truth, just have the data at rest and operate on that where it is. Is a real boon.

I think it's an underrated time saver, a Human cost, both in times of frustration and hours lost chasing data.

**Phil Ewels:** And also monetary cost as well. Data egress is one of these things. It's a hidden cost for cloud computing, but it's where a lot of the charges come from, and it can really trip people up.

**Rob Newman:** Especially with the size of files that we're working with in bioinformatics. They're not just small files. There's multiple, very large files.

**Rob Syme:** Huge files, it's expensive to move around, but little files are expensive to keep track of in terms of Human costs. It works at both ends.

**Phil Ewels:** In that demo, if I remember correctly, it almost looked like you were using a desktop environment within the browser.

**Rob Syme:** Yeah, it's just an Ubuntu desktop running window management. So I could spin up a terminal, a little file browser. It just looked like a normal desktop.

**Rob Newman:** I think you showed IGV as well, right? Pulled up an IGV and did some exploratory analysis.

**Rob Syme:** Exactly, there's no limit. In these custom images, you can install whatever tools you want. And even if you don't have the tool installed, you could just open up a terminal. Run &quot;apt get install&quot; and pull down the latest tool.

So we had an IGV instance running that was loading a huge BAM file. It was sitting on the 1000 Genomes bucket, as if it was a local file, interacting with it and viewing an analysis really quickly and instantaneously.

**Phil Ewels:** We can talk about creating custom Docker images and it sounds a bit intimidating, but the idea of just loading up an Ubuntu desktop and running a regular install as if it was your local laptop is very cool.

**Rob Syme:** So if you want to do your own Docker files, absolutely support that. But if you'd prefer just to have a vanilla container and then do your own custom installs interactively: also supported. So you get to choose what sort of mode you prefer to operate at.

I think that in Seqera, the engineering team and that sort of people on this call, we all come from research backgrounds. We all know that it's a hubris for a provider like us to predetermine what sort of tools and analysis types and how you want to use the products.

We know that research is incredibly variable, and extraordinarily dynamic. So we can give you lots of options and you can explore, and we're excited to see methods and what tools are used by the community.

## Lifetime management

**Phil Ewels:** You said that custom containers is a big thing. It's coming. Is there anything else on the roadmap? Any teasers you can drop for us?

**Rob Newman:** When we've been doing demos with customers, one of the pain points that they also talk about is this idea of having, not runaway compute, but not having a lifespan or lifetime management. component part of these analysis containers or these analysis environments.

So this idea of defining some sort of session limits. So one of the configuration flags would be: &quot;I want to define a finite lifespan for this Data Studio session&quot;.

**Phil Ewels:** You can say &quot;I don't trust myself.&quot; &quot;I don't believe I'm going to remember to shut this thing down&quot;. So just set a limit just in case.

**Rob Newman:** It's also quite complicated to know the difference between: loading in a really big file and maybe you're doing some sort of re- formatting of a data table in Pandas. Or doing some sort of filtering. Knowing what's actually computation, rather than knowing when the user is idle and they've got up and gone to make a cup of coffee, or left for the day on a Friday afternoon.

Defining a lifespan where it says: &quot;Once that computation is finished, is this session now idle?&quot; And if it approaches that lifespan, to provide a prompt to the user to say: &quot;This is now an idle session. Are you still around? Do you still want this running?&quot; And if there's no response to automatically shut that down.

This idea of session management, I think is huge. It's something that many of our customers have said is a real pain point for them right now, when they're managing these independent EC2 instances for doing tertiary analysis.

**Phil Ewels:** So it's not just a simple countdown, it's actually checking what you're doing. It's like when you've been watching Netflix too long and it comes up saying, &quot;Are you still there?&quot; It's probably time to go to bed.

**Rob Newman:** I think what we're really looking for is our users of the Seqera platform to come in and start bashing on this and really seeing what works really well for them, what is challenging for them still, and then to submit feature requests through our feature request platform.

At feedback.seqera.io, there's a category for Data Studios and Data Explorer as well. You can come in and create feature requests. If there are other customers who have already added a feature request for Data Studios, you can upvote that feature request, or add your own comments. Then we take all of that input and use it to help define the product roadmap for Data Studios and actually for any part of the Platform.

**Phil Ewels:** I feel like almost every person I speak to about this: they stay quiet for a little while, you can see the cogs turning. And then it's like &quot;I could use this for...&quot; and then the thing they say next, is just so variable. Everyone has their own use cases.

**Rob Newman:** We're addressing the top of the iceberg. We're maybe addressing a third of what the use cases are. And then there's two thirds that are under the water. We don't even know what people are going to do. Part of the fun and the adventure of building this stuff, is being able to satisfy customers who have these challenging use cases. They're really not satisfied right now, with the processes that they're using to manage tertiary analysis. It's painful for them.

**Phil Ewels:** It's really exciting stuff. I'm really psyched to see where this goes, and sit in the audience myself as a listener to these talks in Barcelona and see what the next live demo will show. It gets better and better every time.

## Wrap up

**Phil Ewels:** Thanks so much for joining me today you guys. It's been an absolute pleasure chatting to you, and actually seeing some of these workflows and some of these ideas that we've been talking about for such a long time. But actually now materializing. I think it's really exciting. I'm sure that folks listening to this are also going to see the application of these things.

So, thanks for taking the time out and thanks to everyone who's been listening to the podcast. We hope that you're as excited as we are!

**Rob Newman:** I would like to add as well, I'm always amazed at what Rob Syme does. It blows me away. We build these things with these use cases in mind. And then Rob and the rest of the SciDev team really push it and show really what matters to people. And, really show the power of the things that we built. It's just amazing to see.

**Rob Syme:** And we're excited to see the researchers do the same thing at another order of magnitude of complexity.

**Phil Ewels:** Are you guys both gonna be in Barcelona?

**Rob Syme:** I hope so.

**Rob Newman:** Yes, I believe so. Yep.

**Phil Ewels:** So I'll see you both there, and also people listening: if you come to Barcelona Summit. Tickets are still on sale right now. You go to summit.nextflow.io. There's an early bird discount running at the moment. Come and see this stuff, chat to us in person, tell us what your use cases are. We'll be keen to hear about it.

**Rob Syme:** Looking forward to those conversations.

**Rob Newman:** Yeah, absolutely.

**Phil Ewels:** All right. Thanks very much, everyone. And until next time.

**Rob Syme:** Thanks Phil. Thanks

**Rob Newman:** Thanks Phil.

:::
