---
title: New syntax in Nextflow 24.04.
episode: 41
description: Workflow output definitions, topics, eval and more.
date: 2024-07-09
type: podcast
subtype: Technical discussion
youtubeid: SFwcOMXKzVE
image: /img/podcast_ep41.png
tags: Nextflow, Seqera
author: Developer advocates
icon: logo_podcast_channels.jpg
---

Join us in exploring the latest Nextflow release, 24.04. [Phil Ewels](https://x.com/tallphil) and Ben Sherman dive deep into the new updates in Nextflow's most recent version. We discuss advanced retry strategies, job arrays, resourceLimits, Singularity’s OCI mode, and the game-changing _Workflow Output Definition_.

<!-- end-archive-description -->

:::tip
We're actively seeking feedback from the community on these new features.
**We would love for you to experiment with the new syntax and let us know your thoughts.**

Please let us know what you think either via [Nextflow GitHub issues](https://github.com/nextflow-io/nextflow/issues) or via the [community forum](https://community.seqera.io/).
:::

:::note
In this episode we used a new tool, [Descript](https://www.descript.com/), to transcribe and edit.
The result is many small edits but fewer _"uhms"_ and _"errs"_.
It also helps with the transcript and summary texts (see below).

Do you like it, or should we go back to how episodes were before?
Let us know either in [Slack](https://nextflow.slack.com/archives/C04E9USS7M2) on in the [community forum](https://community.seqera.io/c/community/13)!
:::

In this episode we refer to a the recent blog post about the release, which you can find here: [Nextflow 24.04 - Release highlights](https://nextflow.io/blog/2024/nextflow-2404-highlights.html).

# Podcast overview

## Intro - tweaks, fixes and new support in 24.04

### Performance and Stability Improvements

Ben kicks off the discussion with key performance and stability improvements. The latest stable release, 24.04, includes numerous bug fixes and performance enhancements. There's a significant focus on closing gaps and adding retry strategies for better stability, especially concerning cloud providers. Notably, if Nextflow’s API calls against cloud providers fail for any server-side reason, it will retry automatically, avoiding pipeline failures caused by temporary server issues.

### Publishing Changes

Previously, if a file failed to publish, Nextflow would only issue a warning. It was possible for a pipeline to complete successfully without noticing missing output files. Now, by default, Nextflow will fail the pipeline if the publishing fails, though there's an option to revert to the old behavior. Additionally, retry strategies for publishing have been introduced to ensure retry attempts if an issue arises.

### Job Arrays

One highly requested feature is the introduction of job arrays. This allows users to submit many jobs as a single submission, alleviating strain on schedulers. The submission happens in a batch, and then the scheduler can process and plan effectively. Once a job is submitted as part of a job array, the jobs run independently. If any child job fails, it is resubmitted without affecting the rest of the array.

### Singularity OCI mode and GA4GH TES

In the past, Singularity supported Docker images by converting them into SIF files, consuming storage and time. Now, both Singularity and Apptainer can run OCI images directly, saving valuable resources.

Additionally, the TES executor has seen significant improvements, now supporting TES 1.1, which brings broader compatibility and integration with existing workflows.

## Major new syntax features

### Topic Channels

A new feature called topic channels offers a more straightforward approach to collecting channel outputs across pipelines. Channels can emit data to a named topic, simplifying the collection and use of version information from various processes within Nextflow pipelines.

### Eval Outputs

Eval outputs simplify the addition of shell commands to tasks. With eval outputs, necessary post-task commands can be defined neatly, avoiding repetitive code within process definitions.

### Workflow Output Definition

The concept of workflow output definitions has been introduced. This new syntax streamlines the publication of files by defining publish targets within workflows. Instead of defining publication behaviors within process definitions, users can now manage them at the workflow level, ensuring better clarity and fewer repetitions in the code.

# Full transcript

## Introduction

**Phil:** Hello and welcome to channels for Nextflow podcast. You're joining us today on episode 41, going out on July the 9th, 2024. My name is Phil Ewels. I'm joining you from Stockholm in Sweden, and I'm a product manager at Seqera and today we've got Ben Sherman also from Seqera joining us to talk all about some of the new Nextflow syntax.

As many of you listening will know, we had a new stable release of Nextflow that went out quite recently. Nextflow's release model has monthly edge releases, and then every six months we do a stable release. So the most recent release went out just at the start of May actually, but it's called 24.04. We've had a couple of patch releases, but all the main stuff is in that stable.
In today's episode we're going to walk through some of the new things that we introduced in that version, some of which are in preview, and we really want to get your feedback about how that syntax works.

The big topics are going to be `topic` channels `eval` outputs. And then the really big one is the new _workflow output definition_. And we're going to dig into those, but to start with, whilst we've all got your attention, we're going to rattle through a couple of the smaller things.

Ben, do you want to lead us off? What have we got first?

## Retry strategies

**Ben:** Yeah. The first thing, I wanted to just go through really quickly is first of all, as usual, the stable release has just a ton of performance, stability, improvements, bug fixes, things of that sort. We're always, trying to close gaps like that.

Notably we've added a bunch of retry strategies to a lot of stuff under the hood. Whenever Nextflow makes API calls against different cloud providers making sure that, If that request fails for some server side reason that we try it. That's because that shouldn't be something that fails your pipeline.

Sometimes the cloud service will just return a 500 or 502 or something like that. And there's nothing you can do about it except retry. So Nextflow should do that for you. Also with publishing, we've changed some things a bit. Actually, two important points I would say is it used to be that if a file was failed to be published Nextflow would give you a warning, but it wouldn't actually fail.

And so what that meant was that it was possible for the pipeline to succeed and to say that everything was fine, but then actually you had some missing output files and I feel if the output file is missing, then you've missed the whole point. And so now we've changed it to by default it will fail the pipeline if the publishing fails and then there's an option in the publishDir that you can do to revert back to the old behavior if you want.

And then alongside that, we added a retry strategy for the publishing. So if the publishing fails for any reason we'll just go ahead and retry it. This also can happen, especially when uploading really large files. For instance like on S3, when you have to publish a very large file, you actually have to do a multi part copy.

S3 does that in the background. And then it has to merge all the chunks. And sometimes, there's a checksum failure or something. And so we just retry that.

**Phil:** Is that retry done by default or is that something that folks need to opt into with the config?

**Ben:** The retry is enabled by default. And then there's a new section of config settings that you can use to control the different parameters of that. So there's a lot of different things in Nextflow that have a retry strategy. So this is like publishing with the executor, when submitting jobs.

Again, with AWS and Google APIs there's config settings for all of those different things to control the retry parameters. So things like how many retries, the retry delay, stuff like that. You can go change all of those, but they're all set to sensible defaults.

So nextflow.publish is the config setting there. We're probably going to change that section to workflow.publish at some point in the future or workflow.output.

**Phil:** One last question, and this is for my benefit when it fails, does it crash the pipeline wherever it is, even if it's halfway through, or does it wait to the end or what mode of failure does that look like?

**Ben:** I believe it will actually fail immediately

**Phil:** okay. So we're more stable. What's next?

## Job arrays

**Ben:** Just a couple quick features enhancements.

First one is job arrays. This is actually a whole new feature, highly requested one too. There's just not much to say about it because there's not much you have to do. It's just a directive called array. And if you want to use job arrays, you just say array and then a number, so like 10 or 100 or whatever.

There's a lot we can say about how it works under the hood. So it's essentially a job array is a way to submit many jobs as a single submission. The reason this is useful is because schedulers can easily get overwhelmed.

So if you submit, thousands or tens of thousands of jobs at a scheduler all at once that can quickly overwhelm the scheduler because it has to process each job, even if it's not going to submit it yet, it has to add it to its queue and all this different stuff. And so it's a lot easier on the scheduler.

If you give it. A single job, which describes, say, a hundred different child jobs, if they all use the same submit script and parameters and things like that then it can process all of them much quicker. And on top of that, it can plan a bit better think if when you get all those jobs in batch, instead of getting, one at a time.

It's like someone texting you like just a little text message every second versus they just, maybe they just give you one big batch text message. And then you can process it how you want. The important thing to note is that. Even though you're submitting all the jobs in a single batch, you're not actually running them as a single job. So once the submission happens the jobs are essentially handled Like independent jobs, just as before. It's just that the submission happens in batch.

And so then after that, the scheduler will run each job independently, as it sees fit again, having the advantage of being able to plan, how it wants to do that. It can do that more efficiently. And in terms of how Nextflow works, if any of those child jobs fails, then they'll just be resubmitted.

So it won't try to like redo the entire array just because one of the jobs failed. And you won't see any difference in terms of the Nextflow, log or anything like that. It'll just look like the same tasks being submitted.

And also the main trade off that you make by doing this is that normally a job might be submitted immediately, as soon as it's created by Nextflow. Now with job arrays, Nextflow will wait until there is a whole array of jobs. Because then it submits the whole array.

if your pipeline is such that, a bunch of jobs typically get created all at once, then you don't really care about that delay.

But if your pipeline is structured in a way where the jobs are trickling in very slowly, such that if you were to use a job array, it would actually significantly delay the submission of those first few tasks, then job arrays might not make sense. Of course, in that situation the job arrays probably don't even help that much because, you might not even be overwhelming the scheduler

**Phil:** What kind of task volumes are we talking about? Presumably you have to be well into the thousands before this really is something anyone needs to worry about.

**Ben:** The actual limits probably depend on the scheduler and computer environment that you're using. I would say thousands is probably a good point at which to do it. And it's a process directive, so you can apply it to specific processes where it's appropriate if you want rather than just applying it across the board.

**Phil:** And this can be set in config files, presumably? So it's not something that needs to be hardwired into a pipeline?

**Ben:** Yeah, I think this is typically something you would, you could put it in your pipeline code, but I don't think you would typically do that. So yeah, you can put it in your config.

It also works on some cloud providers. So there's a bunch of HPC schedulers that support it. You can look at the docs for that. And Google Batch and Amazon Batch both support it as well.

With AWS Batch, a useful thing about this is that when you use job arrays if you cancel your pipeline or your pipeline is terminated for any reason, Nextflow can actually delete the entire job array with a single request. And this is significant because normally Nextflow can only delete those jobs one at a time, and it can only do that so fast.

Say Nexflow is running on a VM, like an EC2 instance, and then, you cancel that VM, it'll basically send a SIG term to the EC2 instance, and then it has, I think, 30 seconds to shut down, and then it'll get killed.

And so if you do the math Nexflow can do about 50 requests a second, against AWS Batch. That means that Nexflow can delete, at most, you know, 1500 jobs if you have to cancel that pipeline. And obviously if you have more than that, then they won't get deleted and they'll just sit around.

This is something that we want to improve, at a more system level, but for now one remedy to that is to use job arrays. Say you use job arrays of. 100 well, now that bumps the number of jobs Nextflow can delete to 150, 000 jobs. And that's probably going to cover you

and with Google batch a nice thing about job arrays Google batch is that.

So Google batch, there's this concept of, the task requirements. So the task can have a certain amount of CPUs and memory that it needs. And then there's also the concept of the VM or machine type, which has a certain number of CPUs in memory. What you can do with job arrays is specify a really large machine type.

Suppose it has eight CPUs and 16 gigs of RAM, and then you specify the tasks as needing, one CPU and two gigs of RAM. And then if you say array eight, then it will pack all of those eight jobs onto the same VM. Which is a nice Google batch gives you more control over how this packing works.

Normally, every task would be assigned to one VM, which would really suck then if you only needed the one CPU and two gigs.

**Phil:** The Google batch API should do some packing itself as well like it should stick more than one task on an instance if there's space?

**Ben:** Yes, but basically what we found was the only way to actually make it do that was through the job arrays feature, essentially. Otherwise, it doesn't really work.

**Phil:** That's good to know that more efficient usage on cloud is never a bad thing. What does Nextflow do if there's 46 tasks for a process?

**Ben:** So Nextflow has the ability to know basically when all the jobs have arrived for a process. It's like a process knows, when the process is &quot;closed&quot;, that means that all of the inputs have arrived. And so once that signal happens the job arrays will essentially be flushed.

So any remaining jobs will be submitted as a partial job array. if any job is retried for any reason, it's just submitted individually.

The reason for this is that, it's common to have things like a dynamic resource requirements. So you want to increase the CPUs or memory on the second attempt.

In that case, it doesn't make sense to submit it as part of a job array because one of the requirements is that all of the submit related directives: so things like CPUs, memory disk, that actually has to be uniform across the job array. Now a lot of process directives can still be different across the tasks. So things like if you're using like conda packages or even the container, those things can be completely different because Nextflow was just found a way to basically, keep that stuff independent. But anything that gets mapped to a submit option.

Like the CPUs and memory, those have to be uniform. And so we just as a simple solution, we just say, anytime a task is retried we'll just submit it separately.

**Phil:** So what happens if I set a memory that's dynamic, based on the input file size , and then try and set array on that as well?

Is Nextflow smart enough to spot that? And not do arrays or do arrays in a different

**Ben:** You know what I actually don't think it would catch that. I think what would happen in that case is basically whichever job gets added to the front of the array, that those settings would be used

**Phil:** there's probably fairly good behavior, to be honest, because I would guess in most pipeline runs, most samples would be fairly homogeneous.

**Ben:** Yeah, it is an edge case that has been in the back of my mind for a while. I would like to find a better way to, detect it. But that's the tricky part just because it's dynamic doesn't mean that it won't work if it's dynamic based on the task attempt then that's fine. But then if it's dynamic based on the inputs, then ideally I would like Nextflow to be able to detect that.

So I think for now we just say in the docs you need to be careful about this. You shouldn't use that level of dynamic dynamicity.

## resourceLimits

**Ben:** Some other features I wanted to mention was we added a new directive called resourceLimits. Anyone at nf-core will know exactly what I'm talking about with the whole check_max function. This is basically a native version of check_max. So you can specify things like max CPUs, memory, and time. And those will represent like a system wide limit. So if you know that you have a cluster and the largest node in that cluster has, 32 CPUs and a bajillion gigs of RAM and your, maybe your queue the max wall time is 72 hours.

Then you can specify all those things as resource limits, and then you can continue to have your dynamic resource strategy where, you increase the resources on each attempt. But if those resources happen to exceed one of those system wide limits, Nextflow will automatically just bump it down. to the system wide limit. And so that way you avoid this issue of the job getting rejected by the scheduler or in the worst case, the job just sitting there forever because the scheduler will never be able to schedule it.

**Phil:** I think this is, this was maybe one of my oldest issue requests. Think I put this one in six years ago or something. That little work around hack that I just did short term, like six years ago, for nf-core has kicked around ever since. I'm very happy to see this come in- that means we can finally get rid of that.

**Ben:** We're also trying to tighten up the config syntax and not just allow people to put whatever they want in there. And so this was an important part of that because you just had this whole groovy function sitting in the config file, and we don't want to really do that anymore.

**Phil:** Also everyone had to copy that into their own config files. Everyone assumed that by increasing max_memory, it would mean that tasks got more memory, which is not what happened. And it's just riddled with problems, so I'm very happy to see the back of it.

## OCI mode for Singularity

**Ben:** And then two more things. So one is OCI mode for Singularity. The main thing I want to get across is just that in the past Singularity, , supported Docker images, but the way it would do it is it would download the Docker image and then convert it to a SIF file. That's the Singularity Image Format, and then it would run it. You would need storage for that SIF file and need time to do the conversion.

Now, Singularity and Apptainer, both have much better native support for OCI images, which are essentially Docker compatible images. So now instead of converting the file, it can actually just run the OCI image directly. And so if you're using Singularity, this can, save a bit of time for you.

**Phil:** And it does that on the task node as well, instead of a head node.

**Ben:** Yeah, so now Singularity can do it when the task is running. And just keep all that local.

**Phil:** Which is good because I had a common complaint from sysadmins is that Nextflow is using lots of CPU and actually it's not Nextflow, it's Singularity. So that's good. I'd like to flag with this one is that this is only really useful if you _only_ have a Docker image. So if you've got native Singularity images, you should _not_ be switching this on.

**Ben:** Yeah, there's no need to use it in that case. You've already done the work.

**Phil:** And I think it forces the use of a Docker image if there's both. So it's better to just stick with the previous settings.

I also want to flag this really fantastic nf-core bytesize talk by Marco who's at Seqera. He really went into all the cutting edge, state of the art for Singularity, including all the stuff we're talking here with loads of examples and all the different settings and how they all work. So if you're interested in Singularity and OCI mode , go and check that out.

It was recorded on YouTube. It's not yet added to the nf-core website, but if you go to the nf-core channel. You will find it there. There you go: containers for HPC Singularity.

**Ben:** Marco knows way more about containers and Singularity than I do. He and Paolo figured all that stuff out with the OCI mode. there was a whole discussion trying to figure out the exact right way to do it.

## TES executor

**Ben:** The last quick thing I wanted to mention was the TES executor. So if you're familiar with GA4GH, they have all these API standards.

And one of them is called TES, which is the task execution service. It's essentially a standard that's meant to standardize interfaces for compute backends. So think like all the different executors that Nextflow supports, the idea is what if there was one API that could wrap around all those different kinds of backends and then as a workflow manager, you don't have to integrate with all these different executors. You just integrate with TES and then it just works.

Of course, I'm always reminded of that XKCD comic about the people making yet another standard, because in reality, we're not just going to drop all the existing executors. We're just going to add TES as another executor, which is exactly what we did.

But there was some implementation gaps, I guess you could say for a long time with Nextflow. There were certain things that Nextflow does, that TES didn't support. The main things were that bin directory. And then also things like glob outputs and directory outputs.

And TES graciously, they actually updated their API standard to accommodate us. We did some work to update our own integration in Nextflow. And so now Nextflow supports TES 1.1. So if you're using a TES 1.1, compliant executor, then you should now be able to run any Nextflow pipeline.

**Phil:** There's another podcast episode where you and Geraldine talk about all the GA4GH standards.

**Ben:** Really interesting space there. They're starting to gain some traction, I think. And there's increasing interest in supporting some of these other standards, things like DRS, for, obtaining data from different storage backends.

**Phil:** I think also the TES just had a publication as well. You can go and dig out and have a read.

Anything else before we get into the big boys?

**Ben:** Nope, I think it's time for the big boys.

So there's two big _topics_, no pun intended, that we wanted to dig into today, which is around new Nextflow syntax. Both of these things I would say are experimental syntax in the sense that they're not finalized.

They work, we've tested them, so they're good in that sense. But we'd really like to get more feedback from the community over how exactly to shape this syntax. And so we wanted to go a bit deeper into it today than maybe what the blog post does, or what Paolo showed at the Boston summit. I've got the blog post here.

**Phil:** So anyone listening, if you go to nextflow.io and go to the blog, Nextflow 24.4 release highlights. So that's what we're looking at here.

## Topic channels

**Ben:** Yeah - first thing I wanted to show here was this new thing called &quot;topic&quot; channels. Paulo and I have a minor philosophical disagreement over whether these should be &quot;topic channels&quot; or &quot;channel topics&quot;. Depends on how you think about it. I'll leave that up to your judgment, but essentially it's a new way to collect values across the pipeline.

One of the problems that we noticed was that the nf-core pipelines they were trying to collect the version of every tool that they use across the different processes.

And so every process has an additional output channel called versions, which simply includes the tool version. And then every time you're invoking a process, you're using the mix operator to collect these version channels into one big channel and then using it at the end.

And so Paulo came up with a two pronged way to simplify this. The one prong is the topic channels, and the other prong is the eval outputs, which I'll get to in a second.

The topic channel is like this global channel that you can use. And so when you have an output channel like this in a process, you can say &quot;send it to this topic&quot;. And this topic can be literally anything. It can be whatever name you want. And then you use that name to group all these different things.

What's happening here is that the foo process and the bar process both have these output channels. These output channels are being mixed under the hood into this one global channel. And then in your pipeline, you can say channel.topic. And then you give that topic and then it will give you essentially a queue channel.

This whole versions thing, which I'll show an example of in a second, is really the only example.

**Phil:** I've got one other, which is MultiQC.

**Ben:** That's the other one I was thinking of,

**Phil:** This is where you've got outputs from nearly every single process in the pipeline and you want to collect them all in at the end to one more process.

## Eval output

**Ben:** The other prong here is this new eval output. If you look at the nf-core processes, the modules, you'll see that every process has this epilogue where it's printing this whole like YAML information into a file. And that's, what's being emitted as the version.

And so this eval output is essentially a shortcut for that. What it allows you to do is to specify shell commands here, and then that shell command will essentially get appended to the end of your task here under the hood. It's essentially a convenience it's a bit easier to read to just have it up here.

You could do a similar thing with an env output. For example, you could have down here, you could assign the output of this command to an environment variable and then export that environment variable here. But again, this is just much simpler syntax wise.

**Phil:** This is what we do in nf-core. At the end of every single module, we have this thing where it's printing a yaml file into a file with this pretty ugly syntax.

**Ben:** You've got a script and a stub, so you also got to put it in the stub. And then also, in some cases, you might have conditional scripts. So you might have some if statements in this section and multiple different scripts.

And if you'll see examples of those in nf-core, that code is being duplicated across each conditional script. And so that's another nice thing about this eval output is that you specify it once. And then no matter what script you end up executing, it will be added.

## Eval + Topics together

**Ben:** Okay, so that's the eval output. And so, the final solution is how these two things come together. And so now what you have is you want to collect a version of all of the tools you've used across your entire pipeline. What you can do is instead of all that YAML boilerplate stuff, just emit an eval output here.

If you look closely at that YAML output, you see what also is being included is the name of the process and the name of the tool itself. And so you can just emit all this information as a tuple and then emit it to this channel topic.

And then down here, I collect this and this versions channel will contain all those output channels across the entire pipeline. And then you can write it to the file however you want.

Now something I didn't include here was the actual mapping of this tuple to the YAML syntax that nf-core uses.

But that would just be a map operator here so that's it.

**Phil:** It doesn't matter if you've got multiple tools within a container here, because you can have as many of those eval statements as you want.

And so if there was a second tool in this process, you could just duplicate that line and then you've got val fastqc, you could call it a different kind of name and then have a different eval statement.

**Ben:** Yeah, it would essentially be just another one of these lines. And that's another nice thing about using the topics, is that you could send both of them to the versions topic and they'll just be added.

You don't have to have a separate emit name or anything like that. It all just goes into the same place.

## Topics and eval in nf-core/modules

**Phil:** Something I wanted to flag at this point is that we're in the early stages of thinking about how to apply some of these new syntaxes across nf-core and nf-core modules.

As you can imagine, it's going to be a bit of a big lift to do this. There's quite a lot of editing code. We're going to automate as much as we can. And also we're going to try and batch together many of these changes in one go because

We're editing the process code Which is going to be in the nf-core modules repo and we're also editing the workflow at the bottom The collection of both channels so that's going to be in each pipeline So this is going to be a bit of a kind of a shift.

You're going to have to update all the nf-core modules. Then you're going to have to pull in the new versions of nf-core modules if you're a pipeline developer. And then you're going to have to update your nf-core template based Pipeline to work with the new syntax.

So it's going to be like a three step thing. we're starting to plan it out. There's a big old issue over on the nf-core modules repo over here, called a batch module updates. And it's got all the different stuff we want to do, all at once to try and rip off the bandaid. You can see here, we've got topic channels and using eval. So if you're really interested in how is this going to roll out, please jump in there. And add your thoughts and take a look at what we're planning. Because it's good to be prepared.

It's going to be very smooth this time around for nf-core there'll be a blog post about it. We're going to do a video walkthrough of updating and it's going to be, we're going to try and make it as simple as possible.

The end result is going to be much clearer syntax. So that should be better pipelines for everybody.

**Ben:** Once you figure out that process, that'll be good because I got more changes coming down the line. So get ready!

## Workflow output definition

**Ben:** Okay - this is the big one, the workflow output definition. We also sometimes call it the output DSL, but I think that's maybe a bit confusing because it's not a new DSL. It's just part of the existing DSL.

The way I think about the workflow outputs is essentially a fulfillment of the DSL2 paradigm, for publishing.

DSL2 was really all about modularity and having your channel logic being defined in workflows and being able to reuse those and all of that.

But we never really. publishDir into that workflow based paradigm, publishing was still being defined at the process level. it was mainly just really clunky. It's just a lot of repetitive code to do simple things.

And it also makes it hard to understand at the workflow level, what is this workflow producing? Because it's all hidden within these processes

nf-core got around that by having this giant modules file that has, all of the publishing defined in one place. It's not pretty to look at, as it's nested within these process selectors. You've got if statements in here because of how some processes are conditionally executed.

You can see you're repeating the same few things like a bajillion times, like the mode and other things like that. Ideally this stuff should, does not need to be nearly as complicated as it is.

So basically the workflow output definition was the result of us rethinking how to define publishing at the workflow level.

And so this snippet captures the main elements of this. The first thing we've added is this output block which goes alongside the entry workflow here. So wherever you have an entry workflow, you can also have an output block. And it essentially defines some configuration settings.

The main one here is the output directory. So this is the base directory.

And normally you would have to specify that as part of the publish path for every single publishdir definition.

**Phil:** You can see that as we've got a publishdir. path, and then you're prefixing it into the string with string interpolation every single time.

**Ben:** So now Nextflow has a built in concept of the base output directory. Now by default, it'll just be a launch directory like normal, and you don't actually have to publish everything to the same directory.

Basically this will be automatically prefixed to the final publishing path and then other published options like mode you can just specify it here and it will be applied across the board. Again, you don't have to specify it every single time.

Most of the publishDir options are supported here. You can look at docs to see which ones are actually supported. We tried to take this opportunity to simplify things a bit. So certain options don't necessarily apply into this paradigm. So we haven't supported them, but most of them are.

And then, in the workflow, we have this new section. So just like you have take, emit, and main, now there's a publish section. And this allows you to define what we call output targets or published targets. For channels. So instead of publishing files out of a process, what we do now is we publish a channel.

Any channel that's defined in this workflow, this could just be a channel you define here, like channel foo, or it could be an output channel of a process like foo or bar. If this is a sub workflow, sub workflow output, any channel that you can reference in this workflow, you can reference here.

We're using this arrow syntax to say, send this channel into this target, this published target. It's essentially a new operator that is only allowed in this published section.

And this target here is also a new concept. Now, this looks like a path, but this target is essentially a name, almost like topics. You can think of it in a similar way as topics. It's not a topic, but it's a name. And that allows you to map many, as many different channels you want to this name.

And then down here in the output block, you can customize the publishing behavior of that name. Including the path, like the actual path that name maps to, but if you don't do any configuration, then it will default to just being the published path.

And so then you just extend that across your entire pipeline. What you essentially have is a bunch of channels here and you all in this one place, you can specify where everything is being published.

And so this is really nice because now you have everything in one place, all the publishing in one place.

**Phil:** I'm just going to walk through a couple of these things here just to try and make it clear because there's some new concepts we're talking about here. So the publish block is new in the workflow. We've got an output channel, which could be anything. What happens if it's a string or something and it's not a file?

**Ben:** So this is the section, by the way in the docs, if you go to the workflows page, then this is the last section on publishing outputs.

So let's say your channel looks like this. This is a typical structure you've got, a map of metadata. And then you've got some files and this is all in a tuple.

So what will happen here is that when you send this channel into a publish target, when this value is published, Nextflow will scan the tuple for any files that it finds, will be published.

And that's essentially how it works. So any files like in this tuple structure. I don't think the map is scanned. So if you have a file like as a value in a map, that's not a very typical pattern. But put it in the tuple, it'll be captured. And of course, if it's just a file by itself, it'll be captured.

**Phil:** Okay, so values are just effectively ignored. So if I did publish the channel where it was just a simple value, just nothing would happen.

Okay. So I've got published my channel here, which is going to, to foo. And if I leave it as it is, then it's going to be published in the directory from the launch directory called results, slash foo.

And if I want to customize that path instead of foo...

**Ben:** So you type the name of the target like that, and then you add some curly braces.

And so this is like a new block. You can think of it as like a target block. So you can say path and then whatever path you want. This is foo by default because the name is foo.

But yeah, you can change it. And so now it would be results/bar. And so this is a relative path against the output directory. And then all the same publishing options that are available in the output block are also available there. Which would allow you to, configure individual targets as you see fit.

**Phil:** And this output block here, does this have to be edited in the pipeline script or can configs as well?

**Ben:** Right now it's only something you can do in the pipeline code. I personally would like to move much of that into the config, because I feel that things like the output directory and the mode are really configuration settings and not really things that are like inherent to the pipeline.

So currently, if you want it to make those settings configurable, like by the user, you would have to make a param and then use that param in the output block, which is essentially what nf-core is already doing with things like params.outdir and params.publish_dir_mode.

**Phil:** params.Outdir and then it's going to be the same behavior as we have in nf-core pipelines.

**Ben:** I'm thinking maybe we could skip a step and just have, maybe an output config scope, or maybe it's a workflow.output config scope, and then you can just specify the directory there. And then you can still have a param if you want and map it into that config setting, but at least that way you're not forced to.

And similarly with all of those options in the output block, if the user, because there's a bunch of them there's stuff around like tags and content type and different things like that. And right now you would have to expose a param for each one of those settings in order for the user to change them without just cloning the pipeline. So I think that's something we would like to do.

**Phil:** I could imagine pipelines getting very bloated if they have to stick in every possible configurable option within the output block and then set every single one to a parameter just in case the end user wants to configure it.

**Ben:** Yeah, it seems unnecessary. And also, I feel like it's an anti pattern to have params that don't pertain directly to the workflow itself. Having these params that are controlling things that are really outside of the pipeline code, I think is something we would like to get away from.

I think in the longterm, maybe those things should just be configuration options. Maybe the outputdir should just be a built in CLI option. So yeah, there's a lot of things to think about with that.

**Phil:** Where should people go with, if they have ideas or thoughts about this?

**Ben:** I guess the easiest thing would be to just submit a Nextflow issue. Of course you can ping us on Nextflow Slack, but I think GitHub issue would probably be the best way to do that

**Phil:** So we've had this syntax has been released in Nextflow for about a month or so Has there been any kind of suggestions or feedback that's come in already?

**Ben:** I think you know what i'm going to say here because when we were working on this feature, Phil said _&quot;You know, people are going to want to do like the whole dynamic output thing&quot;._

And I was like we'll let them come to me if they have that. And then immediately we release it. And the first thing people start asking is, _&quot;can you still do the dynamic published path based on the sample ID and stuff like that?&quot;_ And we were like,

**Phil:** Told you so!,

**Ben:** but it was important for me to, see that community request.

Yeah, the main thing that's missing right now is, with publishDir, you can use saveAs to have a dynamic published path. You can imagine the typical pattern is using the sample ID and you want to have some subdirectory for each sample ID. We didn't support that at first. And now that people have been requesting it, I think I've found a good way to do it. So instead of just adding that saveAs option, when you have that path option in the output target block, will allow you to make that just a closure.

Let's show an example of this. You have the target block and then you do path. So let's say that the value coming in was that structure that I showed before with the metamap and in the file. So what you would do is say something like meta comma, fastq.

Nf-core, this should be familiar this is a common structure. You've got the metamap and then this is just a tuple. Then you would be able to use any of those values to construct the dynamic output path. So you'd have a dynamic string and say, $ braces meta.id slash, whatever.

And this is slightly different because it's not giving you all of the variables that would be available in the process definition, but it's giving you the variables that are available in the channel.

I noticed you were starting to type something like task dot, whatever, if you wanted to use the task variable, like the process name, you would need to include that in the output channel, and then you could use it down here. So that would be an important caveat.

**Phil:** So this will be I guess look out for an edge release in the coming months.

**Ben:** Yeah, so some of the main questions we're seeking feedback on, obviously that dynamic output is one. The config versus pipeline code issue,

do you think it would be useful to configure those individual targets in the config file? So having something like a with target, like we do with like with label for processes. Being able to configure individual targets in the config file and not just the overall config options.

And then also just the use of the publish section in the workflow. One thing I didn't mention was that you can also have a publish section in the process, which is something that we tried to see if we could have things more modular in a sense. But I'm going back and forth on that thinking maybe it's not, A good thing to have the publish section in the process, because what happens then, it's more modular, but then you lose that quality of all of the publishing being defined in one place.

So one thing I'm thinking is maybe it's better to only allow it in the workflow and then that forces you to put everything in the workflow and have it all in one place.

We've given you the ability to put it in the process as well so that you can experiment with it. we want to see what people come up with what edges they run into, what walls they run into. This is your opportunity to propose new syntax or new approach.

**Phil:** So on that last point, nf-core is another good example here because the modules and the pipelines are so separate from one another. So if it can be done at process level, then within the shared nf-core modules, then we can define output, publish, directories in the shared module, so you'd import the module and it would just magically work.

But if it's only a pipeline level, when you have to import the module, it would run, but it wouldn't create any files by default. So after you've imported the module, then you have to define the channels that you want to publish.

I guess the benefit of that then is you have to be a big publish block, which could be very long, but it's going to just be a big list of all the different channels, which the whole pipeline is scooping up and publishing.

**Ben:** Yeah, and there are some shorthands. If you publish the .out of a process or workflow, then it will implicitly publish all of the channels in that .out . So if you structure your workflows a certain way, you can make that publish section much shorter than otherwise.

The topic channels feeds into this because , if you have that extra versions channel on every process and sub workflow, you can't really do that because you don't want to publish those version channels.

But if you move all of that into a topic channel, maybe that's increases the chance that you actually just do want to publish all of the outputs. But even then there's cases where maybe some of the channel outputs are intermediates and some of them are finalized.

This happens to work where you could also publish the .out somewhere and that's a default, you can overwrite specific output channels. So if you wanted to disable certain output channels, you could then specify those afterwards and send them to null. if you publish to null, then it just won't publish it. That gives you the ability to, control the publishing with a param.

**Phil:** So if I'm looking at the the output block, I have foo, I can set the path. I do params.publish, which is the variable. If it's negative or it's false and I set it to null, which then means that the path is null, which means it doesn't get published.

**Ben:** So yeah, that works. And you can also put it in the workflow publish as well. Just grab that ternary expression and plop it into where foo is. That's the syntax for having a conditional publishing, which is also a common pattern in nf-core, where you have lots of intermediate outputs, but you have params to allow people to publish those intermediates so that they can inspect them for debugging and things like that.

**Phil:** And the other thing you said just before that is if I imagine I had a big named sub workflows here and those sub workflows had big emit blocks where they're emitting 20 channels.

I could also have my sub workflow .out to my path. And so that would scoop up every output channel that's emitted from a specific sub workflow.

So if you're smart with grouping stuff by different sub workflows, it makes the final workflow publish block a lot less verbose.

We're publishing everything to mysubrightflow except this specific output, which we're publishing to null, so it doesn't get published.

**Ben:** Can we also show an example here of how the process publish works.

So when you have an output channel like that with an emit name, you can reference that emit name in the publish section. So now I can say bar and then the double arrow.

This is a abstract target, because the process isn't part of any pipeline. So as soon as you call this foo process in here Then that publishing behavior in foo is used as a default and it will be published to my output bar

**Phil:** And I can get rid of the workflow publish block completely.

**Ben:** Yeah. So in this case you wouldn't need it at all, unless you wanted to override the default.

So now if you go into the workflow and add that publish section back in, and let's say, let's publish foo.bar to something else, and so what happens here is you're overriding the default from the process. If this was a sub workflow then you could call that sub workflow in a larger workflow and overwrite the sub workflows. Override so it's this cascading thing.

This is one of the things that I'm not completely sold on and I'm interested on hearing more from the community about this because I feel that the overriding, while, it allows you to separate everything into the processes and not have this giant published section. I'm not convinced that's actually a good thing. I feel like maybe that actually makes it harder to read. Because now, the publishing is split up across all these different components in your pipeline, even though it makes them more reusable. It's split up and now you have to think about how do all these overrides get resolved.

It's the same problem as when you think about process directives, the fact that, process directives can be defined in all these different places. You can have the process definition, but then also in the config file, and you might have many different conflict files being merged in, and then you could have just generic process config, but then you could also be using withName and withLabel.

One of the things we developed was this Nextflow inspect command just to figure out, like, how are all these things being resolved?

We're thinking about doing the same thing for the output to say, how are all these outputs resolved? And, the question that I keep coming back to is maybe we shouldn't even need to do all that in the first place. Maybe we should just make things simpler, just bite the bullet and enumerate all the channels that you want to publish in that workflow.

This is, that's just me. I'd love to hear what other people think about this.

**Phil:** I'm very on the fence on this one. I quite like being explicit about it. One of the complaints we get with nf-core pipelines is that you're chasing this path the whole time trying to track things down, trying to work out where this is being defined.

I feel if we're doing this at process level, it's a bit of that, even if it looks clean, it's actually a pain and it's maybe better to be verbose and stick everything at the end.

On the other hand, I can also see it being annoying that you write a new process, forget to add it into that publish block. Run your pipeline, and then you don't notice that your files just aren't being published until later on. Which would be irritating.

**Ben:** The broader question really is about modularity and reusability, how to make these modules as reusable as possible. Without encoding all this pipeline specific information into the module prematurely.

Think about params, like it's a bit confusing for a process to use a param directly. Because a process isn't attached to any particular pipeline, which means it's not attached to any particular set of params. When you're using a param, you don't really know anything. That param could become anything, depending on what pipeline you plug it into.

But at the same time it's nice to have some way to show: if I incorporate this process into a pipeline, what's an example usage? How would I actually, connect this process to the params on one side and the publishing on the other side.

One idea that I had was maybe in the module, you could have an entry workflow, which gives an example of how to use that. A simple workflow wrapper around the process, and you could include some example params and some example publish targets.

That serves a lot of nice purposes it serves as documentation. So when the user imports the pipeline, they can, think critically about the example and think about which params and publishing you want to incorporate.

But at least that example is there for you to copy and paste and then edit to your needs. It also provides a test case. So it's almost similar to nf-test: you could then run that process directly and you could provide some test input data or whatever, and then it would just work out of the box.

So when you import the process as a module, that workflow would just be ignored. Because in that case, you're not using it as an entry point.

**Ben:** So this is very similar to how Python modules work in your Python module, you can define a, if name equals `__main__` and that defines what this script should execute. If you execute it directly as like your main script.

That allows you to still import the script as a module and then that code doesn't get run. I really like that pattern. That's the way I'm thinking about how to, provide that documentation without having all this params and publishing stuff directly in the process definition.

**Phil:** Another idea that springs to mind while we're talking, which is not mutually exclusive with your suggestion of having workflows at module level. We could have Nextflow basically fail with an error if there's any file outputs coming out the main workflow which don't have a target publish declaration, so then everything has to either be given a publish target directory or be assigned to null. And then if you added a new module and forgot, it would just crash immediately.

**Ben:** I think that's still something where the Nextflow inspect command can be useful. So even if we, do the stuff we've been talking about making the output definition simpler. It would still be useful to have the inspect command just to show you; based on all the published definitions I've set up, what is the final list of things that get published. But yeah, it would still be useful to at least tell the user if something is not being published, even if it's just a warning.

**Phil:** The take-home message for anyone listening then is to go and read the blog post on the NextFlow website, look at the docs for the new workflow output publishing definition, and try it.

Build a workflow, edit an nf-core pipeline, mess around with it, try and break it, try and see what features are missing. Find what syntax is frustrating, find what documentation isn't clear, and then come and tell us about it.

We're going to stamp this as official come Barcelona, in October. So you don't have long to do it. You've got a few months to go and play with this syntax and get all your ideas back. We're looking forward to seeing what people do with this.

Fantastic. Ben, thanks so much! As always, it's been an absolute pleasure. Thanks very much for listening, everyone. We'll be back with you very soon.

**Ben:** Bye everybody.
