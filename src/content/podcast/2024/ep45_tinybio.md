---
title: "tinybio: Nextflow x AI"
episode: 45
description: GenAI for bioinformatics
date: 2024-09-03
type: podcast
subtype: Technical discussion
youtubeid: Qyu8d-Mq5-Q
image: /img/podcast_ep45.png
tags: Nextflow, Seqera
author: Developer advocates
icon: logo_podcast_channels.jpg
---

In episode 45 of the Channels podcast, [Sasha Dagayev](https://www.linkedin.com/in/sasha-dagayev/) from Seqera (formerly [tinybio](https://www.tinybio.cloud/)) shares his journey from working in tech and consumer spaces to entering bioinformatics during the pandemic. Sasha discusses the challenges and innovations behind tinybio, including the evolution towards integrating Generative AI (ChatGPT) to simplify bioinformatics tasks. Key points include:

- Sasha's transition from tech to bioinformatics post-pandemic.
- Initial challenges and pivot towards Generative AI for tinybio.
- Not just ChatGPT: code execution and iteration.
- AI concepts: Foundational models, RAG, Langchain and fine tuning.
- What makes Seqera and GenAI the perfect match.

Stay tuned for more exciting updates soon at the upcoming Nextflow Summit!

- [tinybio](https://www.tinybio.cloud/)
- [Blog: Seqera acquires tinybio to Advance Science for Everyone - Now Through GenAI!](https://seqera.io/blog/tinybio-joins-seqera-to-advance-science-for-everyone-now-through-genai/)
- [Nextflow Summit 2024 Barcelona](https://summit.nextflow.io/2024/barcelona/)

# Podcast overview

## Introduction

Artificial intelligence (AI) is affecting everyone, and bioinformatics is no exception. In this episode, we explore the pioneering work of [tinybio](https://www.tinybio.cloud/). Host [Phil Ewels](https://x.com/tallphil) sat down with [Sasha Dagayev](https://www.linkedin.com/in/sasha-dagayev/) from tinybio, now part of Seqera, to discuss the journey, challenges, and future implications of this merger.

## Sasha's background

Sasha's journey into bioinformatics began after a significant pivot from the tech industry. Initially immersed in consumer tech, Sasha’s transition was catalyzed by the COVID-19 pandemic, propelling him back to academia to study bioinformatics. He quickly recognized the ample opportunities in this burgeoning field and committed fully to it.

Unlike the traditional path of computational science, Sasha's background is in business strategy. His early career as a McKinsey consultant and later at the startup Zocdoc sparked his passion for the tech startup ecosystem. Despite the unconventional start, Sasha self-taught a large part of what he now knows, eventually immersing himself in lab work, where he saw both the potential and shortcomings of bioinformatics.

## The Birth of tinybio

tinybio was born out of Sasha's frustration with the difficulty of running bioinformatics tools. Unlike consumer tech, creating a proof of concept in bioinformatics was challenging. Initially envisioned as an API, tinybio soon pivoted toward leveraging generative AI, significantly after the release of ChatGPT in late 2022.

Targeting efficiency, tinybio developed a chat interface that serves the bioinformatics community by using powerful AI models to generate and execute scripts. This involves not only providing information but also running actual bioinformatics pipelines on a high-performance instance, significantly improving user experience through automation.

Sasha recounts the early partnership with co-founder Vishal Patel, who helped lay the technical foundations of tinybio. Their combined effort culminated in the development of tools that initially summarized documents but evolved significantly with the capabilities of generative AI models like GPT-3 and GPT-4.

tinybio's public launch was marked by high expectations and significant feedback. Despite initial setbacks with features like citing sources — which turned out to be less valuable than anticipated — the platform rapidly iterated to find its footing with features that mattered more to users, like executing scripts directly on integrated computes.

## Integration with Seqera and Future Prospects

Sasha expresses excitement for the future, illustrating how tinybio's integration with Seqera amplifies the utility of both platforms. Seqera’s extensive context management capabilities for bioinformatics workflows complement tinybio’s AI-powered script generation and execution.

A priority for this integration is making Seqera’s advanced features accessible to non-bioinformaticians, such as wet lab scientists, by enabling them to query and interact with bioinformatics data intuitively and effectively. The goal is to democratize access and enhance usability across diverse user groups.

Stay tuned for updates and announcements during the Nextflow Summit in Barcelona and the broader community activities in Boston. This collaboration signifies a transformative step in how bioinformatics workflows are managed and executed, promising to bring about significant advancements in the field.

Join us on this journey into the future of bioinformatics, where AI and seamless workflow integration pave the way for groundbreaking discoveries.

# Full transcript

:::note{title=Transcript}

## Introduction

**Phil Ewels:** Hello and welcome to Channels, the Nextflow podcast. You're joining us for episode 45, going out on September 3rd, 2024. And today I'm really excited to talk all about Nextflow and AI. To do that, we've got Sasha who's joining us from tinybio, but actually very recently now from Seqera.

**Sasha Dagayev:** That's right. We're proudly acquired by Seqera joining the open science mission here.

## Sasha's background

**Phil Ewels:** Fantastic. We've got a lot to talk about today, butlet's start off Sasha by telling us a little bit about who you are and what your background is and how you got into all this in the first place.

**Sasha Dagayev:** Yeah, so I was working in tech for about 11 years, but mostly in the consumer space a lot on the search landscape here in New York City. And I was very happy with it.

And then in 2010, the pandemic happened and II was in a client call and a client was so oblivious to what was happening right outside my window. It felt like an alternate reality. And I was like, okay, something has to change.

So I went back to school to study something different. And I saw bioinformatics and I was like, okay, this is blue ocean field. All I saw was opportunity.

It really reminded me of the search space in the early 2010s where. The people on the inside think it's like, &quot;Hey, it's like over&quot;, while the people on the outside are like, &quot;Oh my God, there's so much to do here&quot;.

So I went both feet in deep side of the pool, started working at a lab. And trying to get more of an understanding of it.

The more I worked with it, the more I saw things that were really good about it and, some of the things that were lacking for lack of a better word.

**Phil Ewels:** So what was your training originally? Did you studied like computational science or something?

**Sasha Dagayev:** No. So I went to school for business strategy. When I was 18, I thought that traveling as a McKinsey consultant was like the best job ever. So I was going to be a,management consultant. That's what I was trying to become.

Then I. Got a job at Zocdoc, which is a online booking service here in New York City that does doctor's appointments, and I was a very early employee there. And I got the startup bug. So self taught myself a lot of the stuff and then just kept going with it.

**Phil Ewels:** Amazing. And when you say you joined the lab, was that a bioinformatics lab or...

**Sasha Dagayev:** Yeah, it was a wet lab that needed a bioinformatician. So it was like the classic church and state type of relationship where, they have a bunch of very smart lab people, but, they have a Zoolander problem, as I like to recall it; where the data is in the computer and they need somebody's help to get it out. And I was trying to be that guy.

**Phil Ewels:** Amazing. And how long did you do that for?

**Sasha Dagayev:** That was about two years. So I was doing it on and off between that and school. I was very fortunate to have sold my last business as well. So I had an opportunity to really Do many things without having to worry about a specific job per se for a while.

So I was doing training both in school and in this lab at the same time and starting to think about the next business, which came out to be tinybio.

## How tinybio started

**Phil Ewels:** Awesome. So, tinybio, tell us about it.

**Sasha Dagayev:** So the way that it started was that when I first started bioinformatics, I was shocked by how hard it was to actually get something to run, just to try it. My previous experience has been a lot with front end technologies like consumer tech. And so if you've ever coded, a consumer application with react or, whatever front end framework, it's so easy to just see it: a proof of concept, right? A proof of concept is so simple to get to.

With bioinformatics, that was really difficult. It was very difficult to just, &quot;Hey, you see a new aligner&quot;. &quot; I'm reading that it works, but I just want to see that it actually does work.&quot; And I saw that just didn't exist.

So I thought from a dev quality of life perspective, this is a service that needed to exist. And initially it was going to be an API, but then November 30th 2022 happened. ChatGPT came out and I was like, okay, like that's now a terrible idea. Like now we're. 100 percent pivoting to generative AI. This is definitely going to be part of the stack that somebody's going to be doing, especially for bioinformatics.

**Phil Ewels:** Cool. So how much have you done by that point? And what did you have to throw away?

**Sasha Dagayev:** So yeah, we had a functional RNA-seq endpoint that you could throw arbitrary data into . To be honest, I thought that, some of the idea was a little bit underbaked, cause the customer research that I was doing was very underwhelming. Like the customer feedback that there was like, &quot;Oh, cool&quot;. Nobody was asking &quot;Hey, when can I try it more?&quot;. So you get polite responses, of people not wanting to hurt your feelings.

But then when we came out with our first demos for tinybio, which was essentially, a service.

If you guys remember when chat GPT came out, the initial demos that were coming out , everything was: &quot; Chat your PDF&quot;. The flow was, Hey, somebody uploads a PDF and then they ask questions against that document.

For us, what we did was we did a slightly different demo, we said, okay what if we did chat? What if we just gave you what a bioinformatician would care about from a paper? So somebody would upload a paper and then we would use GPT three at the time to break it out into blocks, find and do look ups against it to understand Hey, which parts of this paper are actually interesting and which parts aren't.

And when words started to appear on the screen during the demos, people were like, Oh my God. So I felt that I was like, okay, this is clearly the way to go. You know what I mean? The conversations are going much better than they were going before.

## tinybio co-founder: Vishal Patel

**Phil Ewels:** And, it wasn't just you by yourself. You're co founder of tinybio: who else is involved?

**Sasha Dagayev:** Yeah, it's my co founder Vishal Patel. We worked together at a company called Noteworth. This was back in 2014, 2015. And the idea back then was to help folks that were diabetic to become pre diabetic, or if you were pre diabetic to become Not diabetic at all. Using connected devices.

That idea sounds great -did not take off. There was a bunch of reasons why it didn't, but like we stayed in touch and, when I needed some help with configuration of cloud instances and stuff like that, I called up Vishal and then we took off from there.

## tinybio chat public launch

**Phil Ewels:** So then your first productwas to summarize documents basically. Is that right?

**Sasha Dagayev:** Yes. that was the first product. And it got a little bit, like we had really nice traction with VCs. We didn't have a public launch, but we had some nice traction with venture capitalists and being like, Hey, we can build something on this.

And then there was this moment. At that time when essentially everybody was saying that, yeah. Everybody's a Chat GPT wrapper. And they didn't give any credit to the application code around it, but we just kept going. So we were like, okay, we're going to just launch no matter what.

And so on June 6th, 2023, we launched the chat service, for free without any hard caps against it, just to see: Hey, will this thing get any traction whatsoever? And I posted on Biostars and the reception was fantastic initially. There was definitely some feedback that people were getting hallucinated responses, but I think, if they're not totally dumping on you on Biostars, I'll take that as a win for that forum.

And yeah, we were very happy. And that was like the takeoff point for our chat application.

**Phil Ewels:** I feel like there's a healthy dose of skepticism within the bioinformatics community . There's been a definite kind of hype bubble. So to get a) not ignored and b) to not be taken down too hard I think like you say it's probably considered a success

**Sasha Dagayev:** Exactly. So we had a chat service that essentially was like the initial version of the product and at the time chat GPT was running GPT 3. 5 and it did not have default access to the internet. We had Retrieval Augmented Generation pipelines, against a bunch of really popular curated GitHub repos.

And when we tested us versus Chat GPT, Chat GPT native at the time only scored about 16 percent against bioinformatics questions. To be fair, these are like hard bioinformatics questions, things like, &quot;Hey, what is this specific parameter in STAR do?&quot; What is the difference in the ways that, Edge R versus DESeq work in terms of showing your differential expression analysis.&quot; So not out of the box, but things that you would ask on BioStars, which was my kind of bar.

And so vanilla Chat GPT answered at around like a 16 percent and we with our RAG pipeline got a 69. And so I was like, okay, this is much better. Like the quality of the generation is much better. And so we went into that pathway.

And that really worked for us for, like three to four months. We had a bunch of features that we were launching at the time.

## Failed feature launches

**Sasha Dagayev:** But there's like a bunch of failed things that like people just didn't care about that we launched at the time.

**Phil Ewels:** Any examples of those?

**Sasha Dagayev:** Yeah, so we added sources to the generations. Because one of the big criticisms is everybody was like, Oh, there's a hallucination. Like I need to know where the data is from, blah, blah, blah, blah, blah. And we were like, we sunk two months into making very nice sourcing and all this stuff.

At the time it was very difficult to actually display the information properly, and no one cared. Absolutely no one cared, it didn't get any more usage. I posted that Hey, we have sources now. No one cared.

It's like one of those things where it's &quot; do as your customer does, not as your customer says&quot;, type of situations where everybody was just pasting stuff in and using it and if it worked then nobody really actually was that gung ho about the sources. Because one of the things that we realized is that people were trying to solve real problems, nobody was asking meaning of the universe questions.

So do you really care? About the source if I troubleshoot your environment issue correctly?

**Phil Ewels:** What people really want is really good answers. And what they're asking for is sources.

**Sasha Dagayev:** Yeah. So it was very humbling, but a good learning experience to see it. So then our usage flatlined a little bit. And we knew we needed to make a pivot.

## Running the code too

**Sasha Dagayev:** In order to get on the radar and have a bigger release, we knew we needed to not just write this stuff, but we needed to actually execute it.

And so that's when we did this version that you'll see today, where we not only write the pipelines, we can actually run them as well.

**Phil Ewels:** Because it's not just a chat client like GPT. You've got a sidebar there.

**Sasha Dagayev:** So we have this chat client here that looks very much like your typical chat GPT or Gemini or whatever chat you use. But the killer aspect of it is that it's actually hooked up to a high performance, instance that actually has all of the bioinformatics hits written on it.

Essentially gives a GPT access to a compute environment and actually, writes and runs the code.

The cool part about it as you can see today and for those watching the demo is that: it does hallucinate, still even quite a bit on the frontier models. But it self corrects and edits the queries and runs them again and then gives you the proper answers.

**Phil Ewels:** So it's writing a script, running the script, getting an error message, understanding what the error message is, correcting the script itself, rerunning it, and keeps doing that until it gets an answer, basically.

**Sasha Dagayev:** Yeah. And so we can see here, it generated a sample to sample distances, heat map It essentially booted up a instance, loaded up these libraries, and then actually ran this for you.

Our biggest assumption from a competitive perspective was that OpenAI and Google, they're at least two to three years away from letting somebody install DESeq on their machines. So we thought This is a good place to go.

And then once we released this, we saw pretty good public reception, on LinkedIn, on socials. we saw a very high increase in retention rate users that interacted with this feature versus the information only feature. so we knew this was the right thing to go after.

## User trust

**Phil Ewels:** We talked about user trust already with the sources stuff, but you click for a second and brought up an R script. That's something else I really like about this, you can actually see exactly what it's doing and presumably I could grab that R script and run it myself locally.

**Sasha Dagayev:** Exactly. And this was like, our biggest thing: we need to actually show this information as quickly as possible to the user and realistically, like you're not going to do your dissertation using tinybio, but you could really quickly test stuff out.

We are by far the fastest place for you to understand &quot;Hey, does this analysis kind of work right?&quot; Cause it will figure stuff out iteratively all on its own.

**Phil Ewels:** And it gives you a starting point because presumably now I can carry on editing that script and extend it and do whatever else I want to do.

**Sasha Dagayev:** Yeah. We knew that this was going to be the applications killer feature.

And so I reached out to a couple of connections at Seqera cause, I knew that because of NextFlow, Seqera has, a level of trust from the community in terms of actually running the thing.

And I was like, okay, if I'm going to be running a lot of these I'm going to work with people that really understand how to do the thing. Initially I approached, I was like, Hey, if you guys could help us with distribution, we could work something out.

And then I think the conversations become more and more involved. And that's how Vishal and I decided to join the Seqera team.

## Accessing files via chat

**Phil Ewels:** A couple of questions about the chat interface. One, we had the files on the sidebar there. Where are those files? There's an upload button. This isn't just sitting on the instance or is it? How does that work?

**Sasha Dagayev:** Yeah. So this uses Fuse technology. So very similar to Fusion. For each tinybio user, we spin up an S3 bucket, and then we FUSE it to the instance at the time. We save a ton on compute costs, because you don't have to keep the instance running. You can power it up, power it down. And this gives the user some level of continuity in terms of knowing that the data is there.

The longer term idea with this, if we worked with a smaller lab or with somebody. Is that they could bring in their own S3 buckets and we could FUSE them to our compute environment.

Because what was special about tinybio was that nice integration.

## Experiment idea generator

**Phil Ewels:** And when you went in at the start there, there were three different tabs Chat was only one of those tabs. Can you tell us a little bit about what those other kind of branches you took were?

**Sasha Dagayev:** This one is an experiment generator. And essentially what I saw was that a lot of timeswe had access to actually more samples than we can imagine, and we imagine the future where, we had autonomous GPTs that were constantly querying our data.

And I just built this as a proof of concept, for investors and for students to understand, Hey, what can we do with this? And the basic idea was that you would describe your data and it would give you some ideas of what you could run. So in this case, we said, Hey, we have some single cell RNA seq samples of neural stem cells. These samples were collected from a specific region at different developmental stages, and then it gives you some ideas of what you could actually do with these data types. So it says, Hey, you could do differential expression analysis, trajectory analysis, cell type identification, gene co expression network analysis pseudo time analysis and pathway enrichment analysis.

And then for each one, you could hit generate code. And it would actually generate the code for you based on the software packages that you had put in.

And so this was an idea that we had for the step one of an autonomous GPT, which essentially looks at a data repository and then continues to iterate on itwithout break, essentially.

**Phil Ewels:** I can imagine this being quite valuable for interdisciplinary teams as well or anyone coming into bioinformatics who's not already super familiar with what kinds of analysis exist. Because I mean there's a whole bunch of different options there

**Sasha Dagayev:** Yeah. It's very funny, it's when we posted about this to biostars, everybody was like, &quot;Oh, this doesn't need to exist because everybody knows exactly what they want to work on.&quot; And I was like maybe you don't know everything, and so this is like one of those, future ideas that I had.

**Phil Ewels:** I think that's probably like a biased demographic because if you're signed up and reading Biostars then that's true

## Pipeline generator

**Sasha Dagayev:** Exactly. So yeahthis was one of the ideas and the other one was this pipeline generator, which was essentiallyto generate pipelines as quickly as possible based on specific inputs.

Essentially ,the easiest way is that it's just for somebody that one's like a very simple pipeline to do very simple analysis.to get stuff out and to limit the scope of the products.

This tool had some usage and had some product market fit. The good thing about the chat interface is that it allows you to iterate on things. This interface here does not really allow you to iterate. It's just Hey, here's the answer. And you can't adjust it in the interface itself. But the basic idea was to just get you a pipeline as quickly as possible.

**Phil Ewels:** And you had multiple different languages there. The one you're showing here is Slurm and sbatch,

I think. It's quite a structured input form here. What led you to choose these things and narrow it down in this way

**Sasha Dagayev:** Honestly, like the good thing about having the chat interface is that it allowed us to understand how popular something was actually versus how perceived popular something was.

So I think, if you're keeping up with bioRxiv or with what, people are posting about on Twitter. It's all, protein language models and spatial transcriptomics and all that stuff. When you actually talk to people, like what they actually do day to day, it's good old bulk RNA. It's a good old bulk RNA that like, people are still running it, it's still the King Kong of our space.

And yeah, so it's we knew that hey, you want to do, certain specific experiment types. You're going to want to do certain specific workflow managers. When I looked at a lot of enterprise software, I found it to be very overwhelming, especially for a non bioinformatician.

And I was like, okay, we can make something simpler. For people to just get their feet wet and yeah, make it much more Fisher Price.

**Phil Ewels:** I like that .

Maybe a bit of a mean question, but i'm going to ask it. In your experience. How often do these pipelines it generates work?

**Sasha Dagayev:** We don't know, the answer to that question. The ones that I've tested, they've worked. But there's, again, this is like one of the things that we never saw the actual feedback about this.

I do know that when we generate the pipelines in the chat interface and it runs them, It has about a 80 percent success rate in terms of like actually executing.

The biggest thing you have to understand is that in my chat interface you're executing in my environment. I know everything about that environment. So the LLM has perfect information every time about everything that is generating. If you take the script to an HPC, you don't know all the information about your environment. And so it's very difficult to troubleshoot some of these issues.

**Phil Ewels:** That's something I'm very familiar with. The tack we take with Nextflow and Seqera is, we run anywhere and everywhere, which is great. But yeah, you do often end up with what works on my machine.

**Sasha Dagayev:** Exactly. Yeah. It works on local, right?

**Phil Ewels:** Yeah, exactly.

## AI concepts: Foundational models

**Phil Ewels:** Maybe not everyone in the audience is super familiar with AI and we've touched a couple of concepts here: we mentioned RAG earlier, can you take us through just a quick &quot;AI for dummies&quot; about what this is all built on? What the foundations are?

**Sasha Dagayev:** Yeah. So the foundational models themselves are built on a scrape of the internet. So there's this concept for AI called common crawl, which is essentially a crawl of the highest value portions of the internet. So that includes things like PubMed, bioRxiv, Nature papers, things like that.

And those things are what's used for the foundational models, like GPT 3, GPT 4, Lama, Gemini, Clog.

Those models are increasingly getting bigger. Now what that means is that instead of optimizing certain set of neurons inside the model they're optimizing it to such a point where, you have I think the latest LAMA model is 400, 400 billion.

And so it's a language model that can optimize around a very high amount of complexity. And what OpenAI was wise to before everybody else was that actually if you just keep giving it more data, it actually just keeps getting smarter, right? And so that's what the foundational models are really trained on and are really good at.

The reason why they're really good at some of the stuff like bioinformatics is because the data that is used for the inputs that is most valuable is typically open academic journals. They are by far the most trusted sources of information for these models, which is why it does so well on open software and on open scientific concepts. That's why it does such a good job.

**Phil Ewels:** Yet another reason to publish in open access journals.

**Sasha Dagayev:** Yeah. If you want the models to know what you're working on, make sure it's as open as possible. That's why bioinformatics was such a good use case and had much less hallucination than other use cases because it's very well trained on that.

## AI concepts: RAG

**Sasha Dagayev:** The other two things that everyone kind of talks about. Number one is RAG. And what RAG is essentially: showing the model a document alongside a question. So instead of asking, Hey, can you write for me an RNA seq pipeline? You say, Hey, can you write for me an RNA seq pipeline? And you also, at the same time, show it the nf-core/rnaseq pipe.

And so you're giving it the answer and you're using the model's ability to reason. It's essentially you're trying to give it an open textbook quiz, and so that's the beauty of RAG. And that's why a lot of enterprise and a lot of for profit people are using this application.

**Phil Ewels:** That sounds great in concept, but how does the tinybio chat interface, know which things to give to the foundational model?

**Sasha Dagayev:** Yeah. So it uses a concept called cosine similarity. What it does is that it takes all of your collection of documents and it turns them into a giant matrix. So foundationally, what that means is, you know how we can describe colors with RGB. This is this very similar concept, but instead of three numbers, it's four thousand numbers.

So you're describing a document across a four thousand point matrix. And you say, Hey, the question that was asked, which is the closest document to this question? And then you present that document at that specific time. So that's the secret sauce that makes all of these RAG applications work.

**Phil Ewels:** So this is a kind of predefined data set that RAG can work with. It's not going to just the whole open internet in this case. It's going to like a subset.

**Sasha Dagayev:** Yeah. So folks likepopular generative AI applications like perplexity Gemini, they do to for lack of a better word, boil the ocean. They're trying to do this magic trick on billions and billions of documents.

For tinybio, what we found is that there's really like forty thousand core documents that you need that answer most of the questions. By shrinking the size of the search space, you get a much better way of getting to what you need.

**Phil Ewels:** Faster and more accurate, presumably.

**Sasha Dagayev:** Yeah.

**Phil Ewels:** Nice.

## AI concepts: Langchain

**Phil Ewels:** So we talked about models, we talked about RAG tell us about Langchain and fine tuning.

**Sasha Dagayev:** Yeah. So Lang chain is a fantastic framework for anyone trying to get into generative AI to start. The reason is because they have great examples and great plugins for people to use. They have both a Python pip application and a JavaScript NPM package.

So it's very accessible. And it has great getting started packages. The other benefit of using them is that you get access to a whole slew of modules that people have built on top of their applications.

**Phil Ewels:** So what does it do?

**Sasha Dagayev:** Yeah, so what it does is it's like Ruby on rails for folks that, did web development. It obfuscates a lot of the complexity of interacting with a lot of these large language models directly. It makes it very simple to get started and get a RAG application off the ground. So what would be, a hundred lines of code in Langchain is probably only 10 lines of code.

Now you are giving up some control, but you're off the ground much, much faster.

**Phil Ewels:** So for another analogy, it might be like SQL alchemy, providing an interface to all the different SQL databases or something like that

It's like a common interface to many backends.

**Sasha Dagayev:** Yeah. And what it also does is that it makes switching stuff out very simple. So if you want to switch out a language model , let's say you were using OpenAI, but now you want to start using Gemini. With Langchain, that's very simple to accomplish. You Don't have to rewrite the whole application.

**Phil Ewels:** Which is good news, especially when the models are being updated all the time, yeah fast moving space.

## AI concepts: Fine tuning

**Phil Ewels:** The other term was fine tuning. What is that and that fit in?

**Sasha Dagayev:** So fine tuning is essentially putting a set of matrices on top of the foundational set of matrices in order to dramatically change the result.

So what that means is that if you would like your model to never respond in 10 letters or more, and you've tried all you can with the prompt, but it's not happening. You can do something called a fine tune, which essentially takes the matrices that are in the model itself, and augments them with a much smaller matrix on top.

So the way to think about it is like putting an adapter on your electrical plug or something like that, to make it work specifically for a specific application.

The pros and cons Is that it can dramatically improve performance of your model. The con is that it's a little bit more expensive to serve because you're introducing some additional complexity. And that a lot of times you'd be much better off just changing the structure of how you prompt the model.

I think the use cases where fine tuning works really, really well is if you want to tailor your model for a very specific task. For example, code autocomplete, right? For those of you who have used system prompts to try to force a model to not say here's your code, but just give only the code. It can be very frustrating because you'll still fail, one percent of the time, which in engineering will be pretty jarring experience for your user.

**Phil Ewels:** Just imagine being in VS code and hitting tab and then getting this huge explanation about what the code snippet does.

**Sasha Dagayev:** Exactly. And so the way that a lot of these autocomplete applications work is that they fine tune pretty aggressively. So that way it works towards specific languages and it is set up for giving you just the right amount of autocomplete, right?

Because some of these models, they can return, eight pages worth of code. You don't want eight pages. You want, a little code snippet and things like that. So it's really good for very specific applications, but you have to have a tight niche for it. I'd say fine tuning a model just for all of bioinformatics is a bad use of your time.

I think if you wanted to fine tune a model for returning, something like DSL2- great application for it.

**Phil Ewels:** Sign me up

**Sasha Dagayev:** Yeah.

**Phil Ewels:** Cool so that was a really good introdution, thanks. I learned quite a lot there.

## Challenges in bioinformatics AI

**Phil Ewels:** Can you tell me a little bit about some of the specific challenges you found? Thinking with pipeline generation, when I've tried it with just vanilla chat GPT , some of the Nextflow code, it gives me, I get a lot of hallucinations. I get DSL1, I get, it's almost never runnable.

Did you hit anything like that? And what was your approach to it?

**Sasha Dagayev:** Yeah. Tons. You do a tons of it. The highest alpha, like most actionable items are all in the prompt. There's three axes that are really important for your specific application.

So one is defining a character for the LLM to play. So that way it can do that specific task. Number two is providing examples of what a successful generation looks like.

And then number three is be nice. If you're nice to the LLM, I don't know why, but like our future robot overlords produce better code if you say please, thank you and all that stuff.

**Phil Ewels:** That's amazing. It's yeah, doing the right dance before you step into the lab. If you say all the right things to the agent, you'll get a good response.

**Sasha Dagayev:** Yeah. So one of the things that you'll see, for example, for your use case with DSL1 / DSL2 is if you say, generate an RNA-seq pipe for me and you show it let's say ATAC-seq, pipe, but that's written in DSL2, the model knows about DSL2, right? It's just that in its foundational literature, it has more information on DSL1. So when you ask it, it's going to be like, okay, I have, my matrix weights are more towards this format. So that's what I'm going to give you.

**Phil Ewels:** And then by giving it a DSL2 example, it steers it back in the direction you want to go in

**Sasha Dagayev:** Exactly. One of the things that I've seen in examples online is that I see a lot more Flask web applications, which I didn't think was as popular as before chat GPT. But I think that because there were so many tutorials from the mid teens that they've trained on, now it's like, Oh, when you want a web app, here's Flask. So a lot of these older packages are really going to have a much longer shelf life because that's, what's in the models.

**Phil Ewels:** That's fascinating. And again, makes the importance of good documentation twofold. Not only is it important and critical for the users, but it's also critical to give material to the AI agents.

**Sasha Dagayev:** Exactly. Yeah. Yeah. And for everyone that's like making, open source software, I would say make sure to understand that your manuscript and your GitHub will be in GPT5, GPT6, in GPT7.

## AI concepts: Q&amp;A training

**Sasha Dagayev:** And the way that they generate the data is that they ask questions against it. So they will load in your manuscript or they will load in your github and they'll say, Hey, what is this software? What does it do? And then the model will generate training data for itself. Based on those question and answer pairs.

**Phil Ewels:** Is that not circular? it's like a kid marking his own kind of schoolwork and being like, Oh, I did great.

**Sasha Dagayev:** Yeah. So I thought so as well, but that's like how it's done.

OpenAI, it gets a lot of flack for being closed source, but they do publish some great cookbooks and one of their cookbooks on generation, they put out how we taught the models about 2020 Tokyo Olympics.

In there, what they did is they showed GPT a Wikipedia page about it and asked generate 10 questions about this. So generate the 10 questions. And then they said. Okay. Answer these 10 questions. And the questions were who won the a hundred meter dash, who won the 200 meter dash, what was the time and stuff like that.

And then those question answer pairs were used for the inputs for 3. 5, for four and so on and so forth.

**Phil Ewels:** So it's almost like a kind of data sanitation or standardization step.

**Sasha Dagayev:** Yeah.

**Phil Ewels:** Really interesting.

## Motivation for joining Seqera

**Phil Ewels:** We talked at the start a little bit about how you ended up bringing tinybio to Seqera. There are quite a lot of AI startups now, but there's also a lot of tech companies.

You talked about how Seqera is well known for Nextflow. Is there anything else that kind of particularly drew you over two companies together?

**Sasha Dagayev:** I think once I understood that once you sign up for Seqera it brings together the critical pieces for the context of a bioinformatician's workflow all into one place. Seqera is the perfect context manager for all of the things that a bioinformatician would interact with.

The platform takes your code from GitHub. And then it makes it very simple to then run that code on your AWS, your GCP, Azure, whatever cloud. And in doing that, it knows a lot about all of those things. So when things go wrong it can replicate a lot of what tinybio can do with its own home built environment. But on other people's stacks using information from other data sources.

Context for a specific job for a specific run is going to be the difference in those last, points of hallucination, cause you can have GPT, 45, if it doesn't know what your AWS IAM structure is, it's not going to be able to solve your permissions issue if it just doesn't know about it.

**Phil Ewels:** Exactly. You can see that when you go to the credentials page in Seqera Platform, you've got all the different accounts all listed side by side.

**Sasha Dagayev:** Exactly.

## What happens to tinybio now?

**Phil Ewels:** There's a blog post on the Seqera website, which I'll stick the link to in the show notes, which is all about tinybio moving to Seqera. How is that move happening? What happens to tinybio now?

**Sasha Dagayev:** So hopefully in the next couple of months we have feature parity with tinybio in terms of capabilities. And so once we have that finished, we will direct everybody from tinybio to Seqera which I think will be a better experience for almost all of the users and certainly much more functionality for, both user sets.

I would say that tinybio the domain will most likely still be, very functional and be able to service most of what you're doing today, at least through November, December. And then hopefully sooner than later, we have some level of feature parity on Seqera. So that we can offer you similar functionality there.

## Future plans

**Phil Ewels:** Obviously you can't talk about everything we've got planned in the future, but can you give us any little titbits of things that you're excited about going forward? Anything in the Seqera stack that you see as a shining jewel ready to be used?

**Sasha Dagayev:** Yeah. What I'm really excited about is I see Seqera Platform as a motorcycle for very experienced bioinformaticians. I think if you're a very experienced bioinformatician, you get so much out of it. However, I think there's this massive gulf between those people and the wet lab scientists.

What we want to do is make sure that the wet lab researcher can reach in to the Platform and be able to ask a question in English and get some information back about what happened, right? For example if you're working with a researcher and they know that an ATAC-seq run has been done, but they don't know which blacklist you used, for example. To find that information in the code is actually pretty challenging. To understand thatsome applications will use a pre filled blacklist, some will use a custom blacklist, answering specific questions about pipelines and runs that happen on the platform for a hardcore researcher, who's not familiar with the intricacies of the specific pipeline that you're running, I think will unlock so many people to start interacting with the platform.

**Phil Ewels:** That's really cool. I spent many years working in a core bioinformatics group where we're massively overwhelmed with the number of people coming to us with questions and, dealing with wait times and stuff. And I think anyone familiar with that situation will resonate with that.

We run it for you. Here's your data. And then the researcher can ask questions of the data themselves without having to come through us.

**Sasha Dagayev:** Exactly. Yeah. When it's clear that something's wrong, the questions that the researcher will have are going to be extremely time consuming for both parties because the researcher will know something very specific about mitochondria and the bioinformatician might know what mitochondria are but they don't, know the specifics of, the research.

And so it's okay, I have to put a meeting together. Like I got to teach you about this thing. You got to read the documentation. If we can facilitate that whole thing through AI, I think everyone will be significantly happier.

**Phil Ewels:** AI has got an unfair advantage that it can be an expert in everything.

**Sasha Dagayev:** But that's the point of providing the context, is that we don't think that AWS will ever provide you the context about the code that it ran, right? Specifically for bioinformatics. We think that's a great place to start for what we need to do.

**Phil Ewels:** Makes a lot of sense.

## AI and open-source

**Phil Ewels:** And last thing before I wrap up here, many people listening to this podcast are going to be heavily involved in our open source community. It's something that Seqera has grown from, something where we put a lot of time into. How does this interact with all of that? Is it a paid service for starters? That's an obvious question. And how do you see this interfacing with our community of Nextflow developers?

**Sasha Dagayev:** So for Nextflow developers, I think initially, the usage cap will be very generous. There's the realities of the fact that the foundational models that are made by, the closed providers, meaning Anthropic, Google and OAI, are still better than the open model.

So if you want better results, you have to use those models. And so we have to use those models as well. It's very likely that we'll be providing a service that has a much higher cap than, those direct providers to frontier models for this specific use case.

And I will say that we will definitely be the best place to generate pipeline code that is what you're looking for, right? So that way you don't have to craft that specific prompt. We already assume that, okay, if you're on Seqera, when you're talking about Picard, we know that you're not talking about Star Trek. We'll definitely be the best place to get the least hallucinated type of, information back.

**Phil Ewels:** And that's for code generation and doing analysis. Is it also going to help people who are trying to learn Nextflow and kind of interaction with those kinds of questions?

**Sasha Dagayev:** Our goal is to answer every question in an eighth grade level. So it should be very accessible. It's a 13 year old, right? so that one of the things that I've learned is that if you want the most people to understand what you're saying without losing the core of it, is that perfect amount. I find that even a lot of times that, we give explanations for what we work on in bioinformatics, we use words like &quot;elucidate&quot; or &quot;imputed&quot;, and we know that, for many of our colleagues, English is not their first language. So this can be communicated so much simpler. And reach so many more people.

So yeah, we will definitely have like the easiest way to QA your things because it's going to be very focused on bioinformatics and especially focused on Nextflow. So when you ask about channels on Seqera, we know you're not talking about TV, and it'll be a much more targeted generation of content.

**Phil Ewels:** Very cool

## Find more updates

**Phil Ewels:** So Anyone who's listening to this podcast and they're excited about this. Where can people find out more? Where should people be looking for updates and news?

**Sasha Dagayev:** So we're aiming to have our first release for the Summit in Barcelona. So tune in there. It's not done yet, but we're confident that we'll have something to show for Barcelona and then have a much bigger release for Nextflow summit in Boston, but we're hoping to be much more iterative. Especially given that, the team is in the AI space and everything's, constantly changing.

**Phil Ewels:** Fantastic.

And you're going to be at the Summit, presumably Sasha?

**Sasha Dagayev:** Yep. I'll be there. So yeah, Sasha GPT will be there.

**Phil Ewels:** Going to join the nf-core hackathon or anything and see if we can do some cool stuff?

**Sasha Dagayev:** Oh, yeah, absolutely.

## Wrap up

**Phil Ewels:** Fantastic. Thanks so much for joining me today. Absolute pleasure to talk to you always, and actually really informative for me. I learned quite a lot about AI and I'm really excited about some of the use cases we've been discussing here.

I'm really excited to do this again in a year. See where we get to. But yeah, thanks very much and thanks everyone for listening.

**Sasha Dagayev:** Thanks guys.

:::
