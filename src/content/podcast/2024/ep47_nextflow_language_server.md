---
title: New Nextflow language server & VS Code extension
episode: 47
description: Take your developer experience to the next level.
date: 2024-10-05
type: podcast
subtype: Technical discussion
youtubeid: ajGX3Tp_iBg
image: /img/podcast_ep47.png
tags: Nextflow
author: Developer advocates
icon: logo_podcast_channels.jpg
---

In this episode, [Phil Ewels](https://github.com/ewels/) and [Ben Sherman](https://github.com/bentsherman/)
discuss the launch of the new Nextflow language server, a significant upgrade providing advanced code intelligence features such as code completion and error hints for VS Code users.

They detail the benefits of formalizing Nextflow as its own programming language, reducing reliance on Groovy, and improving error messages and code clarity. They also cover new documentation, plans for future features like type annotations, and encourage community feedback and adoption.

<!-- end-archive-description -->

- Download the extension
  - [Microsoft Visual Studio Marketplace](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow)
  - [Open VSX Registry](https://open-vsx.org/extension/nextflow/nextflow)
- Announcements
  - [Paolo & Ben's talk at the Nextflow Summit](https://www.youtube.com/watch?v=yAMvfs9gw4o)
  - [Blog post with announcement](https://seqera.io/blog/modernizing-nextflow-developer-experience/)
- Getting help
  - [Community forum `vscode` tag](https://community.seqera.io/vscode)
  - [New VS Code extension docs](https://www.nextflow.io/docs/latest/vscode.html) _(includes troubleshooting guide)_
  - [New language syntax reference](https://www.nextflow.io/docs/latest/reference/syntax.html)
- Source code
  - [VS Code extension GitHub repo](https://github.com/nextflow-io/vscode-language-nextflow/)
  - [Language server GitHub repo](https://github.com/nextflow-io/language-server)

**Key features include:**

- **Enhanced Code Intelligence:**
  - With code completion, hover hints, snippets and syntax highlighting, streamline your workflow coding like never before.
- **Code Navigation & Hygiene:**
  - Effortless navigation using control-click and Peek features allow you to fly around your code base at speed.
- **Code Formatting:**
  - Solve your indentation woes automatically on save, with auto code formatting.
- **Errors and Warnings:**
  - See problems with your code inline as you type, with immediate feedback
- **Future-Proof Language Evolution:**
  - Introduction of Nextflow as its own programming language, promising a more strict and streamlined syntax with less Groovy magic.

## Podcast overview

## Introduction

Welcome to episode 47 of The Nextflow Podcast, where Phil Ewels, Product Manager for Open Source at Seqera, and Ben Sherman, Senior Software Engineer at Seqera, take us through a significant milestone in the development of the Nextflow language. This episode is packed with exciting news from the recent Nextflow Summit held in Barcelona and introduces groundbreaking advancements, including a newly enhanced Nextflow language server. This update promises to transform the way developers interact with Nextflow code, simplifying complexities and improving the overall development experience.

## Evolution of The Nextflow Podcast

Before diving into the latest updates, a quick note for our keen listeners: the podcast has been rebranded as "The Nextflow Podcast," dropping the "Channels" from its name. This change reflects our intent to streamline the podcast's identity, making it easier to find on platforms like Spotify and simplifying access via seqera.io—another step in Seqera's strategic efforts to unify its digital presence.

## Nextflow Language Server – A Game-Changer

The highlight of the episode is the announcement of the new Nextflow language server, which brings sophisticated code intelligence features to the Nextflow VS Code extension. This innovation upgrades the existing extension that currently provides basic syntax highlighting, by adding capabilities like code completion, hover hints, and inline error detection.

Phil and Ben provide insights into the development journey of this feature, a year-long endeavor driven by feedback from the Nextflow community. This new language server addresses a longstanding request for improved error messages and developer experience. It represents a leap forward by considering Nextflow as its own first-class programming language, independent of its Groovy roots.

## Beyond Syntax: A Vision for Nextflow as a Standalone Language

Ben discusses the shift of Nextflow from being just a Groovy DSL, evolving into a standalone language. This allows for a more specific and precise syntax, facilitating better code intelligence and clearer documentation. This requires balancing the freedom traditionally associated with Groovy while introducing structured guidelines for Nextflow scripts.

## Code Interactions with VS Code and The Language Server

The new language server is designed to offer an interactive experience, running in the background as you edit your Nextflow scripts in VS Code. It rapidly parses code, providing immediate error feedback, including variable references and argument count mismatches, greatly enhancing the development workflow.

## Community Involvement and Future Directions

The podcast encourages community involvement, inviting developers to extend support for additional text editors like IntelliJ, NeoVim, and Emacs. Phil and Ben stress the importance of community feedback to refine these new tools, with robust support resources available through Seqera’s community forum and GitHub repositories.

## Upcoming Innovations

Looking ahead, exciting developments are on the horizon, including potential type-checking capabilities, which could provide more robust error detection and lead to more reliable pipeline executions. The conversation closes with an enthusiastic outlook on the language server’s future role in Nextflow's evolution.

## Conclusion

This episode underscores a pivotal moment in Nextflow’s evolution, driven by the demands of its vibrant community and the ingenuity of its developers. With enhanced features coming to the VS Code extension and broader language improvements, the future promises to simplify and empower the Nextflow development experience. As Phil and Ben wrap up this insightful episode, they express eager anticipation for user feedback and collaboration in shaping the future of Nextflow. Stay tuned for further developments, and happy coding!

# Full transcript

:::note{title=Transcript}

## Welcome and introductions

**Phil Ewels:** Hello, and welcome to The Nextflow Podcast, supported by Seqera. You're joining us on episode 47, which is going out in November, just after the Nextflow Summit, being held in Barcelona.

My name is Phil Ewels. I'm Product Manager for Open Source at Seqera. And today I'm joined by Ben Sherman, Senior Software engineer at Seqera.

Ben, thanks very much for joining!

**Ben Sherman:** Hey, Phil, thanks for having me on.

**Phil Ewels:** It's a pretty exciting episode today. I've got to say, we've been building on this for quite a while and we've got some really big new stuff to show to the world. I hope everyone else is excited as I am about this.

## Podcast rebranding

**Phil Ewels:** Before I get into it. Any of you listening with keen ears might have noticed I introduced the podcast slightly differently to normal.

We've dropped the &quot;Channels&quot; part of the podcast name. So we're simplifying things a little bit and now just calling it &quot;The Nextflow Podcast&quot;. Should make it a little bit easier to find if you're searching for it on Spotify or anywhere else. A little bit clearer. And the other thing is anyone looking for the podcast webpages, they're going to be moved over to seqera.io.

Because of the way that our company grew and all the different projects started as their own little stand alone things, we have many websites scattered all over the place. We're trying to consolidate everything a little bit under one roof and have a global search, which will work for everything. So this is just one step. Everything will be redirected, so you shouldn't really notice much, just a shiny, nice new design.

## Nextflow Summit announcements

**Phil Ewels:** Anyway, let's kick off. Got to confess, we are actually recording this a little bit early. It's a hectic few weeks for us. Firstly, Ben, congratulations on your demo. It was fantastic, I'm sure. You knocked it out of the park.

**Ben Sherman:** Thank you. It was, and if it wasn't, I had nothing to do with it. Future Ben's fault.

**Phil Ewels:** Future Ben, exactly. Yeah, no, tell us about it. What are we going to launch at the Summit then? Give us a quick intro.

**Ben Sherman:** So, if you're watching this now, we will have announced a new Nextflow language server, which comes with a sort of formalized definition of the Nextflow language, as well as the configuration files. And we will have given a demo, of using the language server in the VS Code extension. We already have an excellent VS Code extension that provides basic syntax highlighting, but now we've supercharged it with this new language server that's going to give it all the sorts of code intelligence that you expect from a programming language like Python, things like Code completion and hover hints and inline errors, all that good stuff.

**Phil Ewels:** Just, I can't emphasize how happy this makes me. I'm also a bit jealous because I think past me, would have absolutely loved this for all the years been writing Nextflow code. This has got to be one of the most requested features for Nextflow, right? Every time we run a community survey , we hear error messages and developer experience.

## Language server background

**Phil Ewels:** So Ben, we've been working on this for how long now, like a year or so, tell us a little bit about where we came into it and the motivation and some of the past work.

**Ben Sherman:** I guess I've been working on it full steam ahead for about six months. Really been thinking about it for a couple of years at this point. I was a user before I was a developer. I experienced all these same issues with writing code and debugging.

After I joined Seqera, I created a GitHub discussion. I think it's just called Nextflow language improvements. It was basically a catch all where I said, Hey, everybody, just bring all of your grievances all in one place so we can triage them as best we can. We got a ton of great feedback, and we actually were able to fix a lot of things on the margins. You know, We clarified some documentation, fixed a couple of errors, trying to improve certain errors where we could, but we hit a certain limitation with many errors where in order to improve them, we would have had to go a lot deeper, rip out entire parts of the Nextflow runtime and rebuild them ourselves.

At the time I thought that would be completely untenable, because Groovy is such a complicated language in terms of its implementation. There's a lot of things wired together in it. And I wasn't sure that we would be able to really support that well.

**Phil Ewels:** Like you say, because Nextflow is this Groovy DSL goes through this compilation step, everyone I know who's looked into this problem has come up with one of two things.

Either, they've come up with a tool which uses regular expressions and just looks at the Nextflow code as if it's just one big string and tries to do stuff with it. That's how the existing VS Code plugin works, right? The syntax highlighting works with TextMate grammar with regular expressions. That's how the nf-core code linting works, a combination of that plus the Nextflow config output.

And then the other approach is we had, Eric Danielson, who was an intern at SciLifeLab working on nf-core, and I was like, right, he seems to know about this stuff. I'm going to throw him at this problem. And he spent a good few weeks looking at it, and eventually wrote up, basically, this is impossible unless we build a language server. And I think many people came to the same conclusion.

**Ben Sherman:** I actually just listed to the podcast we did a little over a year ago on the same topic where we talked about a lot of these different things. And it was so strange listening to it, considering where I am now, I feel like so much has happened since then. Not that I've aged 10 years, but like in terms of the stuff that I've learned, I feel like I've lived half a lifetime.

At the time I was not very confident that we could really do it well, the best I could think of was, let's just take the groovy compiler, stick it into a language server, and at least then you could get some amount of feedback in the editor. Instead of having to run the script, the errors would actually be put into the editor. Then you could iterate a little bit more quickly.

But that still wasn't solving the problem. The groovy compiler does not understand the Nextflow language. The Nextflow language consists of all these different little parts, sort of duct taped together. That has made it difficult to provide really good code intelligence, things like, understanding the structure of a process or a workflow and knowing, what are the available functions and built in variables whenever I'm writing code at any given point.

**Phil Ewels:** There's the integration with the editor, but there's also the specificity of where the error message is, right? You have an error message, which says something like compilation error workflow, and you're like, well, that workflow might be a thousand lines of code, and there's an error in there somewhere,

**Ben Sherman:** Right.

**Phil Ewels:** because it doesn't know where the error is actually coming from.

**Ben Sherman:** I actually brought this up with some of the Groovy developers on the Groovy Slack. Saying, Hey, we get these weird error messages sometimes with our little DSL and they basically said, yeah, the syntax is really complicated and sometimes the parser isn't able to pinpoint exactly where the error happened.

So even they were saying, yeah, it's a really difficult problem.

## Nextflow as a programming language

**Ben Sherman:** And the solution that I found was, we have to implement our own parser for Nextflow. So bypassing the Groovy parser, and we implement all the syntax that we want to be, supported in the Nextflow language.

When you do that, you're able to cut out a lot of syntax in Groovy that we don't need or typically use. That allows the parser to be a lot more specific in a lot of these cases. There's still some difficult edges here and there, but again , even in those cases, because we have the language server and the editor feedback, you can still pull yourself out a lot more easily than when you were having to run Nextflow just to get the error feedback.

**Phil Ewels:** Does this mean that Nextflow is not going to be a Groovy DSL anymore? Are you saying we're ditching groovy? How does this work?

**Ben Sherman:** That's ultimately where this is leading. Implementing our own parser led me into this philosophical question of, should Nextflow just be a superset of Groovy or should we actually, treat Nextflow as its own programming language?

Now there's a lot of, implications of that. If you want to have your own first class programming language, you have to be able to provide really good tooling and really good developer experience. But if you can do that, that actually provides a lot of clarity.

In the documentation now, we can say, here is the full specification of what is allowed in the Nextflow syntax. You can do this, you can't do that, so on and so forth, rather than just saying, whatever you can do in groovy, you can do in Nextflow.

This has also required us to think through, what syntax do we want to support? Because, just about any kind of Groovy syntax, you can find some example out in the wild of somebody using that in, in their Nextflow code, because it's there and maybe it works in their situation.

But what we're trying to do now is restrict it down a bit to the syntax that is actually intended to be used in Nextflow scripts. And in cases where maybe somebody was using some esoteric Groovy syntax that we're not supporting anymore, we can evaluate, okay, what were they trying to do at a higher level? And how can we support that?

And as long as you have that clarity of, this is what's allowed, this is what isn't allowed. And lots of examples of here's how to do , pattern A, pattern B, pattern C. I think that will be a huge benefit overall. All the groovy wizards out there, they'll find a way to trim down their code where they need to. The vast majority of cases, whatever code you are already writing will just work as before.

**Phil Ewels:** There's definitely a pattern I feel like I've seen. You start writing Nextflow code, and then you get to a certain point , and you start digging around in a groovy documentation and you basically end up learning two programming languages, instead of just one.

Also a bit more Pythonic, dare I say it, one way to do everything.

**Ben Sherman:** Yes, we're certainly trying to trim down the many ways to skin a cat kind of problem, where people are writing things three to five different ways. By the time you guys are listening to this, the documentation will have been updated with, some new pages on things like syntax description.

And in the future, we intend to go even further documenting things like the standard library. The most common things you might use: lists, strings, maps, instead of just giving a few code examples, actually going through here, all of the methods that are available for those kinds of data structures. Going a lot deeper into the actual libraries that are available. While still being able to access the full Java, Groovy standard libraries.

So if there are classes that you use, for HTTP requests or JSON parsing, YAML parsing. You'll still be able to use. It's just the syntax, will be a little bit more strict.

## Building a language server

**Phil Ewels:** How does this tie into the language server and then later VS Code.

**Ben Sherman:** So first of all, AST stands for abstract syntax tree. You can think of it like a machine readable data structure representation of your code. So when you write code, you're just writing text. And then when Nextflow wants to execute that code, it first has to translate that text into some sort of data structure that it can manipulate.

And so the approach that basically every programming language uses is the abstract syntax tree, where you represent your code as some sort of tree structure.

With Nextful, you can think of that as something where, you know, the top node of the tree is a script. And then a script might contain any number of declarations. Processes, workflows, functions, include statements. And so that would be one node with a bunch of children. Then with each one of those, you might have further children.

Of course, you can imagine these things can get very complicated very quickly. But the nice thing is that because it's a tree, the child of a tree can also be a tree. You can have this infinitely hierarchical structure.

And so that's how we represent code. Now, the Groovy AST already exists. So what I've been able to do is to basically extend the Groovy AST with custom nodes, like process node, workflow node, include node, script node.

Then represent everything in those terms, but reusing as much of the groovy code as we can. Because this is all a huge amount of work to build from scratch. So one of the reasons why we're able to do this at all, why I've been able to do it on my own is because we're reusing a ton of the Groovy parser logic, and even a ton of, language server stuff, it already exists and we're able to build on that.

So you also mentioned language servers. Basically you think of Nextflow itself is a compiler and interpreter, it takes code, parses it into an AST and then executes it. And if it runs into errors, it just spits out the errors and then quits. The language server is designed to be more of an interactive experience.

When you spin up your editor and you have the language server installed, it will run as a background process. Sort of like a &quot;Nextflow server mode&quot;, except it doesn't execute your code. It only knows how to parse code and report errors and do certain kinds of queries.

So when you open a Nextflow script, right out the gate, it will parse that script and it will actually scan your entire workspace and find all the Nextflow files, all the config files. Parse all of them and then maintain an internal model of those. And then as you're editing code, it will continually update the parts that need to be updated.

When you do things like hover over a symbol or try to navigate through your code. The language server will act like a server would. The client will send requests like, Hey, send me a list of completions that I can do right here, or send me a hover hint for this symbol the user just hovered over it and the server will send that back. And it's essentially just querying all the ASTs that it's gathered from your workspace.

**Phil Ewels:** And where is this code? Is this part of Nextflow itself? Is it open source? If people are into this kind of thing, where can they go?

**Ben Sherman:** So at the moment this new parser has been implemented entirely in the language server. And the language server is a Java library, similar to Nextflow. So it's a separate project.

At this moment in time, while I'm talking, it's the repo is still private, but when you guys are listening, we will have published the repo. So you'll be able to see all the code.

Long term, what we want to do is move all of this parsing code into Nextflow itself. So that whether you're using the language server or running a Nextflow or script on the command line. You'll get that same level of parsing and code intelligence.

Then from then on when we want to implement new language features, we have to implement that in Nextflow itself.

**Phil Ewels:** While they're separate, it means you can get different error messages between the language server and Nextflow. You can get one error message in VS Code and then you run your pipeline and you get a different error message.

**Ben Sherman:** This is something probably worth warning people about. For the most part, the language server will essentially be a stricter version of Nextflow. It should give you all the errors that Nextflow would give you. And maybe a few more to boot. And then once we incorporate all that parsing code into the Nextflow itself, then all of that will be unified and you'll get the exact same experience with both.

## VS Code interactions with the language server

**Phil Ewels:** How does VS Code know when to use the Nextflow language server?

**Ben Sherman:** It's basically file extensions. When you open a file that's .nf or .config that will trigger the Nextflow VS Code extension. It will spin up that Nextflow language server if it's not already running.

So it's a very thin client, basically just start the server, connect it to VS Code, and then VS Code and language server will handle the rest. The VS Code extension itself is actually, pretty simple. And this is the intention of the language server protocol is that you write, all of the meat of the code intelligence and parsing in a language server, and then, you can hook it up to whatever editor you're using by writing a very thin client that just spins up that server. So for now we just have VS Code support, but in the future it would not be very difficult at all to add that same support to IntelliJ, NeoVim, Emacs, and so on.

**Phil Ewels:** Brilliant. So this is where we appeal to the community then, people can go out and write their own IDE plugins that use that language server.

**Ben Sherman:** I know there's at least a few people using NeoVim and Emacs, so I won't be surprised if they figure it out. Yeah, I hope over time, whether by us or through the community, we'll be able to extend that support to all the major editors.

## VS Code demo

**Phil Ewels:** Should we have a look at it?

**Ben Sherman:** I suppose we can!

**Phil Ewels:** Let's see what all the fuss is about. So what have we got here?

**Ben Sherman:** So I've got open, the nf-core/rnaseq pipeline. You can also see we're right in the thick of it because this is literally my local development build that we're running. So let's, pull up a Nextflow scripts and see what happens.

So first of all, you'll see, down here, it's initializing the workspace. Depending on the size of your project, it could take, one or two seconds or, 10 seconds.

## Errors + warnings tray

**Ben Sherman:** Once that's done, very likely you will see some errors. How many errors you see will depend on how idiomatically you write your code. RNA Seq, we've actually trimmed down a lot of the errors in RNA Seq, but there's still some, but that's fine because I actually wanted to show some of the common types of errors you might run into.

Going from the informal DSL 2 that we do, to the new Nextflow language specification.

**Phil Ewels:** So you just brought up a kind of a tray along the bottom of VS Code saying problems. I think in my VS code it shows the number of errors and warnings, bottom left as well, and I can just click that, and it shows the same panel.

**Ben Sherman:** So if you go in here, you'll see we've got about a hundred or so errors and warnings across all our Nextflow scripts and conflict files.

You can search for different kinds of things. We can do something like this and it would pull up, all of these annoying, unexpected input errors. .

Another thing I want to point out is if you go to output. And then you go to Nextflow language server. This is the language server logs, the internal logs. You might not see any by default because most of them are debug logs that are disabled by default. If you ever run into some wonky behavior and you're not sure what's going on, you can always pull this up and see maybe you've got an error message. That's very helpful to me when you're reporting issues.

## Errors: Level 1

**Ben Sherman:** So this is the top level main script for RNA Seq. We've got a couple errors here. Let's try to fix some of these errors so that you see how to deal with them.

Hopefully the error messages you get will be helpful on their own, but some of them are difficult to handle. So in this case, we've got these includes and it's saying module could not be parsed. What this means is that the file exists, but there was some kind of syntax error in that file. And so it wasn't able to parse the file, which means it wasn't able to pick up these definitions.

I point this out because if you had something like an invalid file here, I just added an S here. It would say invalid include source. So that's how you tell whether it's an issue with the file not existing or an issue with the code itself being wrong.

Now, a nice thing about this is that all of these include sources are links. You can see they're underlined. I'm going to control click on this to navigate to it. You could also hit this, follow link to, to get, go to it that way.

So now we've opened up this other file, the easiest way to look at this, you can either see on my sidebar here, there's some red, so we can click on that to go on here to see what the error is.

You can also go to it this way. So basically, filter out the file you're on, and then you can click through the errors here. Either way, however you do it. We can see now that there's a syntax error here.

This is like the level one error, where the parser hits some sort of character. And it doesn't know what to do with it. And it's so bad that it just quits entirely. So it's not able to just skip over it and keep going. It's basically said unexpected input. I'm giving up.

So this will happen if you have some kind of syntax that isn't recognized. In this case, it looks like somebody thought they were in Python and tried to do some kind of named argument thing, which is not how it works in Groovy, but this happens to be valid Groovy syntax.

This is basically an assignment. So it's going to assign that variable to itself. So it's just a no op. This is one of those examples of things where people might've done something that was invalid from the Nextflow point of view, but because Groovy has so much syntax sugar, it happens to pass the Groovy parser. And so you don't even notice.

This error is harmless, it doesn't really do anything. But you can imagine other things where maybe it's a valid syntax, but it causes some weird runtime issue down the line.

**Phil Ewels:** It's underlined the exact line, the exact character, where the parser hit trouble. Because even though the error message here is not very informative. The fact that it's on the equals symbol the exact character means, you don't really need the error message. You can guess what's going on here

**Ben Sherman:** Right. You can play with the code and get immediate feedback. If you encounter one of these errors, right after you type something, you probably know that the thing you typed is related to that error.

So let's just take those out. And once you do that, you'll see a few other errors pop up here.

## Errors: Level 2

**Ben Sherman:** Now Nextflow has been able to parse the entire file and it's giving us some more errors.

First of all, you've got this import statement. Now, this is something from Groovy, that we've decided not to support in the formal Nextflow language, because we already have these include statements and it'd be confusing to have both. But this is typically how people will bring in classes from the Groovy or Java standard library.

We may add some kind of syntax for this in the future, where you include groovy JSON, JSON slurper instead of import. And at least in the language is consistent, but for now you're going to get an error for this.

Actually, we're in luck. It looks like they're not even using it in this file.

**Phil Ewels:** Orphaned at some point

**Ben Sherman:** Wherever you were doing something like new JSON slurper, you can just remove the import and then use this fully qualified name, inline.

And then let's move back down to these other errors. So it says, statements cannot be mixed with script declarations . What this is saying is that- The distinction between top level declarations, like processes, workflows, functions, includes , feature flags. And then there's statements, which are things like this assignments, function calls, if else statements, try catch, process calls, operator calls. Statements can only exist inside a declaration. in the workflow, you would put statements.

This is a carry over from DSL1, because it used to be that you would have your process definitions and then all of your workflow logic would just be scattered spaghetti style all throughout your script. And then in DSL2 we introduced workflows and as a best practice, you should move all of that extra code into your workflow so that it's well packaged, but these sort of top level statements were still allowed as a convenience.

Nextflow language server is going to put its foot down and say, you can't mix stuff like this anymore. And so what we can do instead is, just move this into a workflow.

You can do stuff where you have code like this, just out in the wild. And as long as you don't have things like processes and workflows, you can This is still valid. A code snippet.

**Phil Ewels:** You've just written println &quot;hello!&quot; there. There's no process. There's no workflow. There's no scoping.

So this is the kind of thing that you might see in the Nextflow docs, for example, little 2-3 line snippets and examples.

**Ben Sherman:** Yeah. It's mainly used for documentation training. If I try to use a script declaration, like a workflow, now it's going to say, okay, now you've got a workflow in here. So you need to do something with this statement. Anyway, that's just an aside.

## Errors: Level 2 (3?)

**Ben Sherman:** Let's go back in here. So now the sidebar is really lighting up. So this is like. Level two error checking, you've entered the next level, which is name checking. So now it's going to do things like check all of your process workflow definitions and all your variable declarations and make sure that everything is linked up properly.

Here, this warning is saying this variable was declared, but not used. It's just not an error like this, there's nothing wrong with doing that, but it's just giving you a warning to let you know, and you can see that it's because it's, it wasn't referenced anywhere in this workflow.

## To def or not to def

**Ben Sherman:** So there's some errors here around the variable being declared without def. This is just another one of those little nitpicks that the language server is a little more strict about.

**Phil Ewels:** This one's come up before, the def thing, and it gets asked quite a lot on Slack. Should you use def? Should you not? What's the differences? I can never really remember the definitive answer.

**Ben Sherman:** This is actually a really important thing to bring up because this can actually lead to a race condition. So what's happening here is that we're in a function. I declare this variable. Because I declared it without def, in groovy land, that's essentially treated as a global variable. This looks like it's a helper function. That's maybe only called once in the workflow, but you could imagine, Nextflow is a parallel concurrent system, right? So imagine if this function was being called many different places all at the same time, they would actually start clobbering each other because they would all be writing to the same variable.

And so it's in functions in particular where the language server will make a fuss. In things like workflows and processes, this is not as big of a deal. So for instance, this summary params is being declared without def. That's fine because Nextflow actually is able to scope it to just that workflow. So you won't hit a race condition here.

The only case where, declaring a variable without def has a use, is in a process , So I'm going to pull up MultiQC, let's say that I want to declare a variable down in the script body, and then use it in an output. This is a very common thing. In order for the output block to see that variable. I actually have to declare that variable without def. This is basically the one case where it's useful to declare without def. Everywhere else I would use def just to be safe.

**Phil Ewels:** And that's basically because you're trying to get outside the scope of the script in this case.

## Variable references

**Ben Sherman:** Anyway, some other things I wanted to talk about. So we mentioned variable references. So here we've got this status variable. Let's say I make a typo . It's going to say, first of all, you declared this variable, but you didn't use it. Second of all, you've got this variable that doesn't exist, right? So you get that instant feedback. This is probably one of the quickest things language server will help with is just catching typos like this.

**Phil Ewels:** Probably one of the most important as well. I love this.

**Ben Sherman:** Yeah, because you'd like to be able to catch these errors in the editor and not when you're halfway through a pipeline run.

**Phil Ewels:** Exactly, it shortens that feedback loop, otherwise I would save this, do a load more work, then spend five to ten minutes running a test workflow, and then find some kind of error, and then have to trace back, whereas just being able to see it light up like that straight away , it's so easy to make those typos.

## Number of arguments

**Ben Sherman:** Another thing in this level two checking is when you call a workflow or process, it will compare the number of arguments you gave to the number of inputs that are actually defined. And so we can see it says and this is probably the case where you really want it. It's like you gave me 24, but actually there needs to be 25. And you're probably not going to be able to tell just from looking at it, whether there's 24 or 25.

**Phil Ewels:** Kind of an extreme example, but.

## Hover hints on workflows and processes

**Ben Sherman:** Yeah. And then also when you hover over the workflow process, you get a little hover hint here. So it's actually showing those inputs and what they're called. That can be. Nice to have side by side.

**Phil Ewels:** That hover hint isn't showing the entire workflow script. It's actually showing a cut down summary of that workflow or that process.

**Ben Sherman:** So it's looking at the definition and it's just taking out the inputs and outputs.

## Code navigation

**Ben Sherman:** I feel like you and I both were saying this, in the process of developing this language server, we've been discovering more about VS Code in general. So VS code has all these weird features that I didn't realize exist.

One of them being, I can control click on a workflow or process, and it will actually take me to the definition.

**Phil Ewels:** And yeah, this works in Python as well. . I'm like, ah, these cool features I never knew existed.

**Ben Sherman:** Yeah, this has been just as much an exercise in learning how to use VS Code as it has been learning how to implement language features for Nextflow. So it's been a lot of fun.

## Output hints

**Ben Sherman:** One more thing here is the outputs. We've got this dot out syntax for accessing named outputs. If you have a typo it will let you know. And again this is where the hover hint comes in handy because we can go down here to the emits and see, okay, was it version or versions? Okay. It's versions.

**Phil Ewels:** That's the kind of thing where you could remove an output from a process and forget that it's being referenced somewhere else in the workflow. But now you're going to see that error pop up,

## Code peeks

**Phil Ewels:** Can you show the code peak? Cause that was another VS Code feature that I just didn't know existed.

**Ben Sherman:** Yeah. So let's go in here. So right click on the symbol again, there's go to definition, which is the same thing as the control click. And then there's also this Peek. So if you don't want to pull up a whole nother file we can also just do this peak and it'll pull up a sort of an inline view of the definition.

And then also let's go over here and go to Peek References. It's showing you all the places where it's used. You've included it here. You called it here. And then now you're also referencing all these outputs. So you can see literally every time that symbol is being used.

## Getting help

**Phil Ewels:** Just before we forget, if people try this out, where do they go if they find some weird behavior? How should they tell us?

**Ben Sherman:** I think the main place we should go is the community forum. So community.seqera.io there is now a VS Code tag, where you can post questions, bug reports, whatever you want. And then we'll get those sent to the right people. Mainly me. There's also GitHub repos. There's a GitHub repo for the VS Code plugin, and then another repo for the language server.

**Phil Ewels:** There's a couple of things I want to touch on. There's three I want to talk about.

## Code formatting

**Phil Ewels:** One, I, I've said this before, but it's burnt into my memory, of spending half a day formatting someone's DSL1 script with proper code indentation. A little bird tells me that my days of hitting tab 10, 000 times are over. Can you talk a little bit about code formatting?

**Ben Sherman:** Yeah. Once you fix all your syntax errors and you fix all your name errors, um, that opens up a couple of features to you. One of them is formatting.

We can go to the command palette and either do format document. Or on Windows, you can do Shift Alt F, so I'm going to do that and then it will format your code.

Now, the formatter is not perfect. So we formatted it. You can see what it's doing, it puts everything in a certain order. You'll have your includes and then your params and then your entry workflow, and then your named sub workflows, processes, functions, and then everything is formatted as you would expect.

There's also a neat little option for, Harshil alignment, which is highly popular in nf-core. It'll try to align things like your includes, you can see the includes are aligned here. It's not a hundred percent. I'd say the formatting is like 98 percent there.

There still might be some weird cases that don't quite work. This is something that I think is still open to discussion. We can talk about what kind of white space conventions that people like and people disagree on. Maybe we can add some options or maybe we can just say, Nope, you got to do it this way. I think the language server currently will give a good starting point for us to move towards a more unified formatting convention.

**Phil Ewels:** And the plan is to bring that into a command line tool as well, right? So you don't have to be in VS Code to run an auto formatter

**Ben Sherman:** Yeah, it would be nice to have something like a &quot;nextflow format&quot; command where, every time you cut a release on your Nextflow pipeline, you just run format. And you don't have to worry about it. It just happens automatically.

**Phil Ewels:** I'm gonna be there on every commit

## Nextflow schema params

**Phil Ewels:** Nice and then a couple more things. In nf-core, we developed the Nextflow schema to describe and define parameters in a workflow. Can you tell me a little bit about that?

**Ben Sherman:** Yeah. We would actually like to incorporate this into Nextflow itself so that you don't have to have the nf-schema plugin. As long as this file is here, Nextflow can pick it up and use it to validate the parameters natively. You don't have to do it yourself.

So we should be able to see that here in the language server because this file is defined alongside main script. In the entry workflow, the synonymous workflow. That's where this JSON schema is visible. And so it will load that schema and use it to provide all sorts of validation.

So I can I can hover over this param, and it's given me the docs from the JSON schema. Similarly, if I type it wrong, it'll say, Hey that parameter is not defined in the schema. And it should give you code completion. So you can see the different parameters that are available in the documentation.

**Phil Ewels:** Very cool, especially on bigger workflows.

**Ben Sherman:** Yeah. Now, one thing I'll note is that this sort of validation is only available in the entry workflow, because that's where parameters are meant to be used. I know it's common for a lot of people to use parameters in a sub workflow. This is also a a carryover from DSL1, where everything was just one big script and parameters could be used wherever you wanted.

But again, we're trying to, restrict it down a bit. This is not something we're going to force on you. You can still use params wherever you want. What we're going to encourage is a best practice, try to move those params into the entry workflow. And if you need to pass it into a sub workflow, just pass it in as an input.

## Future warnings

**Phil Ewels:** It seems like a good time to mention that scary checkbox in the plugin settings.

**Ben Sherman:** Oh boy. Okay. Yeah, we, I guess we can show them this. So,

**Phil Ewels:** Let's go for it.

**Ben Sherman:** uh, future warnings. Paranoid slash deprecation warnings. They're turned off by default, but basically, if you want to get more hints about how the language is going to change in the future, you can turn that off and then it's going to recompile your project.

And then we should see a whole lot more yellow. So for instance, I mentioned this thing, it's variable being declared without def. Now it's going to tell you about that. And then also using params. Anywhere outside of the entry workflow , it's going to tell you about that saying, Hey, this might not work in the future.

I don't want to show this to scare people. Whatever transitions or breaking changes we do, we want to make them smooth. So we're not going to remove support for these kinds of things anytime soon. But as we evolve the language, if you want to be ahead of the curve , you can turn on this and go ahead and start trying to adapt your code to get rid of these.

But if you don't want to, you don't have to. We'll have plenty of time to introduce these changes. And of course that will come with documentation and training and clear communication on what we're trying to do.

**Phil Ewels:** actually think this is a really powerful feature. We're a big community with thousands of users and developers and an increasingly large team working on Nextflow itself. There's so much activity, it's just a huge volume of information to consume and you've got to be able to pick out the things which are relevant for your workflow. Syntax X or Y is going to stop working in six months time.

So I think this is a great way to pull that signal, out of that noise. If you have that magic tick box unchecked, you're going to update your VS code extension and suddenly you'll see a few new warnings saying, you've got a while to do this, but get ahead of the curve now and you won't have a nasty shock when you update Nextflow in six months.

**Ben Sherman:** Yeah. We hope it'll be a good way for us to communicate upcoming changes, without you having to come to us.

## Preview DAG

**Phil Ewels:** Yeah, there's one more sparkly magic jewel in the crown, which I thought would be nice to finish on actually. Just above workflow it says preview DAG.

**Ben Sherman:** Yeah. Let's just, uh, let's click on it. See what it does. So yeah, this is showing a mermaid diagram, a preview of this workflow, very similar to the DAG that Nextflow generates, same sort of mermaid syntax.

But what you'll note is that first of all, you don't have to run the pipeline to get it.

Second of all, it's only showing this workflow. It's calling some sub workflows here, but it's not expanding them out. It's just showing them. Then it's a lot more compact. So this channel that empty would normally show up as a node on this graph, but it's cut out a lot of stuff. It's only showing process calls, workflow calls, inputs and outputs. And so we can just go and click on any of these and pull them up and look at them all.

## Fixing new syntax errors

**Phil Ewels:** Is there anywhere that folks can go to get help clearing up errors? Some of them, it seemed like there was a pattern to them, standard things that might need fixing all over the place .

**Ben Sherman:** We've tried to collect all the most common issues that I've seen. Common patterns that emerge because they are allowed in Groovy, but we want to be a little bit more strict in the Nextflow language. And so there's going to be a documentation page somewhere in the docs, and we'll probably have a blog post as well. And there will be a page where we go through a lot of those common patterns, you might have been doing this before, but it's no longer allowed so do this instead, and it'll work. And especially with the nf-core pipelines, going through that guide should cover just about all of the issues you run into. We've also been updating the nf-core pipeline template. So , the basic template boilerplate you should see has fewer errors now.

## Config files

**Ben Sherman:** I guess there's one last thing we should show really quickly, which is the config files. The config syntax in particular is a lot more strict than before, because it used to be, you could just have arbitrary groovy code.

Now it's a little bit more of a, just a strict declarative syntax. I won't get into all that right now, just to say, it's a lot of the same code intelligence. So things like you have documentation for all the different config options.

**Phil Ewels:** We didn't mention this before, but those tooltips also, they say read more underneath when you hover over it there. That takes you to the Nextflow docs, is that right?

**Ben Sherman:** Yeah. I won't pull it up now because you won't see it anyway, but this is also a nice way to first get a quick blurb on what a thing is. And then if you need more details, the link is right there and you don't have to be constantly searching the docs for things.

**Phil Ewels:** That's something I've really got to enjoy the last kind of month or so of trying this out. I find myself clicking those links to the docs a

**Ben Sherman:** those little time savings are going to really add up. And I also think in the config, this is where the code completion really comes in. So we didn't really show it off in the scripts, but, there's hundreds of config options available and remembering how to type it out every which one.

Instead, maybe we can just, we can just go through the different scopes, and then now we can see, all the different options that are available. It's pulling up a description of each one, as you type, you can quickly see what options are available. This was actually one of the original motivations for me for the language server.

## Updating the plugin

**Phil Ewels:** So this is all just ready to go today, right? People can just update their VS Code extension for Nextflow. And maybe if you're not already using VS Code, you might feel motivated to give it a go.

It's on the Microsoft, marketplace . It's also in the Open VSX marketplace for open source distributions of VS code.

**Ben Sherman:** Right. Just pull it up in your editor and update. Once you've got that version 1. 0. That's the new language server.

## Roadmap: Type checking

**Phil Ewels:** We talked a little bit about a couple of things which are coming down the line. And dare I ask about, variable type annotations.

**Ben Sherman:** Yes. I don't know how long it's going to take. My goal is to have a prototype in the next six months to a year, we'll see how it goes.

But this is one of those things where having a custom language, grammar, really comes in handy because now we can add type annotations, our own way. Not just however Groovy does it. That'll be sort of the level three, right? I mentioned level one, level two, level zero checking. Level three checking will be type errors. So once you fix all those name errors, just wait, there'll be a whole new slew of errors that you'll start to get.

But those errors in particular, it's gonna give you a much greater level of confidence. That everything is working, all your channels and processes are being piped into each other in the right way. And you're not calling a function with the wrong type of input or something like that.

It'll be very much closer to this idea of , if it compiles, it will run correctly. Whereas right now, it might compile correctly, but then you get some kind of runtime error halfway through your pipeline run. Those kinds of errors, the most annoying ones.

**Phil Ewels:** It's weird to be excited about getting more errors uh, I am. I think it's going to save a lot of time.

**Ben Sherman:** You get the same amount of errors. It's just whether you get them now or later after you've spent a hundred dollars.

## Conclusion

**Phil Ewels:** All right, Ben, I think I've taken enough of your time. Thank you for joining me today and running through all this stuff. It's a lot of work. It was clear from the outset. It's a huge project. It's been really difficult for me in Slack to not let on that we're working on this. I'm really excited to see what the community makes of this. We absolutely want as much feedback as possible. There, there will be rough edges. It's a new thing.

**Ben Sherman:** Yeah. And this is just the beginning. It's crazy. As much time as I've spent working on it, this is just the beginning.

**Phil Ewels:** Right. Thanks very much. We'll speak soon. And,

**Ben Sherman:** All right. Bye everybody.

**Phil Ewels:** Thanks everyone.

:::
