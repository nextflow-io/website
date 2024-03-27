---
title: Nextflow's colorful new console output
date: 2024-05-05
type: post
description: Bringing a splash of color to your pipeline monitoring.
image: img/blog-2024-03-14--share.jpg
tags: nextflow
status: published
author: Phil Ewels
icon: phil.jpg
---

Nextflow is a command-line interface (CLI) tool that runs in the terminal. Everyone who has launched Nextflow from the command line knows what it’s like to follow the console output as a pipeline runs: the excitement of watching jobs zipping off as they’re submitted, the satisfaction of the phrase _“Pipeline completed successfully!”_ and occasionally, the sinking feeling of seeing an error message.

Because the CLI is the primary way that people interact with Nextflow, a little bit of polish can have a big effect. In this article, I’m excited to describe an upgrade for the console output that should make monitoring workflow progress just a little easier.

The new functionality is available in `24.02-0-edge` and will be included in the next `24.04.0` stable release. You can try it out now by updating Nextflow as follows:

```bash
NXF_EDGE=1 nextflow self-update
```

## Background

The Nextflow console output hasn’t changed much over the 10 years that it’s been around. The biggest update happened in 2018 when “ANSI logging” was released in version `18.10.0`. This replaced the stream of log messages announcing each task submission with a view that updates dynamically, giving an overview of each process. This gives an overview of the pipeline’s progress rather than being swamped with thousands of individual task submissions.

<figure>
  <img src="/img/blog-nextflow-colored-logs/nextflow_log_with_without_ansi.png" alt="Nextflow console output with and without ANSI logging">
  <figcaption>
    ANSI console output. Nextflow log output from running the nf-core/rnaseq pipeline before (Left) and after (Right) enabling ANSI logging.
  </figcaption>
</figure>

I can be a little obsessive tool user interfaces. The nf-core template, as well as MultiQC and nf-core/tools all have coloured terminal output, mostly using the excellent [textualize/rich](https://github.com/Textualize/rich). I’ve also written a couple of general-use tools around this such as [ewels/rich-click](https://github.com/ewels/rich-click/) for Python CLI help texts, and [ewels/rich-codex](https://github.com/ewels/rich-codex) to auto-generate screenshots from code / commands in markdown. The problem with being surrounded by so much colored CLI output is that any tools _without_ colors start to stand out. Dropping hints to the Nextflow team didn’t work, so eventually I whipped up [a proposal](https://github.com/nextflow-io/nextflow/issues/3976) of what the console output could look like using the tools I knew: Python and Rich. Paolo knows me well and [offered up a bait](https://github.com/nextflow-io/nextflow/issues/3976#issuecomment-1568071479) that I couldn’t resist: _“Phil. I think this a great opportunity to improve your Groovy skills 😆”._

## Showing what’s important

The console output shown by Nextflow describes a range of information. Much of it aligns in vertical columns, but not all. There’s also a variety of fields, some of which are more important than others to see at a glance.

<figure>
  <img src="/img/blog-nextflow-colored-logs/nextflow_coloured_logs.png" alt="New coloured output from Nextflow">
  <figcaption>
    Introducing: colored console output. Output from running nf-core/rnaseq with the new colors applied (nf-core header removed for clarity).
  </figcaption>
</figure>

With some judicious use of the `dim` style, we can make less important information fade into the background. For example, the “stem” of the fully qualified process identifiers now step back to allow the process name to stand out. Secondary information such as the number of tasks that were cached, or the executor that is being submitted to, are still there to see but take a back seat. Doing the reverse with some `bold` text helps to highlight the run name – key information for identifying and resuming pipeline runs. Using color allows different fields to be easily distinguished, such as process labels and task hashes. Greens, blues and reds in the task statuses allow a reader to get an impression of the run progress without needing to read every number.

Probably the most difficult aspect technically was the `NEXTFLOW` header line. I knew I wanted to use the “Nextflow green” here, or as close to it as possible. But colors in the terminal are tricky. What the ANSI standard defines as `green`, `black` and `blue` can vary significantly across different systems and terminal themes. Some people use a light-color scheme and others run in dark mode. This hadn’t mattered much for most of the colors up until this point, I could use the [Jansi](https://github.com/fusesource/jansi) library to use named colors and they should look ok. But for the specific RGB of the _Nextflow Green_ I had to [hardcode specific ANSI control characters](https://github.com/nextflow-io/nextflow/blob/c9c7032c2e34132cf721ffabfea09d893adf3761/modules/nextflow/src/main/groovy/nextflow/cli/CmdRun.groovy#L379-L389). But it got worse - it turns out that the default Terminal app that ships with OS X only supports 256 colors, so I had to find the closest match (_“light sea green”_ if you’re curious). Even once the green was ok, using `black` as the text color meant that it would actually render as white with some terminal color themes and be unreadable. In the end, the header text is a very dark gray.

<figure>
  <img src="/img/blog-nextflow-colored-logs/testing_terminal_themes.png" alt="Testing many horrible terminal themes">
  <figcaption>
    Testing color rendering across a wide range of themes in the OS X Terminal app.
  </figcaption>
</figure>

## More than just colors

Whilst the original intent was focused on using color, it didn’t take long to come up with a shortlist of other niggles that I wanted to fix. I took this project as an opportunity to address a few of these, specifically:

- Make the most of the available width in the terminal
  - Redundant text is now cut down when the screen is narrow. Specifically the repeated `process >` text, plus other small gains such as replacing the three `...` characters with a single `…` character. The percentage-complete is removed if the window is really narrow. These changes happen dynamically every time the screen refreshes, so should update if you resize the terminal window.
- Be more selective about which part of process names are truncated
  - There’s only so much width that can be saved, and fully qualified process names are long. The current Nextflow console output truncates the end of the identifier if there’s no space, but this is the part that varies most between pipeline steps. Instead, we can truncate the start and preserve the process name and label.
- Don’t show all pending processes without tasks
  - The existing ANSI logging shows _all_ processes in the pipeline, even those that haven’t had any tasks submitted. If a pipeline has a lot of processes this can push the running processes out of view.
  - Nextflow now tracks the number of available rows in the terminal and hides pending processes once we run out of space. Running processes are always printed.

The end result is console output that makes the most of the available space in your terminal window:

<figure>
  <img src="/img/blog-nextflow-colored-logs/nextflow_console_varying_widths.png" alt="Nextflow console output at different terminal window widths">
  <figcaption>
    Progress of the nf-core/rnaseq shown across 3 different terminal-width breakpoints, with varying levels of text truncation.
  </figcaption>
</figure>

## Contributing to Nextflow

Despite building tools that use Nextflow for many years, I’ve spent relatively little time venturing into the main codebase myself. Just as with any contributor, part of the challenge was figuring out how to build Nextflow, how to navigate its code structure and how to write tests. I found it quite a fun experience, so I described and demoed the process in a recent nf-core Bytesize talk titled “[Contributing to Nextflow](https://nf-co.re/events/2024/bytesize_nextflow_dev)”. You can watch the talk on [YouTube](https://www.youtube.com/watch?v=R0fqk5OS-nw), where I explain the mechanics of forking Nextflow, enhancing, compiling, and testing changes locally, and contributing enhancements back to the main code base.

## But wait, there’s more!

I’m happy with how the new console output looks, and it seems to have been well received so far. But once the warm glow of the newly merged pull-request started to subside, I realized that there was more to do. The console output is great for monitoring a running pipeline, but I spend most of my time these days digging through much more verbose `.nextflow.log` files. Suddenly it seemed a little unfair that these didn’t also benefit from a similar treatment.

This project was a little different because the logs are just files on the disk, meaning that I could approach the problem with whatever code stack I liked. Coincidentally, [Will McGugan](https://github.com/willmcgugan) (author of [textualize/rich](https://github.com/Textualize/rich)) was recently [writing about](https://textual.textualize.io/blog/2024/02/11/file-magic-with-the-python-standard-library/) a side project of his own: [Toolong](https://github.com/textualize/toolong). This is a terminal app built using [Textual](https://www.textualize.io/) which is specifically aimed at viewing large log files. I took it for a spin and it did a great job with Nextflow log files right out of the box, but I figured that I could take it further. At its core, Toolong uses the [Rich](https://github.com/textualize/rich) library to format text and so with a little hacking, I was able to introduce a handful of custom formatters for the Nextflow logs. And voila, we have colored console output for log files too!

<figure>
  <img src="/img/blog-nextflow-colored-logs/nextflow_logs_side_by_side.png" alt="Formatting .nextflow.log files with Toolong">
  <figcaption>
    The tail end of a `.nextflow.log` file, rendered with `less` (Left) and Toolong (Right). Try finding the warning log message in both!
  </figcaption>
</figure>

By using Toolong as a viewer we get much more than just syntax highlighting too - it provides powerful file navigation and search functionality. It also supports tailing files in real time, so you can launch a pipeline in one window and tail the log in another to have the best of both worlds!

<figure>
  <video controls>
    <source src="/img/blog-nextflow-colored-logs/nextflow_logs_best_both_worlds.mp4" type="video/mp4">
    Your browser does not support the video tag.
  </video>
  <figcaption>
    Running nf-core/rnaseq with the new Nextflow coloured console output (Left) whilst simultaneously tailing the `.nextflow.log` file using `nf-core log` (Right).
  </figcaption>
</figure>

This work with Toolong is still in two [open](https://github.com/Textualize/toolong/pull/47) [pull-requests](https://github.com/nf-core/tools/pull/2895) as I write this, but hopefully you’ll soon be able to use the `nf-core log` command in a directory where you’ve run Nextflow, and it’ll launch Toolong with any log files it finds.
