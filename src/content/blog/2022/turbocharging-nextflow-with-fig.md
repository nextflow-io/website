---
title: Turbo-charging the Nextflow command line with Fig!
date: 2022-09-22
type: post
description: Powerful autocomplete for Nextflow in your command line interface with Fig.
image: /img/turbocharging_nextflow_with_fig.jpg
tags: nextflow,development,learning
author: Marcel Ribeiro-Dantas
icon: marcel.png
---

Nextflow is a powerful workflow manager that supports multiple container technologies, cloud providers and HPC job schedulers. It shouldn't be a surprise that wide ranging functionality leads to a complex interface, but comes with the drawback of many subcommands and options to remember. For a first-time user (and sometimes even for some long-time users) it can be difficult to remember everything. This is not a new problem for the command-line; even very common applications such as grep and tar are famous for having a bewildering array of options.

![xkcd charge making fun of tar tricky command line arguments](/img/xkcd_tar_charge.png)
https://xkcd.com/1168/

Many tools have sprung up to make the command-line more user friendly, such as tldr pages and rich-click. [Fig](https://fig.io) is one such tool that adds powerful autocomplete functionality to your terminal. Fig gives you graphical popups with color-coded contexts more dynamic than shaded text for recent commands or long blocks of text after pressing tab.

Fig is compatible with most terminals, shells and IDEs (such as the VSCode terminal), is fully supported in MacOS, and has beta support for Linux and Windows. In MacOS, you can simply install it with `brew install --cask fig` and then running the `fig` command to set it up.

We have now added Nextflow for Fig. Thanks to Figs open source core we were able to contribute specifications in Typescript that will now be automatically added for anyone installing or updating Fig. Now, with Fig, when you start typing your Nextflow commands, you’ll see autocomplete suggestions based on what you are typing and what you have typed in the past, such as your favorite options.

![GIF with a demo of nextflow log/list subcommands](/img/nxf-log-list-params.gif)

The Fig autocomplete functionality can also be adjusted to suit our preferences. Suggestions can be displayed in alphabetical order or as a list of your most recent commands. Similarly,  suggestions can be displayed all the time or only when you press tab.

The Fig specification that we've written not only suggests commands and options, but dynamic inputs too. For example, finding previous run names when resuming or cleaning runs is tedious and error prone. Similarly, pipelines that you’ve already downloaded with `nextflow pull` will be autocompleted if they have been run in the past. You won't have to remember the full names anymore, as Fig generators in the autocomplete allow you to automatically complete the run name after typing a few letters where a run name is expected. Importantly, this also works for pipeline names!

![GIF with a demo of nextflow pull/run/clean/view/config subcommands](/img/nxf-pull-run-clean-view-config.gif)

Fig for Nextflow will make you increase your productivity regardless of your user level. If you run multiple pipelines during your day you will immediately see the benefit of Fig. Your productivity will increase by taking advantage of this autocomplete function for run and project names. For Nextflow newcomers it will provide an intuitive way to explore the Nextflow CLI with built-in help text.

While Fig won’t replace the need to view help menus and documentation it will undoubtedly save you time and energy searching for commands and copying and pasting run names. Take your coding to the next level using Fig!
