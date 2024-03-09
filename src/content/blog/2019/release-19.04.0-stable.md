---
title: Nextflow 19.04.0 stable release is out!
date: 2019-04-18
type: post
tags: nextflow,release,stable
author: Paolo Di Tommaso
icon: paolo.jpg
---

We are excited to announce the new Nextflow 19.04.0 stable release!

This version includes numerous bug fixes, enhancement and new features.

#### Rich logging

In this release, we are making the new interactive rich output using ANSI escape characters as the default logging option. This produces a much more readable and easy to follow log of the running workflow execution.

<script type="text/javascript" src="https://asciinema.org/a/IrT6uo85yyVoOjPa6KVzT2FXQ.js" id="asciicast-IrT6uo85yyVoOjPa6KVzT2FXQ" async></script>

The ANSI log is implicitly disabled when the nextflow is launched in the background i.e. when using the `-bg` option. It can also be explicitly disabled using the `-ansi-log false` option or setting the `NXF_ANSI_LOG=false` variable in your launching environment.

#### NCBI SRA data source

The support for NCBI SRA archive was introduced in the [previous edge release](/blog/2019/release-19.03.0-edge.html). Given the very positive reaction, we are graduating this feature into the stable release for general availability.

#### Sharing

This version includes also a new Git repository provider for the [Gitea](https://gitea.io) self-hosted source code management system, which is added to the already existing support for GitHub, Bitbucket and GitLab sharing platforms.

#### Reports and metrics

Finally, this version includes important enhancements and bug fixes for the task executions metrics collected by Nextflow. If you are using this feature we strongly suggest updating Nextflow to this version.

Remember that updating can be done with the `nextflow -self-update` command.

### Changelog

The complete list of changes and bug fixes is available on GitHub at [this link](https://github.com/nextflow-io/nextflow/releases/tag/v19.04.0).

### Contributions

Special thanks to all people contributed to this release by reporting issues, improving the docs or submitting (patiently) a pull request (sorry if we have missed somebody):

- [Alex Cerjanic](https://github.com/acerjanic)
- [Anthony Underwood](https://github.com/aunderwo)
- [Akira Sekiguchi](https://github.com/pachiras)
- [Bill Flynn](https://github.com/wflynny)
- [Jorrit Boekel](https://github.com/glormph)
- [Olga Botvinnik](https://github.com/olgabot)
- [Ólafur Haukur Flygenring](https://github.com/olifly)
- [Sven Fillinger](https://github.com/sven1103)
