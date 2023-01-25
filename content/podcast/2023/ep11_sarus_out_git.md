title=The Sarus container engine, staging outputs and the GITHUB_TOKEN
episode=11
description=The Sarus container engine, staging outputs and the GITHUB_TOKEN.
date=2023-01-25
type=podcast
subtype=News and Views
youtubeid=WHmRP9ixLqg
image=img/podcast_ep11.jpg
tags=news and views,opensource,community
status=published
author= Developer advocates
icon=logo_podcast_channels.jpg
~~~~~~

In this News and Views episode, [Chris Hakkaart](https://twitter.com/chris_hakk) and [Marcel Ribeiro-Dantas](https://twitter.com/mribeirodantas) discuss the hottest topics in the Nextflow world.

<!-- end-archive-description -->

#### The Sarus container engine

- Recently [Sarus](https://sarus.readthedocs.io/en/stable/) was [added to Nextflow](https://github.com/nextflow-io/nextflow/pull/3470) as a container engine.
- Sarus is an OCI-compatible container engine for HPC.
- Sarus is an alternative container runtime to Docker.
    - Sarus works by converting Docker images to a common format that can then be distributed and launched on HPC systems.
- Similar deployment options to other container engines using Nextflow. 

#### New stageOutMode options

- Two new `stageOutMode` options (`rclone` and `fcp`) have been added to Nextflow.
    - The `stageOutMode` directive defines how output files are staged-out from the scratch directory to the process work directory.
- [rclone](https://rclone.org/)
    - Has been described as "The Swiss army knife of cloud storage".
    - A command-line program to manage files on cloud storage.
    - Over 40 cloud storage products support `rclone`.
- [fcp](https://github.com/Svetlitski/fcp)
    - A significantly faster alternative to the classic Unix `cp` command.
    - Handle the most common use cases of cp with much higher performance.
    - `fcp` is optimized for systems with an SSD.
- The new options will improve file movement as there are new ways for the files to be staged from the scratch directory in different circumstances.

#### GITHUB_TOKEN

- Git has become the de-facto standard for source code version control systems.
- Nextflow provides built-in support for Git and the most popular Git hosting platforms.
- Nextflow does not require any special configuration to access public repositories but requires repository credentials to access private repositories.
- The [recent Nextflow commit](https://github.com/nextflow-io/nextflow/commit/c4d39384c3536e6777bc0edcd9fd33f147aec462) adds the ability to use the `GITHUB_TOKEN` environment variable.
    - The environment variable is used when no credentials are found in the source code management configuration file (`$HOME/.nextflow/scm`).
    - It is especially useful when accessing pipeline code from a GitHub action.
    Read more about the token authentication in the [GitHub documentation](https://docs.github.com/en/actions/security-guides/automatic-token-authentication).

#### Upcoming events

- There's a [nf-core/bytesize talk on January 31](https://nf-co.re/events/2023/bytesize_transcripts) by [Franziska Bonath](https://github.com/FranBonath) about her work to generate transcripts of bytesize talks and what these might be used for in the future.
- The [Nextflow / nf-core training](https://nf-co.re/events/2023/training-march-2023) is being held March 13-16.
    - Training will be presented in different languages.
    - Registration is now open - head over to the [nf-core event](https://nf-co.re/events/2023/training-march-2023) page to sign up.
- The [nf-core hackathon](https://nf-co.re/events/2023/hackathon-march-2023) is being held March 27-29.
- Hackathon 27-29 March.
    - We will also support local hubs where local community members are hosting other community members.
    - Watch out on Slack and Twitter for an announcement that registration is open. 
