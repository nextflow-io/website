title=S3 file names, Spack package manager support and local hackathon sites.
episode=12
description=S3 file names, Spack package manager support and local hackathon sites.
date=2023-02-01
type=podcast
subtype=News and Views
youtubeid=PZGhZ9LZM9w
image=img/podcast_ep12.jpg
tags=news and views,opensource,community
status=published
author= Developer advocates
icon=logo_podcast_channels.jpg
~~~~~~

In this News and Views episode, [Chris Hakkaart](https://twitter.com/chris_hakk) and [Marcel Ribeiro-Dantas](https://twitter.com/mribeirodantas) discuss the hottest topics in the Nextflow world.

<!-- end-archive-description -->

#### Improved handling of S3 paths with spaces

- AWS jobs would fail to resume a job when an `S3` file path contained a space.
    - The issue didn't occur on local servers where spaces were.
- A better implementation of the `S3Serializer` was added and has resolved this problem.
- The Nextflow documentation was also [updated](https://github.com/nextflow-io/nextflow/commit/2fe2322773dc943d52160e8fa8331b91ba74268c) to make it clear that `toUriString`, not `toString`, should be used to retain the protocol schema for path strings.
    - The protocol schema is important when using `S3` paths.

#### The Spack  package manager

- Support for the [Spack](https://spack.io) package manager has been [added](https://github.com/nextflow-io/nextflow/pull/3580) to Nextflow.
- Spack is a package manager for supercomputers, Linux, and macOS.
    - It makes installing software easy and supports multiple versions, configurations, platforms, and compilers.
- Nextflow now has built-in support for Spack that allows the configuration of workflow dependencies using Spack recipes and environment files.
    - This allows Nextflow applications to build packages from source on the compute infrastructure whilst taking advantage of the configuration flexibility provided by Nextflow.
- It is used in a similar way to using conda environments with its own scope and directives.

#### Local sites for the nf-core hackathon

- Hacking with others is a lot of fun and local sites are our way of replicating the fun we had at the in-person hackathon prior to the 2022 Nextflow Summit.
- Everything will still be online in Gather so you can join from anywhere.
- There are more than [15 local sites](https://nf-co.re/events/2023/hackathon-march-2023#local-events) that are being hosted by community members.
    - Thanks to funding from the Chan Zuckerberg Initiative and Seqera Labs all sites are being supported with goodies.
    - We are especially excited to have a site hosted by Google at the [Google Academy](https://nf-co.re/events/2023/hackathon-march-2023/uk-google.md) in London.
    - Just like the last hackathon, an [nf-core training event](https://nf-co.re/events/2023/training-march-2023) is being held prior to the hackathon for new community members to upskill.

#### Upcoming events

- The next [nf-core/bytesize talk on February 7](https://nf-co.re/events/2023/bytesize_precommit).
    - [Matthias Hörtenhuber](https://github.com/mashehu) is going to explain the use of the newly added nf-core `pre-commit` tool.
    - `pre-commit` was developed to inspect the snapshot that is about to be committed and helps to check the code before adding it to the repository.
- The Nextflow and nf-core training is being held March 13-16.
    - Training will be presented in different languages.
    - Registration is [open now](https://nf-co.re/events/2023/training-march-2023).
- The nf-core hackathon is being held March 27-29.
    - Registration is [open now](https://nf-co.re/events/2023/hackathon-march-2023).
