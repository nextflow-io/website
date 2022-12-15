title=End of DSL1 and AI pipelines
episode=9
description=Updates to the Nextflow language and putting ML pipelines to good use..
date=2023-01-11
type=podcast
subtype=News and Views
youtubeid=6WSz15_aFfQ
image=img/podcast_ep09.jpg
tags=news and views,opensource,community
status=published
author= Developer advocates
icon=logo_podcast_channels.jpg
~~~~~~

In this News and Views episode, [Phil Ewels](https://twitter.com/tallphil), [Chris Hakkaart](https://twitter.com/chris_hakk), and [Marcel Ribeiro-Dantas](https://twitter.com/mribeirodantas) discuss the hottest topics in the Nextflow world.

<!-- end-archive-description -->

#### New "fair" process directive

- New Nextflow directive added! See [nextflow-io/nextflow@`60d34cf`](https://github.com/nextflow-io/nextflow/commit/60d34cfddba721a9544908f294d008f607a9071f)
    - > This commit adds the process 'fair' directive. When fair is set to true the process outputs are guaranteed
to me emitted in the same sequence as the inputs where received instead of the first-completed-first-output semantic that's usually used by nextflow tasks

#### The end of DSL1 support!

- DSL1 now not supported in the latest edge release - `22.12.0-edge`
- You can still run old pipelines with older versions of Nextflow
    - Use `NXF_VER` before commands, this is good practice for reproducibility anyway
    - eg. `NXF_VER=22.10.4 nextflow run [...]`
- All documentation and training should now be about DSL2. Should be less confusing...
- Never a better time to convert to DSL2!

#### ChatGPT

- The community is putting it to good use ;)
    - [_""Write a song about Nextflow using "New York, I Love You but You're Bringing Me Down" as a template"_](https://twitter.com/LonsBio/status/1600266542876610560)
    - If you don't know the song, it's by LCD Soundsystem ([have a listen here](https://youtu.be/-eohHwsplvY))
- Careful with trusting it too much! It also comes out with a lot of rubbish.

#### Stable diffusion ML pipeline

- New pipeline for [Stable diffusion](https://stability.ai/blog/stable-diffusion-public-release): a [deep learning](https://en.wikipedia.org/wiki/Deep_learning), [text-to-image model](https://en.wikipedia.org/wiki/Text-to-image_model) released in 2022.
- It is primarily used to generate detailed images conditioned on text descriptions, though it can also be applied to other tasks such as [inpainting](https://en.wikipedia.org/wiki/Inpainting), outpainting, and generating image-to-image translations guided by a [text prompt](https://en.wikipedia.org/wiki/Prompt_engineering).
- Evan wrote a Nextflow pipeline ([evanfloden/stable-diffusion-nf](https://github.com/evanfloden/stable-diffusion-nf)) that incorporates Stable diffusion and uses Tower to run on AWS.
- Great example of using Tower and Nextflow to execute AI on the cloud.
- Expanding this type of application in the future.

#### Whisper pipeline

- [OpenAI Whisper](https://openai.com/blog/whisper/) is a _"neural net that approaches human level robustness and accuracy on English speech recognition"_
- Marcel wrote a Nextflow pipeline ([mribeirodantas/nf-whisper](https://github.com/mribeirodantas/nf-whisper)) to use Whisper pre-trained models to generate transcriptions / translations of audio content.
- Can fetch YouTube videos and generate transcriptions in a few minutes using GPU
- Runs on [Nextflow Tower](https://cloud.tower.nf/) and can use [Wave containers](https://seqera.io/wave/)

#### Upcoming events

- There's a [nf-core/bytesize talk on January 17th](https://nf-co.re/events/2023/bytesize_nf-core-taxprofiler) by Sofia Stamouli about the nf-core/taxprofiler pipeline
    - > A bioinformatics best-practice analysis pipeline for taxonomic profiling of shotgun metagenomic data. It allows for in-parallel profiling with multiple profiling tools against multiple databases, produces standardised output tables.
- Please sign up for the March 2023 [Nextflow / nf-core Training](https://nf-co.re/events/2023/training-march-2023)!
- Please sign up for the March 2023 [nf-core hackathon](https://nf-co.re/events/2023/hackathon-march-2023)!
