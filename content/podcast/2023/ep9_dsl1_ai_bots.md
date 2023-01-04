title=End of DSL1, AI and chatting to bots
episode=9
description=End of DSL1, AI and chatting to bots
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

#### The new 'fair' directive

- Addition of the new 'fair' directive.
- Not the traditional use of FAIR.
- When `fair` is set to true the process outputs are ordered in the same sequence as the inputs were received instead of the first-completed-first-output semantic that's usually used.

#### The end of DSL1

- DSL1 is now not supported in the latest edge release - 22.12.0-edge.
- Old pipelines can still be run with older versions of Nextflow.
- `NXF_VER` can be used before commands to set the version of Nextflow, good practice for reproducibility.
- All documentation and training will be about DSL2.
- Dropping DSL1 support should be less confusing for everyone.
- It has never been a better time to convert to DSL2!

#### Chatting with bots

- Applications of AI are appearing more and more in social media and research.
- [ChatGPT](https://openai.com/blog/chatgpt/).
  - ChatGPT is a chatbot launched by OpenAI.
  - Nextflow was featured in a tweet from [Andrew Lonsdale](https://twitter.com/LonsBio/status/1600266542876610560) showing a ChatGPT generated song about Nextflow using "New York, I Love You but You're Bringing Me Down" as a template.
- [Stable diffusion](https://stability.ai/blog/stable-diffusion-public-release).
  - Stable Diffusion is a deep learning, text-to-image model launched by Stability AI.
  - [Evan Floden](https://twitter.com/evanfloden) wrote [stable-diffusion-nf](https://github.com/evanfloden/stable-diffusion-nf), a proof of concept Nextflow pipeline for submitting stable-diffusion jobs to the cloud using Tower.
- [OpenAI Whisper](https://openai.com/blog/whisper/).
  - OpenAI Whisper is a tool for generating transcripts from audio content launched by OpenAI.
  - [Marcel](https://twitter.com/mribeirodantas) wrote [nf-whisper](https://github.com/mribeirodantas/nf-whisper), a Nextflow pipeline for downloading and transcribing audio content from YouTube videos using OpenAI Whisper.
- Many other great examples of AI and Nextflow exist and will likely emerge in the future.
