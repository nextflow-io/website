title=Clarification about the Nextflow license
date=2018-07-20
type=post
tags=nextflow,gpl
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~


Over  past week there was some discussion on social media regarding the Nextflow license 
and its impact on users' workflow applications. 

<blockquote class="twitter-tweet tw-align-center" data-lang="en"><p lang="en" dir="ltr">… don’t use Nextflow, yo. <a href="https://t.co/Paip5W1wgG">https://t.co/Paip5W1wgG</a></p>&mdash; Konrad Rudolph 👨‍🔬💻 (@klmr) <a href="https://twitter.com/klmr/status/1016606226103357440?ref_src=twsrc%5Etfw">July 10, 2018</a></blockquote>
<script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

<blockquote class="twitter-tweet tw-align-center" data-lang="en"><p lang="en" dir="ltr">This is certainly disappointing. An argument in favor of writing workflows in <a href="https://twitter.com/commonwl?ref_src=twsrc%5Etfw">@commonwl</a>, which is independent of the execution engine. <a href="https://t.co/mIbdLQQxmf">https://t.co/mIbdLQQxmf</a></p>&mdash; John Didion (@jdidion) <a href="https://twitter.com/jdidion/status/1016612435938160640?ref_src=twsrc%5Etfw">July 10, 2018</a></blockquote>
<script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>


<blockquote class="twitter-tweet tw-align-center" data-lang="en"><p lang="en" dir="ltr">GPL is generally considered toxic to companies due to fear of the viral nature of the license.</p>&mdash; Jeff Gentry (@geoffjentry) <a href="https://twitter.com/geoffjentry/status/1016656901139025921?ref_src=twsrc%5Etfw">July 10, 2018</a></blockquote>
<script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>


### What's the problem with GPL?

Nextflow has been released under the GPLv3 license since its early days [over 5 years ago](https://github.com/nextflow-io/nextflow/blob/c080150321e5000a2c891e477bb582df07b7f75f/src/main/groovy/nextflow/Nextflow.groovy). 
GPL is a very popular open source licence used by many projects 
(like, for example, [Linux](https://www.kernel.org/doc/html/v4.17/process/license-rules.html) and [Git](https://git-scm.com/about/free-and-open-source)) 
and it has been designed to promote the adoption and spread of open source software and culture. 

With this idea in mind, GPL requires the author of a piece of software, *derived* from a GPL licensed application or library, to distribute it using the same license i.e. GPL itself.

This is generally good, because this requirement incentives the growth of the open source ecosystem and the adoption of open source software more widely. 

However, this is also a reason for concern by some users and organizations because it's perceived as too strong requirement by copyright holders (who may not want to disclose their code) and because it can be difficult to interpret what a \*derived\* application is.  See for example 
[this post by Titus Brown](http://ivory.idyll.org/blog/2015-on-licensing-in-bioinformatics.html) at this regard. 

#### What's the impact of the Nextflow license on my application?

If you are not distributing your application, based on Nextflow, it doesn't affect you in any way. 
If you are distributing an application that requires Nextflow to be executed, technically speaking your application is dynamically linking to the Nextflow runtime and it uses routines provided by it. For this reason your application should be released as GPLv3. See [here](https://www.gnu.org/licenses/gpl-faq.en.html#GPLStaticVsDynamic) and [here](https://www.gnu.org/licenses/gpl-faq.en.html#IfInterpreterIsGPL).

<b>However, this was not our original intention. We don’t consider workflow applications to be subject to the GPL copyleft obligations of the GPL even though they may link dynamically to Nextflow functionality through normal calls and we are not interested to enforce the license requirement to third party workflow developers and organizations. Therefore you can distribute your workflow application using the license of your choice. For other kind of derived applications the GPL license should be used, though.
</b>

### That's all? 

No. We are aware that this is not enough and the GPL licence can impose some limitation in the usage of Nextflow to some users and organizations. For this reason we are working with the CRG legal department to move Nextflow to a more permissive open source license. This is primarily motivated by our wish to make it more adaptable and compatible with all the different open source ecosystems, but also to remove any remaining legal uncertainty that using Nextflow through linking with its functionality may cause. 

We are expecting that this decision will be made over the summer so stay tuned and continue to enjoy Nextflow.
