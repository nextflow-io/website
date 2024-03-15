---
title: Innovation In Science - The story behind Nextflow
date: 2015-06-09
type: post
tags: innovation,science,pipelines,nextflow
author: Maria Chatzou
icon: maria.png
---

Innovation can be viewed as the application of solutions that meet new requirements or
existing market needs. Academia has traditionally been the driving force of innovation.
Scientific ideas have shaped the world, but only a few of them were brought to market by
the inventing scientists themselves, resulting in both time and financial loses.

Lately there have been several attempts to boost scientific innovation and translation,
with most notable in Europe being the Horizon 2020 funding program. The problem with these
types of funding is that they are not designed for PhDs and Postdocs, but rather aim to
promote the collaboration of senior scientists in different institutions. This neglects two
very important facts, first and foremost that most of the Nobel prizes were given for
discoveries made when scientists were in their 20's / 30's (not in their 50's / 60's).
Secondly, innovation really happens when a few individuals (not institutions) face a
problem in their everyday life/work, and one day they just decide to do something about it
(end-user innovation). Without realizing, these people address a need that many others have.
They don’t do it for the money or the glory; they do it because it bothers them!
Many examples of companies that started exactly this way include Apple, Google, and
Virgin Airlines.

### The story of Nextflow

Similarly, Nextflow started as an attempt to solve the every-day computational problems we
were facing with “big biomedical data” analyses. We wished that our huge and almost cryptic
BASH-based pipelines could handle parallelization automatically. In our effort to make that
happen we stumbled upon the [Dataflow](http://en.wikipedia.org/wiki/Dataflow_programming)
programming model and Nextflow was created.
We were getting furious every time our two-week long pipelines were crashing and we had
to re-execute them from the beginning. We, therefore, developed a caching system, which
allows Nextflow to resume any pipeline from the last executed step. While we were really
enjoying developing a new [DSL](http://en.wikipedia.org/wiki/Domain-specific_language) and
creating our own operators, at the same time we were not willing to give up our favorite
Perl/Python scripts and one-liners, and thus Nextflow became a polyglot.

Another problem we were facing was that our pipelines were invoking a lot of
third-party software, making distribution and execution on different platforms a nightmare.
Once again while searching for a solution to this problem, we were able to identify a
breakthrough technology [Docker](https://www.docker.com/), which is now revolutionising
cloud computation. Nextflow has been one of the first framework, that fully
supports Docker containers and allows pipeline execution in an isolated and easy to distribute manner.
Of course, sharing our pipelines with our friends rapidly became a necessity and so we had
to make Nextflow smart enough to support [Github](https://github.com) and [Bitbucket](https://bitbucket.org/) integration.

I don’t know if Nextflow will make as much difference in the world as the Dataflow
programming model and Docker container technology are making, but it has already made a
big difference in our lives and that is all we ever wanted…

### Conclusion

Summarising, it is a pity that PhDs and Postdocs are the neglected engine of Innovation.
They are not empowered to innovate, by identifying and addressing their needs, and to
potentially set up commercial solutions to their problems. This fact becomes even sadder
when you think that only 3% of Postdocs have a chance to become PIs in the UK. Instead more
and more money is being invested into the senior scientists who only require their PhD students
and Postdocs to put another step into a well-defined ladder. In todays world it seems that
ideas, such as Nextflow, will only get funded for their scientific value, not as innovative
concepts trying to address a need.
