---
title: Application of Nextflow and nf-core to ancient environmental eDNA
date: 2024-04-17
type: post
description: In this blog post, James goes through a workshop he organized to demonstrate the efficiency and reproducibility Nextflow and nf-core can bring to analyze ancient environmental DNA.
image: /img/blog-2024-04-17--share.png
tags: nextflow,nf-core,workshop,ambassador_post
status: published
author: James Fellows Yates
icon: james.jpeg
community_post: true
ambassador_post: true
---

Ancient environmental DNA (eDNA) is currently a hot topic in archaeological, ecological, and metagenomic research fields. Recent eDNA studies have shown that authentic ‘ancient’ DNA can be recovered from soil and sediments even as far back as 2 million years ago<sup>1</sup>. However, as with most things metagenomics (the simultaneous analysis of the entire DNA content of a sample), there is a need to work at scale, processing the large datasets of many sequencing libraries to ‘fish’ out the tiny amounts of temporally degraded ancient DNA from amongst a huge swamp of contaminating modern biomolecules.

<!-- end-archive-description -->

This need to work at scale, while also conducting reproducible analyses to demonstrate the authenticity of ancient DNA, lends itself to the processing of DNA with high-quality pipelines and open source workflow managers such as Nextflow. In this context, I was invited to the Australian Center for Ancient DNA (ACAD) at the University of Adelaide in February 2024 to co-teach a graduate-level course on ‘Hands-on bioinformatics for ancient environmental DNA’, alongside other members of the ancient eDNA community. Workshop participants included PhD students from across Australia, New Zealand, and even from as far away as Estonia.

<div style="margin-top: 2rem; margin-bottom: 2rem;">
    <img src="/img/blog-2024-04-17-img1a.jpg" alt="Mentor compliment about new module added" />
    © Photo: Peter Mundy and Australian Center for Ancient DNA
</div>

We began the five-day workshop with an overview of the benefits of using workflow managers and pipelines in academic research, which include efficiency, portability, reproducibility, and fault-tolerance, and we then proceeded to introduce the Ph.D. students to installing Nextflow, and configure pipelines for running on different types of computing infrastructure.

<div style="margin-top: 2rem; margin-bottom: 2rem;">
    <img src="/img/blog-2024-04-17-img1b.jpg" alt="Review comment in GitHub" />
    © Photo: Peter Mundy and Australian Center for Ancient DNA
</div>

Over the next two days, I then introduced two well-established nf-core pipelines: [nf-core/eager](https://nf-co.re/eager)<sup>2</sup> and [nf-core/mag](https://nf-co.re/mag)<sup>3</sup>, and explained to students how these pipelines can be applied to various aspects of environmental metagenomic and ancient DNA analysis:
nf-core/eager is a dedicated ‘swiss-army-knife’ style pipeline for ancient DNA analysis that performs genetic data preprocessing, genomic alignment, variant calling, and metagenomic screening with specific tools and parameters to account for the characteristics of degraded DNA.
nf-core/mag is a best-practice pipeline for metagenomic de novo assembly of microbial genomes that performs preprocessing, assembly, binning, bin-refinement and validation. It also contains a specific subworkflow for the authentication of ancient contigs.
In both cases, the students of the workshops were given practical tasks to set up and run both pipelines on real data, and time was spent exploring the extensive nf-core documentation and evaluating the outputs from MultiQC, both important components that contribute to the quality of nf-core pipelines.

The workshop was well received by students, and many were eager (pun intended) to start running Nextflow and nf-core pipelines on their own data at their own institutions.

I would like to thank Vilma Pérez at ACAD for the invitation to contribute to the workshop as well as Mikkel Pedersen for being my co-instructor, and the nf-core community for continued support in the development of the pipelines. Thank you also to Tina Warinner for proof-reading this blog post, and I would like to acknowledge [ACAD](https://www.adelaide.edu.au/acad/), the [University of Adelaide Environment Institute](https://www.adelaide.edu.au/environment/), the [Werner Siemens-Stiftung](https://www.wernersiemens-stiftung.ch/), [Leibniz HKI](https://www.leibniz-hki.de/), and [MPI for Evolutionary Anthropology](https://www.eva.mpg.de) for financial support to attend the workshop and support in developing nf-core pipelines.

---

<sup>1</sup> Kjær, K.H., Winther Pedersen, M., De Sanctis, B. et al. A 2-million-year-old ecosystem in Greenland uncovered by environmental DNA. Nature **612**, 283–291 (2022). [https://doi.org/10.1038/s41586-022-05453-y](https://doi.org/10.1038/s41586-022-05453-y)

<sup>2</sup> Fellows Yates, J.A., Lamnidis, T.C., Borry, M., Andrades Valtueña, A., Fagernäs, Z., Clayton, S., Garcia, M.U., Neukamm, J., Peltzer, A.. Reproducible, portable, and efficient ancient genome reconstruction with nf-core/eager. PeerJ 9:10947 (2021) [http://doi.org/10.7717/peerj.10947](http://doi.org/10.7717/peerj.10947)

<sup>3</sup> Krakau, S., Straub, D., Gourlé, H., Gabernet, G., Nahnsen, S., nf-core/mag: a best-practice pipeline for metagenome hybrid assembly and binning, NAR Genomics and Bioinformatics, **4**:1 (2022) [https://doi.org/10.1093/nargab/lqac007](https://doi.org/10.1093/nargab/lqac007)
