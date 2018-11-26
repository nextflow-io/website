title=Industrial Personalised Immunotherapy Pipeline Development with Nextflow
date=2018-10-11
type=col8
tags=nextflow,nfhack18,workshop
status=published
~~~~~~

## Industrial Personalised Immunotherapy Pipeline Development with Nextflow

### Luke Goodsell
*Team Leader, Bioinformatics Software Development, Achilles Therapeutics, UK* 

Cancer immunotherapy is a promising means of treatment that uses the patient's own immune system to destroy cancer cells. For over three decades, labs have been exploring ways to identify tumour specific surface proteins (neoantigens) and use them to induce the patient's immune cells to specifically target their cancer. This process is labour intensive and has highly variable results, in part due to tumour heterogeneity.
 
[Achilles Therapeutics](https://achillestx.com/) is developing a process for delivering a patient-specific immunotherapy that exploits recent research into tumour evolution to reliably identify and target the neoantigens that occur early in a tumour's development and thus are present in every tumour cell (i.e. clonal). At the heart of this process is a bioinformatic pipeline that takes DNA- and RNA-sequencing data and reports clonal, actionable neoantigens. Achilles inherited an R-based pipeline suited for answering cutting-edge research questions from one of its academic founding labs, but adapting this for industrial use posed many challenges including: adaptations for high-throughput operation, guaranteeing reproducibility, portability to other compute environments (particularly AWS), conformance to required standards for clinical use, end-to-end automation and quality control, "fail-safe" operation and rapid prototyping of new features.
 
This presentation will describe how the Achilles bioinformatics team used Nextflow (and other tools) to adapt a research-focussed sequence-processing pipeline into an efficient, reliable, portable and risk-managed system suitable for delivering a clinical personalised immunotherapy. This will include tips and tricks for developing modular Nextflow workflows and for ensuring early detection of errors.

### Deck

<a href='/misc/nfhack18/luke.pdf'><img src='/img/deck.png' width='45pt' /></a>

### Bio 

[Luke Goodsell](https://www.linkedin.com/in/luke-goodsell-910a2793/) leads the Bioinformatics Software Development team at Achilles Therapeutics, where he has worked for 2 years. Before that he worked as a Computational Biologist at Oxford Gene Technology after completing his PhD in Structural, Computational and Chemical Biology. He is interested in applying data analysis techniques and software engineering experience to identify answers to bioinformatics questions and implementing robust, dependable software that uses those findings.

### More information 

The event program is available at [this link](https://github.com/nextflow-io/nf-hack18/blob/master/schedule.md). For registration and other information check it out [this page](http://www.crg.eu/en/event/coursescrg-nextflow-reproducible-silico-genomics-0).