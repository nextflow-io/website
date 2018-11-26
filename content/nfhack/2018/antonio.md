title=Nextflow for supporting the European Reference Laboratory for GMO Food and Feed
date=2018-11-15
type=col8
tags=nextflow,nfhack18,workshop
status=published
~~~~~~

## Nextflow for supporting the European Reference Laboratory for GMO Food and Feed

### Antonio Puertas Gallardo
*Bioinformatician, Directorate General Joint Research Centre, European Commission* 

The European Union Reference Laboratory for Genetically Modified Food and Feed (EU-RL GMFF), hosted by the Joint Research Centre (JRC) of the European Commission, must, pursuant to Article 32 of Regulation (EC) N. 882/2004, provide National Reference Laboratories with reference methods and tools for GMO analysis. Commission Regulation (EC) No 641/2004 requires the EU-RL GMFF to maintain a database containing GMO events sequence. The polymerase chain reaction (PCR) has proven to be the most accurate and reliable technique available for GMO detection, identification and quantification and is applicable to a wide range of samples, from seeds to highly processed food and feed. The validation of PCR-based detection methods is a fundamental task that requires integrating the combination of both the experimental approach and of bioinformatics analyses for sequence similarity searches. The EU-RL GMFF is requested to perform in silico validation of the proposed detection methods, with respect to specificity.

In order to fulfil this duty, the JRC has implemented a tool called METSCAN that allows the performance of complex bioinformatics analyses on the specificity of PCR-based detection methods. METSCAN relies on the power of a High Performance Computing cluster which runs HTCondor scheduler, through a simple user interface. METSCAN has the objective to make in silico predictions on the detection methods  specificity, to direct recommendations on the need or not for experimental testing of method specificity and, in case, to define what source of DNA, e.g. vector, plant species or other GMOs can cross-react with each detection method.

We are now exploring with Nexflow a new framework based on new workflows in order to improve MetScan and make the complete tool/pipeline independent of the local HW/SW platform. This idea might be a starting point for assessing the validation of Bioinformatics pipelines under the European Regulation, as part of the concept of Regulatory Bioinformatics, which is proposed to be applied into other fields like clinical trials and medicine workflows.

### Deck

<a href='/misc/nfhack18/antonio.pdf'><img src='/img/deck.png' width='60pt' /></a>

### Bio 

Born and raised in Argentina, Antonio Puertas Gallardo has studied Electronic Engineering.
He arrived to the European Commission in 2003 as IT System Manager. He was responsible for Joint 
Research Centre (JRC) data centre since 2005 to 2008. Then he was in charge for JRC Central Storage infrastructure since 2006 to 2012. He was also IT System manager and responsible for HPC in the JRC 
since 2008 to  2013.

In 2014 he moved to the JRC Bioinformatics team. Today he is  dealing with HPC
support to Bioinformatics team, doing working on metagenomics analysis and lately
he started to work with Natural Language Processing and Machine Understanding.



### More information 

The event program is available at [this link](https://github.com/nextflow-io/nf-hack18/blob/master/schedule.md). For registration and other information check it out [this page](http://www.crg.eu/en/event/coursescrg-nextflow-reproducible-silico-genomics-0).