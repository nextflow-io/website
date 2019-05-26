title=Pay As You Go Cloud Bioinformatics for Pathogens
date=2019-05-28
type=col8
tags=nextflow,nfcamp,2019,workshop
status=published
~~~~~~

## Pay As You Go Cloud Bioinformatics for Pathogens

### Anthony Underwood
*Bioinformatics Implementation Manager, Centre for Genomic Pathogen Surveillance, Wellcome Trust Sanger Institute, UK* 

### Ben Taylor
*Senior Software Developer, Centre for Genomic Pathogen Surveillance, Wellcome Trust Sanger Institute, UK* 

Nextflow provides a mechanism to develop high throughput parallelizable pipelines on a small desktop machine and then scale out to larger infrastructures such as HPC clusters. If you do not have a cluster Nextflow allows the same pipeline to run in the cloud such as on AWS Batch where hundreds or thousands of jobs can run in parallel. However: 

1. it is not trivial to set up
2. users need to be comfortable on the command line
3. and there’s a risk you could rack up large bills

We have developed a CloudFormation template and web application to make scaling up Nextflow-based pipelines via deployment in AWS far simpler.  Using this approach:
 
1. a new pipeline can be deployed with  just a few clicks
2. end users can start and monitor the pipeline using a web page
3. limits can be enforced to prevent unexpected bills. 

By using AWS Lambda to provide the services that bind everything together, you only pay to store your data and for the EC2 compute used by batch when the pipeline is actually running.
 
We will talk about why and how we did this and ask for feedback on how our approach could be improved.

### Bio

Anthony Underwood is the bioinformatics implementation manager for the NIHR Global Health Research Unit in Genomic Surveillance of Antimicrobial Resistance (https://ghru.pathogensurveillance.net/), a project that aims to provide pathways for Low and Medium income countries to access genomic technologies for pathogen surveillance. His background is in bioinformatics for public health microbiology.

Ben Taylor is a senior software developer in the Centre for Genomic Pathogen Surveillance developing software that optimises common bioinformatics processes and delivers them through user-friendly interfaces. He’s previously worked for the UK Government and private companies to make it easier to use Cloud Infrastructure.

### Registration 

To attend Nextflow Camp 2019 register at [this link](https://www.crg.eu/en/event/coursescrg-nextflow-2019).