title=Robust and reproducible pipelines to support routine clinical diagnostic and research projects in oncology
date=2019-05-28
type=col8
tags=nextflow,nfcamp,2019,workshop
status=published
~~~~~~

## Robust and reproducible pipelines to support routine clinical diagnostic and research projects in oncology

### Matthias Monfort
*Research Engineer, Institut Curie, France* 

Institut Curie is a European Comprehensive Cancer Center which comprises a research center and three cancer hospitals. The bioinformatics platform support both research and patient care in particular by providing different kinds of analysis pipelines to address different use cases. 

On the one hand, diagnostic pipelines are used by hospital teams in their daily clinical routine to process patient samples. As the therapeutic decisions rely on their results,  the reproducibility and robustness of these pipelines are a major concern. 

On the other hand, research pipelines are more exploratory programs with frequent changes (tool versions, state-of-the-art for the domain) requiring more flexibility and a high level of customization (overriding default settings, adjusting some tools parameters) ; such pipelines have to be easily tunable, versatile, portable, and sometimes sharable with other stakeholders. 

To reach these objectives, the bioinformatics platform has made the choice of Nextflow as a common workflow manager for implementing these pipelines which raised technical and organisational challenges that we propose to expose and discussed. We will focus on how: 
- we defined  common guidelines both for Nextflow usage but also for more generic concerns (language independent, e.g. algorithmic) in order to migrate historical pipelines composed by thousand lines of legacy code, accumulated for years, involving different languages 
- we have turned the progressive deployment of Nextflow in our technical environment into an opportunity to generalize and improve our skills and practices regarding some well known best practices in many different areas such as containerisation, reproducibility, modularity and so on. 

To address some of the concerns mentioned above (portability, scalability, reproducibility, robustness, etc.) we have implemented a MPI (Message Passing Interface) version of the first steps of the variant calling NGS pipelines. Those steps are the alignment, the sorting and the marking of duplicates. We will provide a feedback experience on how we integrated programs optimized with MPI  within Nextflow and see how much effort it takes and how it could affect the performances.   

### Bio

Matthias validated the bioinformatics master of the University of Bordeaux by doing a 6 months internship in Eileen Furlong's laboratory at EMBL Heidelberg under the joint supervision of Charles Girardot. He then joined the Furlong team right after as a bioinformatician and software developer. Dynamix, his master project, is a genome browser plugin for automated genome browsing and was released in Bioinformatics (2017). 

Matthias joined Emmanuel Barillot's team, in the Cancer System Biology Unit at the Curie institute of Paris, in november 2018 and currently works for the Bio-IT/HPC platform under Philippe Hupé's supervision. The platform picked Nextflow to industrialise data processing and analysis workflows, used both for research purposes as well as assisted medical diagnosis.


### Registration 

To attend Nextflow Camp 2019 register at [this link](https://www.crg.eu/en/event/coursescrg-nextflow-2019).