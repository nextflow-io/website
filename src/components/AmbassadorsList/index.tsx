import React, { useState, useMemo } from "react";
import AmbassadorCard from "../AmbassadorCard/index.tsx";
import AmbassadorFilter from "../AmbassadorFilter/index.tsx";

interface Ambassador {
  name: string;
  img: string;
  country: string;
  github?: string;
  linkedin?: string;
  twitter?: string;
  mastodon?: string;
  bluesky?: string;
  title?: string;
  children?: string;
}

const ambassadors: Ambassador[] = [
  {
    name: "Abhinav Sharma",
    img: "abhinav.jpg",
    country: "za",
    github: "abhi18av",
    linkedin: "abhi18av",
    twitter: "abhi18av",
    title: "Nextflow Ambassador",
  },
  {
    name: "Adam Talbot",
    img: "adam_talbot.jpeg",
    country: "gb",
    github: "madamrtalbot",
    linkedin: "adam-talbot-9040a826",
    twitter: "adamrtalbot",
    title: "Nextflow Ambassador",
  },
  {
    name: "Adolf Oyesigye Mukama",
    img: "adolf-mukama.png",
    country: "ke",
    github: "adolfmukama",
    linkedin: "oyesigye-mukama-adolf-571b35296",
    twitter: "AdolfOyesigye",
    title: "Nextflow Ambassador",
    children:
      'Adolf is an Msc Bioinformatics Student at <a href="https://www.pu.ac.ke/" target="_blank"> Pwani University </a> and a research fellow at <a href="https://kemri-wellcome.org/" target="_blank"> Kemri-wellcome Trust Program</a>, Passionate about developing efficient and scalable pipelines for pathogen genomics.',
  },
  {
    name: "Abdoallah Sharaf",
    img: "abdoallah.png",
    country: "de",
    github: "abdo3a",
    linkedin: "abdoallah-sharaf-4a52614a",
    twitter: "abdo3a",
    title: "Nextflow Ambassador",
    children:
      'Abdoallah is a senior bioinformatician at SequAna - Sequencing Analysis Core Facility <a href="https://www.biologie.uni-konstanz.de/sequana/sequana/" target="_blank" > University of Konstanz, Germany</a>, Passionate about developing efficient and scalable pipelines to tackle complex bioinformatics challenges.',
  },
  {
    name: "Agyekum Richard",
    img: "ragyekum.jpeg",
    country: "gh",
    github: "QuadjoLegend",
    linkedin: "richardagyekum",
    title: "Nextflow Ambassador",
    children:
      "Richard is a Senior Medical Laboratory Scientist who has developed keen interest in Bioinformatics and Computational Biology. Richard currently mentors interns at HackBio.",
  },
  {
    name: "Alex Valcourt Caron",
    img: "alexvc.jpeg",
    country: "ca",
    github: "AlexVCaron",
    linkedin: "alex-valcourt-caron-545bb72b7",
    twitter: "AlexVCaron",
    title: "Nextflow Ambassador",
    children:
      'Alex Valcourt Caron is a PhD Candidate at the <a href="https://scil.usherbrooke.ca" target="_blank" >Sherbrooke Connectivity Imaging Laboratory</a >. He develops nextflow components for Magnetic Resonance Imaging (MRI) and is co-leading the <a href="https://nf-co.re/special-interest-groups/neuroimaging" target="_blank" > Neuroimaging Special Interest Group.</a>',
  },
  {
    name: "Ali Nawaz",
    img: "AliNawaz.jpeg",
    country: "pk",
    github: "AlleyNawaz",
    linkedin: "AlleyNawaz",
    twitter: "AlleyNawaz",
    title: "Nextflow Ambassador",
    children:
      'Ali Nawaz is a GeoInformatics Engineering student at <a href="https://nust.edu.pk" target="_blank">NUST</a>. He is the Community Organizer of <a href="https://linktr.ee/TFUGIslamabad" target="_blank">TensorFlow User Group Islamabad</a> and as Community Outreach for <a href="https://gdg.community.dev/gdg-cloud-islamabad/" target="_blank">GDG Cloud Islamabad</a>.',
  },

  {
    name: "Amofa Justice Ohene",
    img: "justice.JPG",
    country: "gh",
    github: "Iamamofa",
    linkedin: "justice-ohene-amofa-349b44173",
    twitter: "I_am_amofa",
    children:
      'Justice is a Bioinformatician and a ML/AI engineer at the <a href="https://noguchi.ug.edu.gh" target="_blank" >Noguchi Memorial Institute For Medical Research</a >. He is currently on the <a href="https://pangens.org" target="_blank" >PANGenS Project.</a>',
  },
  {
    name: "Amrei Binzer-Panchal",
    img: "Amrei.png",
    country: "se",
    github: "amrei-bp",
    children:
      'Amrei is a Bioinformatician at the <a href="https://www.slubi.se" target="_blank" >Bioinformatics Infrastructure</a > of the <a href="https://www.slu.se/" target="_blank" >Swedish University for Agricultural Sciences</a >, where she supports researchers and teaches reproducible bioinformatics, including Nextflow.',
  },
  {
    name: "Anabella Trigila",
    img: "anabella.jpeg",
    country: "ar",
    github: "atrigila",
    linkedin: "anabellat",
  },
  {
    name: "Andre Fonseca",
    img: "oandrefonseca.jpg",
    country: "br",
    github: "oandrefonseca",
    linkedin: "oandrefonseca",
    children:
      "Andre is a Bioinformatician Scientist with solid expertise in cancer biology, tumor immunology, and machine learning applications. He has worked in multiple institutions, including the prestigious MD Anderson Cancer Center.",
  },
  {
    name: "Antoine Buetti-Dinh",
    img: "AntoineBuettiDinh.JPG",
    country: "ch",
    github: "antoine-buetti",
    linkedin: "antoine-buetti-dinh",
    children:
      "Antoine is a bioinformatician with experience spanning human genetics, cancer biology, and microbial ecology. He loves tackling biological questions through math and computational modeling.",
  },

  {
    name: "Arthur Gymer",
    img: "ArthurGymer.jpg",
    country: "au",
    github: "awgymer",
    linkedin: "awgymer",
    children:
      'Arthur is a bioinformatics developer with experience developing high-throughput pipelines for long read cancer data. He is also an active <a href="https://nf-co.re/" target="_blank" >nf-core</a> member.',
  },
  {
    name: "Ashley Dederich",
    img: "AshleyDederich.png",
    country: "us",
    github: "ashdederich",
    linkedin: "ashley-dederich-bioinformatician",
    children:
      'Ashley is a Scientific Consultant at the <a href="https://chpc.utah.edu/" target="_blank">CHPC.</a> She has developed image analysis and genomic assembly pipelines and now uses her expertise to consult researchers on their scientific computing requirements.',
  },
  {
    name: "Bulut Hamali",
    img: "BulutHamali.png",
    country: "us",
    github: "buluthamali",
    linkedin: "bulut-hamali",
    twitter: "BioinfUniverse",
    children:
      'Bulut Hamali is a Bioinformatician at the <a href="https://med.uc.edu/depart/cancer-biology/" target="_blank" >UC Cincinnati</a >, studying HER2-positive breast cancer mechanisms. He holds a PhD from the Medical University of Vienna.',
  },
  {
    name: "Charalampos (Harris) Lazaris",
    img: "harris.jpg",
    country: "us",
    github: "chlazaris",
    linkedin: "chlazaris",
    twitter: "chlazaris",
    children:
      'Charalampos is a Senior Data Scientist in Oncology at Novartis Biomedical Research in Cambridge, MA. He earned his PhD from <a href="https://www.nyu.edu/">NYU.</a> He specializes in gene dysregulation in cancer and how to exploit it therapeutically.',
  },
  {
    name: "Chiachun Chiu",
    img: "chiachun.jpg",
    country: "tw",
    github: "godkin1211",
    twitter: "ChiachunChiu",
    linkedin: "michael-nostalgie-57630a61",
    children:
      'Chiachun is a bioinformatician at the <a href="http://sub.chimei.org.tw/55480/index.php/english/english03/13-english" target="_blank" > Center for Precision Medicine of Chi Mei Medical Center</a >. He\'s also the community organizer of <a href="https://www.facebook.com/groups/446434039038963" target="_blank" >Taipei Bioinformatics Omnibus</a >.',
  },
  {
    name: "Clément Igiraneza",
    img: "ClementIgiraneza.jpg",
    country: "rw",
    github: "igiraclement",
    linkedin: "igiraneza-clement-56441897",
    twitter: "IGIRANEZACLEME1",
    children:
      "Clement Igiraneza is a molecular biologist with expertise in infectious disease genomics, and bioinformatics data analysis. He has worked at the Rwanda Biomedical Center in the Molecular and Genomic Unit and focuses on genomic surveillance, sequencing technologies, diagnostic innovation, and capacity building for malaria drug-resistance monitoring.",
      },
  {
    name: "Cris Tuñi",
    img: "cristuni.jpeg",
    country: "es-ct",
    github: "ctuni",
    linkedin: "cristina-tuñí-i-domínguez-75a053145",
    twitter: "c_tunyi",
    children:
      '<a href="https://ctuni.dev" target="_blank">Cris</a> is a bioinformatics scientist and Ph.D. candidate at <a href="https://www.flomics.com/" target="_blank" >Flomics Biotech S.L.</a > They develop Nextflow RNA-Seq analysis pipelines to advance the field of early cancer diagnostics.',
  },
  {
    name: "Daniel Lundin",
    img: "daniel.jpeg",
    country: "se",
    github: "erikrikarddaniel",
  },
  {
    name: "Davi Marcon",
    img: "davi.jpeg",
    country: "br",
    github: "Mxrcon",
    linkedin: "davi-marcon-2088a722b",
    twitter: "mxrcon_",
  },
  {
    name: "Edmund Miller",
    img: "edmundmiller.png",
    country: "us",
    github: "Emiller88",
    linkedin: "edmund-miller-01974a105",
    twitter: "E_Miller88",
    children:
      'Edmund Miller is a Ph.D. Candidate in the Functional Genomics Lab at the <a href="https://www.utdallas.edu" target="_blank" >University of Texas at Dallas</a >. He\'s a <a href="https://nf-co.re/governance#maintainers" target="_blank" >maintainer</a > in the <a href="https://nf-co.re" target="_blank">nf-core</a> project involved in an eclectic group of nf-core projects.',
  },
  {
    name: "Edoardo Giacopuzzi",
    img: "edoardogiacopuzzi.png",
    country: "it",
    github: "edg1983",
    linkedin: "edoardo-giacopuzzi-a1654823",
    twitter: "supergecko",
    bluesky: "https://bsky.app/profile/edg1983.bsky.social",
    children:
      'Edoardo is a Senior Bioinformatician at <a href="https://humantechnopole.it/en" target="_blank" >Human Technopole</a> research institute in Milan, Italy. He works in population genetics using Nextflow to bring order to genomic data.',
  },

  {
    name: "Evangelos Karatzas",
    img: "evangelos.jpg",
    country: "gr",
    github: "vagkaratzas",
    linkedin: "vagkaratzas",
    children:
      'Evangelos is a Research Fellow in <a href="https://www.ebi.ac.uk/" target="_blank" >EMBL-EBI</a >. He is currently developing Nextflow/nf-core pipelines to generate and annotate metagenomics derived protein families.',
  },
  {
    name: "Fadinda Shafira",
    img: "fadindas.jpg",
    country: "id",
    github: "fadindashafira",
    linkedin: "fadindashafira",
    children:
      'Fadinda is an AI Engineer in Bioinformatics at <a href="https://www.kalbe.co.id" target="_blank" >Kalbe Farma</a >. She is part of Kalbe Digital Lab, focusing on digital biology. She promotes bioinformatics by developing courses for <a href="https://kdu.kalbe.co.id" target="_blank" >Kalbe Digital University</a >.',
  },
  {
    name: "Felipe Almeida",
    img: "felipe-almeida.jpg",
    country: "dk",
    github: "fmalmeida",
    linkedin: "almeida-fm",
    twitter: "fmarquesalmeida",
    children:
      'Felipe is a bioinformatician at <a href="https://www.zs.com/careers/where-we-work/europe/copenhagen" target="_blank" >ZS</a >. He is involved in projects related to pipelines, promoting and facilitating the use of Nextflow with guidance and trainings, when fit. Also an active <a href="https://nf-co.re/" target="_blank" >nf-core</a> member.',
  },
  {
    name: "Firas Zemzem",
    img: "FirasZemzem.jpg",
    country: "tn",
    github: "Zemzemfiras1",
    linkedin: "firaszemzem",
    twitter: "ZemzemFiras",
    children:
      "Firas is PhD student at the laboratory of Cytogenetics, Molecular Genetics, and Reproductive Biology at CHU Farhat Hached Sousse, focuses on unraveling the genetic mechanisms responsible for genetic disorders.",
  },

  {
    name: "Florian Heyl",
    img: "heylf.jpg",
    country: "de",
    github: "heylf",
    linkedin: "florian-heyl",
    children:
      "Florian Heyl is a Bioinformatician at DKFZ, Germany. He is bridging between nf-core, GHGA, Galaxy, and the single cell community. He holds a PhD in Bioinformatics from Freiburg University, Germany.",
  },
  {
    name: "Francesco Lescai",
    img: "francesco_nxf.png",
    country: "it",
    github: "lescai",
    linkedin: "francescolescai",
    twitter: "tokybo",
    children:
      'Francesco Lescai leads the <a href="https://lescailab.unipv.it" target="_blank" >Computational Genomics Lab</a > at the <a href="https://dbb.dip.unipv.it/en" target="_blank" >University of Pavia</a >. He teaches bioinformatics and he\'s a developer of <a href="https://nf-co.re" target="_blank" >nf-core</a> pipelines.',
  },
  {
    name: "Francisco Martin Garcia",
    img: "FranciscoMGarcia.jpeg",
    country: "ar",
    github: "garciafranciscomartn",
    linkedin: "garciafranciscomartin",
    children:
      'Francisco Martin Garcia is a Bioinformatician at <a href="https://www.garrahan.gov.ar/" target="_blank" >Garrahan Paediatric Hospital</a >, working on genomic variants, structural bioinformatics, and multi-omics approaches. He promotes the adoption of Nextflow in clinical and research settings.',
  },
  {
    name: "Franz AKE",
    img: "Franz_ake.png",
    country: "fr",
    github: "Franzx7",
    linkedin: "franz-arnold-ake-3a657b11b",
    twitter: "franz_ake",
    children:
      'Franz AKE is a Bioinformatician specializing in single-cell transcriptomics. He completed his PhD in Bioinformatics at <a href="https://p-cmrc.cat/research/plass-group/" target="_blank" > IDIBELL Institute</a >. He advocates for workflow automation and the integration of Nextflow in omics studies.',
  },
  {
    name: "Franziska Bonath",
    img: "franziska.jpeg",
    country: "se",
    github: "FranBonath",
    linkedin: "franziska-bonath-a0827a7a",
    twitter: "FranBonath",
    children:
      'Fran got her PhD in Developmental Biology at the University of Manchester, then moved to Stockholm where she eventually landed at the NGI, one of the founding institutions of nf-core. She is a <a href="https://nf-co.re/governance#core-team" target="_blank" >core team</a > member of <a href="https://nf-co.re/" target="_blank" >nf-core.</a>',
  },
  {
    name: "Friederike Hanssen",
    img: "friederikehanssen.png",
    country: "de",
    github: "FriederikeHanssen",
    linkedin: "friederike-hanssen",
    twitter: "RikeHanssen",
  },
  {
    name: "Georgie Samaha",
    img: "georgie.jpg",
    country: "au",
    github: "georgiesamaha",
    linkedin: "georgie-samaha-95095b230",
    children:
      'Georgie leads the bioinformatics group at the <a href="https://www.sydney.edu.au/research/facilities/sydney-informatics-hub.html" target="_blank" >Sydney Informatics Hub</a >, University of Sydney. She is working toward making bioinformatics more accessible by developing public digital infrastructure with the <a href="https://www.biocommons.org.au/" target="_blank" >Australian BioCommons</a >.',
  },
  {
    name: "Gisela Pattarone",
    img: "GiselaPattarone.png",
    country: "ar",
    github: "gpattarone",
    linkedin: "giselapattarone",
    children:
      'Gisela is a Bioinformatician at <a href="http://www.cnea.gov.ar/" target="_blank" >CNEA</a >, where she designs and executes a meta-analysis project focused on gene expression profiles to identify diagnostic and therapeutic biomarkers in pediatric gliomas.',
  },
  {
    name: "Graeme Grimes",
    img: "ggrimes.jpg",
    country: "gb",
    github: "ggrimes",
    linkedin: "graeme-grimes-a753743a",
    twitter: "bioggrimes",
    children:
      'Graeme is a Bioinformatician at the <a href="https://www.ed.ac.uk/institute-genetics-cancer" target="_blank" >Institute of Genetics & Cancer</a > at the <a href="https://www.ed.ac.uk" target="_blank" >Univeristy of Edinburgh</a> in Scotland, he is also the lesson maintainer of the Carpentries <a href="https://carpentries-incubator.github.io/workflows-nextflow/" target="_blank" >Introduction to Bioinformatics workflows with Nextflow and nf-core</a >.',
  },
  {
    name: "Hemanoel Passarelli",
    img: "hemanoel-passarelli.jpg",
    country: "br",
    github: "Passarelli-bio",
    linkedin: "hemanoel-passarelli",
    twitter: "he_passarelli",
    children:
      'Hemanoel is Bioinformatics Engineer at <a href="https://gen-t.science/" target="_blank" >Gen-t</a >, where he develops Nextflow pipelines to analyze human genetics data, enhancing Precision Medicine efforts in Brazil.',
  },
  {
    name: "Houcemeddine Othman",
    img: "houcem.png",
    country: "tn",
    github: "hothman",
    linkedin: "houcemeddine-othman-5502b541",
    twitter: "Houcemeddi61361",
    children:
      "Houcemeddine Othman is an Assistant Professor of Bioinformatics from Tunisia, working on the development of workflows for clinical genomics.",
  },
  {
    name: "Hyun-Hwan Jeong",
    img: "hjeong.jpg",
    country: "us",
    github: "hyunhwan-bcm",
    linkedin: "hyunhwan-jeong",
    twitter: "hyunhwan_jeong",
    children:
      'Hyun-Hwan Jeong is an Assistant Professor of Pediatrics at Baylor College of Medicine. He is also an open-source contributor to multiple projects, including <a href="https://github.com/hyunhwan-jeong/CB2" target="_blank" >CB<sup>2</sup></a> and <a href="https://github.com/LiuzLab/AI_MARRVEL" target="_blank" >AI-MARRVEL</a >.',
  },
  {
    name: "Ira Iosub",
    img: "iraiosub.jpeg",
    country: "gb",
    github: "iraiosub",
    linkedin: "ira-iosub-618254276",
    twitter: "IraIosub",
  },
  {
    name: "Isha Parikh",
    img: "ishaparikh.png",
    country: "us",
    github: "isha2106",
    linkedin: "isha2106",
  },
  {
    name: "Jacques Dainat",
    img: "jacques-dainat.jpg",
    country: "fr",
    github: "Juke34",
    linkedin: "jacques-dainat-02257376",
    mastodon: "https://genomic.social/@jacquesdainat",
    bluesky: "https://bsky.app/profile/jacquesdainat.bsky.social",
    children:
      'Currently Bioinformatician at <a href="https://en.ird.fr" target="_blank" > IRD </a> and part of the <a href="https://bioinfo.ird.fr" target="_blank" > i-Trop platform</a >, Jacques enjoys simplifying complex analyses through automation, minimizing technical barriers and enabling biologists to make the most of the power of bioinformatics.',
  },
  {
    name: "James Fellows Yates",
    img: "james.jpeg",
    country: "de",
    github: "jfy133",
    linkedin: "james-fellows-yates-999859179",
    twitter: "jfy133",
    children:
      'James has a PhD in Bioinformatics from <a href="https://www.uni-jena.de" target="_blank" >FSU Jena</a> and <a href="https://www.eva.mpg.de" target="_blank" >MPI-EVA</a >. He is currently a PostDoc at <a href="https://www.leibniz-hki.de/en/" target="_blank" >Leibniz-HKI</a >, working in metagenomics with a specialism in ancient DNA. He is also a <a href="https://nf-co.re/governance#core-team" target="_blank" >core team</a > member of <a href="https://nf-co.re/" target="_blank" >nf-core</a> <a href="https://nf-co.re/governance#core-team" target="_blank" >core team.</a>',
  },
  {
    name: "Javier Carpinteyro Ponce",
    img: "javi.jpg",
    country: "mx",
    github: "javibio-git",
    linkedin: "javier-carpinteyro-ponce",
    children:
      'Javier (Javi) is a bioinformatics research associate at <a href="https://carnegiescience.edu" target="_blank" >Carnegie Science</a >. He currently provides bioinformatics support to bench-based researchers helping them with their data analysis needs. He is also interested in developing reproducible and scalable workflows for population genomics.',
  },
  {
    name: "Jeferyd Yepes",
    img: "jeferyd.jpg",
    country: "ch",
    github: "jeffe107",
    linkedin: "jeferyd-yepes-garcía",
    twitter: "JeferydY",
    children:
      '<a href="https://www.jeferydyepes.com" target="_blank" >Jeferyd</a > is currently a PhDc in bioinformatics at the University of Fribourg and SIB. He is the developer of <a href="https://f1000research.com/articles/13-640" target="_blank" >MAGFlow/BigMAG</a >, a Nextflow-based tool to compare metagenomics pipelines.',
  },
  {
    name: "Jehee Lee",
    img: "jehee.png",
    country: "kr",
    github: "jhlee0637",
    linkedin: "jehee-lee-202002",
    children:
      "Jehee is a bioinformatician from Seoul, South Korea. He aims to become an expert in harmonizing AI, BI, and Cloud technologies to accelerate biological research.",
  },
  {
    name: "Jelena Pejovic Simeunovic",
    img: "jelenap.jpg",
    country: "rs",
    github: "JelPej",
    linkedin: "jelena-pejovic-simeunovic-ph-d-bb300468",
    children:
      "Jelena is an Advanced Bioinformatics Engineer at Loka, where she builds bioinformatics workflows and ensures seamless integration and optimal performance using Nextflow. She holds a PhD in Nano Biotechnology.",
  },
  {
    name: "Jennifer (Jen) Reeve",
    img: "jenniferreeve.jpg",
    country: "nz",
    github: "jen-reeve",
    linkedin: "jennifer-reeve-microbio",
    children:
      "Jen is a Software Developer with NetValue Ltd in Hamilton, New Zealand. Her current work uses Nextflow in biomedical contexts, but her scientific background spans microbiology, biochemistry and geochemistry.",
  },
  {
    name: "Jimmy Trace Lail",
    img: "jimmylail.jpg",
    country: "us",
    github: "tracelail",
    linkedin: "tracelail",
    children:
      "Trace is a graduate student in the Master of Bioinformatics program at Northeastern University. He brings industry experience in gene therapy and cell culture, with interests in conservation and environmental biotechnology applications.",
  },
  {
    name: "John Vusich",
    img: "johnvusich.png",
    country: "us",
    github: "johnvusich",
    linkedin: "vusich",
    twitter: "johnvusich",
    children:
      "John is a PhD candidate in the Andrechek Lab at Michigan State University. His research focuses on integrating multi-omics to characterize the genomic and epigenomic landscape of breast cancer.",
  },
  {
    name: "Jose Espinosa-Carrasco",
    img: "joseespinosacarrasco.png",
    country: "es",
    github: "JoseEspinosa",
    linkedin: "jose-espinosa-carrasco",
    children:
      'Jose is a Bioinformatician at the CBCRG group at the <a href="https://www.crg.eu" target="_blank" >CRG</a >, the laboratory where Nextflow was born. As an early adopter of Nextflow, he is also a member of the <a href="https://nf-co.re/" target="_blank" >nf-core</a > <a href="https://nf-co.re/governance#core-team" target="_blank" >core team.</a>',
  },
  {
    name: "Janan Gawra",
    img: "Janan.jpeg",
    country: "es",
    github: "JGawra",
    linkedin: "janangawra",
    children:
      'Dr. Janan Gawra is a postdoctoral researcher at the <a href="https://www.idaea.csic.es/person/janan-gawra/" target="_blank" >IDAEA-CSIC</a >. He is part of the <a href="https://epiboost.web.ua.pt/?page_id=6445" target="_blank" >EPIBOOST</a> project, studying the epigenetic and gene expression changes in zebrafish and seabass in response to environmental pollutants.',
  },
  {
    name: "João Vitor Cavalcante",
    img: "joaocavalcante.png",
    country: "br",
    github: "jvfe",
    linkedin: "joao-vitor-cavalcante",
    children:
      'João is a MSc student affiliated with the <a href="https://dalmolingroup.imd.ufrn.br/" target="_blank" >Dalmolin Systems Biology Group</a >, at the Bioinformatics Multidisciplinary Environment in Natal, Brazil. His interests are in metagenomics, neurogenetics and scientific software development.',
  },

  {
    name: "Julia Apolonio de Amorim",
    img: "juliaapolonio.jpeg",
    country: "br",
    github: "juliaapolonio",
    linkedin: "juliaapolonio",
    children:
      'Julia is a MSc student at the <a href="https://www.ufrn.br/en" target="_blank" >Rio Grande do Norte Federal University</a > in Natal, Brazil. She is interested in GWAS and post-GWAS analyses, functional genomics and genetics of psychiatric and neurodegenerative diseases.',
  },
  {
    name: "Júlia Mir Pedrol",
    img: "juliamirpedrol.png",
    country: "es-ct",
    github: "mirpedrol",
    linkedin: "juliamirpedrol",
    twitter: "juliamirpedrol",
    children:
      'Júlia is a bioinformatician at the Computational Biology and Health Genomics group from the CRG. She is also an nf-core tools developer and a member of the <a href="https://nf-co.re/" target="_blank" >nf-core</a > <a href="https://nf-co.re/governance#core-team" target="_blank" >core team.</a>',
  },

  {
    name: "Kimberly Christine Coetzer",
    img: "KCCoetzer_NF_ambassador.png",
    country: "za",
    github: "Kimmiecc19",
    linkedin: "kimberly-christine-coetzer-hugo-a02049151",
    children:
      "Kimberly is a PhD candidate in Bioinformatics and Computational Biology at Stellenbosch University, South Africa. She has experience in clinical proteomics and human genetics, specifically in the field of rare diseases.",
  },
  {
    name: "Kim Huy Vo",
    img: "kimhuyvo.png",
    country: "vn",
    github: "vkhuy",
    linkedin: "kimhuyvo", 
    children:
      'Kim Huy is a bioinformatician at <a href="https://www.ktest.vn/" target="_blank">KTest Science Co. Ltd</a>, with hands-on experience in developing robust analysis pipelines for both long-read and short-read sequencing data.',
  },
  {
    name: "Kinley Tenzin",
    img: "kinleytenzin.jpg",
    country: "bt",
    github: "tkinley",
    title: "Nextflow Ambassador",
    children:
    "Kinley Tenzin is a PhD student at Kansas State University working in microbial bioinformatics, with interests in genome assembly, mobile genetic elements, and scalable workflow development with Nextflow.",
  },
  
  {
    name: "Kobe Lavaerts",
    img: "kobelavaerts.png",
    country: "be",
    github: "kobelavaerts",
    linkedin: "kobe-lavaerts-170489191",
    children:
      'Kobe Lavaerts is a staff bioinformatician at the genomics core facility <a href="https://nucleomicscore.sites.vib.be/en" target="_blank" >Nucleomics Core</a > in <a href="https://vib.be/en" target="_blank" >VIB</a >, Belgium. He develops and maintains the Nextflow pipelines used in the core facility. He also teaches the <a href="https://github.com/vib-tcp/nextflow-workshop" target="_blank" >VIB Nextflow training</a > twice a year.',
  },
  {
    name: "Kristina K. Gagalova",
    img: "KristinaGagalova.jpg",
    country: "au",
    linkedin: "kristina-gagalova",
    github: "KristinaGagalova",
    children:
      'Kristina is a Bioinformatics Scientist at the <a href="https://www.curtin.edu.au/" target="_blank" >Curtin University</a >/<a href="https://www.ccdm.com.au" target="_blank" >Centre for Crop and Disease Management</a >. She develops innovative algorithms and pipelines for pangenomics & structural bioinformatics, advancing crop disease research.',
  },
  {
    name: "Kübra Narcı",
    img: "kubranarci2.jpg",
    country: "de",
    github: "kubranarci",
    mastodon: "https://www.ghga.de/about-us/team-members/narci-kuebra",
    linkedin: "kubranarci",
    twitter: "kubranarci",
    children:
      "Kübra Narcı is a bioinformatician at DKFZ, Germany. She is bridging between nf-core and GHGA. She holds a PhD in Bioinformatics from METU, Turkey.",
  },

  {
    name: "Lara Ianov",
    img: "Ianov.jpg",
    country: "us",
    github: "lianov",
    linkedin: "lara-ianov",
    children:
      'Lara Ianov is the Co-Director of the <a href="https://www.uab.edu/cores/ircp/bds">UAB Biological Data Science Core</a >. She directs the development of transcriptomics pipelines including the nf-core/scnanoseq pipeline.',
  },
  {
    name: "Lili Andersson-Li",
    img: "lilianderssonli.jpg",
    country: "se",
    github: "LilyAnderssonLee",
    mastodon:
      "https://ki.se/forskning/forskningsomraden-centrum-och-natverk/forskargrupper/evolution-och-epidemiologi-for-hiv-och-enterovirus-jan-albert-grupp",
    linkedin: "lili-andersson-li-78565b134",
    children:
      "Lili Andersson-Li is a bioinformatician at KI, Sweden. She develops pipelines for metagenomics, viral typing and viral drug resistance testing. She holds a PhD in population genetics.",
  },
  {
    name: "Lynn Langit",
    img: "lynnlangit.png",
    country: "us",
    github: "lynnlangit",
    linkedin: "lynnlangit",
    children:
      'Lynn is an independent cloud AI architect building tools for bioinformatics teams world-wide. She is also part of <a href="https://dimi-lab.github.io/team/" target="_blank" >DIMI Lab</a >.',
  },
  {
    name: "Louis Le Nézet",
    img: "louislenezet.png",
    country: "fr",
    github: "louislenezet",
    linkedin: "louis-le-nezet",
    children:
      'Louis is a PhD student in bioinformatics at the <a href="https://igdr.univ-rennes.fr/en" target="_blank" >IGDR</a >. His research focuses on developing and applying computational methods to study the genetics basis of the canine hip dysplasia.',
  },
  {
    name: "Luca Cozzuto",
    img: "lcozzuto.jpeg",
    country: "es",
    github: "lucacozzuto",
    linkedin: "cozzuto",
    twitter: "lucacozzuto",
    children:
      'Luca is a Senior Bioinformatician at <a href="https://biocore.crg.eu/" target="_blank" > the Bioinformatics technology platform </a> of <a href="https://www.crg.eu/" target="_blank" > the Centre for Genomic Regulation.</a > He has a PhD in Molecular Medicine from <a href="https://www.semm.it/" target="_blank" > The European School of Molecular Medicine.</a>',
  },
  {
    name: "Luke Pembleton",
    img: "luke.jpg",
    country: "au",
    github: "lpembleton",
    mastodon: "https://genomic.social/@lwpembleton",
    linkedin: "lpembleton",
    twitter: "lwpembleton",
    children:
      'Luke has a PhD in Molecular Genetics from <a href="https://www.latrobe.edu.au/" target="_blank" >La Trobe University</a > and <a href="https://agriculture.vic.gov.au/about/our-research" target="_blank" >AVR</a >, Australia. He is currently a Genomic Breeding Scientist and Global Strategic Science Manager at <a href="https://www.barenbrug.com" target="_blank" >Barenbrug</a >.',
  },
  {
    name: "Mahesh Binzer-Panchal",
    img: "mahesh_binzer-panchal.png",
    country: "se",
    github: "mahesh-panchal",
    mastodon: "https://genomic.social/@mbp",
    linkedin: "mahesh-binzer-panchal-79a726a2",
    twitter: "arcane_mahesh",
    children:
      'Mahesh is a bioinformatician for <a href="https://www.nbis.se/" target="_blank" >NBIS</a >, <a href="https://www.scilifelab.se/" target="_blank" >SciLifeLab</a >, the Swedish node of <a href="https://elixir-europe.org/" target="_blank" >Elixir</a >. He supports research groups, performing de novo genome assembly, or workflow development. He is also an nf-core maintainer.',
  },
  {
    name: "Mahima Sanjay Gomladu",
    img: "mahima.jpeg",
    country: "no",
    github: "mohe1linux",
    linkedin: "mahima-sanjay-gomladu-ms-a52867176",
    twitter: "GomladuMahima",
    children:
      '<a href="https://www4.uib.no/en/find-employees/Mahima.Sanjay.Gomladu%2C.MS" target="_blank" >Mahima Sanjay Gomladu</a >, is a Senior Bioinformatics Engineer at the <a href="https://www.uib.no/en/clin2/genomics" target="_blank" >Genomics Core Facility</a >, University of Bergen, Norway. She designs modular, scalable pipelines using Nextflow to drive large-scale sequencing analyses (WES, WGS, RNA-seq, single-cell, spatial omics). As an ELIXIR collaborator, she promotes reproducible sciencand FAIR data principles.',
  },
  {
    name: "Margherita Mutarelli",
    img: "margherita.png",
    country: "it",
    github: "daisymut",
    linkedin: "margherita-mutarelli",
    children:
      'Margherita is a researcher at the <a href="https://www.isasi.cnr.it/" target="_blank" > Institute of Applied Sciences and Intelligent Systems (CNR-ISASI)</a > with a Ph.D. in Computational Biology and long term experience in NGS analysis. Happy member of nf-core.',
  },
  {
    name: "Marie Lataretu",
    img: "marie_lataretu.jpg",
    country: "de",
    github: "MarieLataretu",
    linkedin: "marie-lataretu-b55103218",
    children:
      'Marie is a bioinformatician at the <a href="https://www.rki.de/EN/Home/homepage_node.html" target="_blank" >Robert Koch Institute</a >. She mainly develops pipelines for NGS data and viroinformatics.',
  },
  {
    name: "Martin Beracochea",
    img: "martin_beracochea.jpg",
    country: "gb",
    twitter: "_martinbc",
    github: "mberacochea",
    linkedin: "martin-beracochea-87802b44",
    children:
      'Martin is the Production lead of <a href="https://www.ebi.ac.uk/metagenomics/" target="_blank" >MGnify</a> at <a href="https://www.ebi.ac.uk/" target="_blank" >EMBL-EBI</a >. He oversees the development and maintenance of web services and bioinformatics pipelines.',
  },
  {
    name: "Mateus Falco",
    img: "mateus.jpeg",
    country: "br",
    github: "mlfalco-bioinfo",
    linkedin: "mateusfalco",
    children:
      'Mateus is a PhD student on Pharmaceutical Sciences at <a href="https://www3.unicentro.br/ppgcf/" target="_blank" >Unicentro</a >, working at the facility <a href="https://ipec.org.br/" target="_blank" >IPEC</a> in wetlab and bioinformatics on Program Genomas Paraná.',
  },
  {
    name: "Matthias De Smet",
    img: "matthias_de_smet.jpeg",
    country: "be",
    github: "matthdsm",
    linkedin: "desmetmatthias",
  },
  {
    name: "Matthias Hörtenhuber",
    img: "matthias.jpeg",
    country: "se",
    github: "mashehu",
    twitter: "mashehu",
  },
  {
    name: "Maxime Garcia",
    img: "maxime.jpeg",
    country: "se",
    github: "maxulysse",
    linkedin: "maxugarcia",
    twitter: "gau",
    children:
      'Maxime, bioinfomagician, develops nf-core pipelines, mainly Sarek, and is a member of the <a href="https://nf-co.re/" target="_blank" >nf-core</a > <a href="https://nf-co.re/governance#core-team" target="_blank" >core team.</a>',
  },
  {
    name: "Michael Heuer",
    img: "michael_heuer.png",
    country: "us",
    github: "heuermh",
    twitter: "michael_l_heuer",
    children:
      'Michael Heuer (he/him/his) is a Staff Engineer at <a href="https://mammoth.bio/" target="_blank" >Mammoth Biosciences</a >. Michael contributes to several nf-core workflows and is a member of the nf-core <a href="https://nf-co.re/governance#safety" target="_blank" >safety team.</a>',
  },

  {
    name: "Muntadher Jihad",
    img: "munji.png",
    country: "gb",
    github: "muntajihad",
    linkedin: "muntadherjihad",
    twitter: "munta04",
    children:
      'Muntadher is a Bioinformatician at <a href="https://www.cruk.cam.ac.uk/" target="_blank" >CRUK CI</a >, develops statistical models to identify cellular crosstalk in pancreatic cancer and develops Nextflow pipelines for scRNA-seq variant calling.',
  },
  {
    name: "Nicholas Owen",
    img: "nowen.png",
    country: "gb",
    github: "nicholas-owen",
    linkedin: "nicholas-owen-3ab03820",
    twitter: "dr__no",
    children:
      'Nick is a Principal Research Data Steward / Bioinformatician at <a href="https://www.ucl.ac.uk/" target="_blank" >UCL Advanced Research Computing</a >. He has a strong background in splicing and disease, developing pipelines for NGS analysis using Nextflow.',
  },

  {
    name: "Nicolas Vannieuwkerke",
    img: "nicolas.jpeg",
    country: "be",
    github: "nvnieuwk",
    linkedin: "nicolas-vannieuwkerke-316874163",
    children:
      'Nicolas Vannieuwkerke is a Bioinformatician at <a href="https://www.cmgg.be/" target="_blank" >CMGG</a >. He mainly develops Nextflow pipelines for in-house deployment using the nf-core standards. He also is an avid believer in "Learning on the job" and is always eager to learn',
  },

  {
    name: "Nicolas Servant",
    img: "nservant.jpg",
    country: "fr",
    github: "nservant",
    linkedin: "nicolas-servant-0577471",
    children:
      'Nicolas is the co-director of the <a href="https://institut-curie.org/plateforme/curiecoretech-bioinformatics-cubic" target="_blank" > bioinformatics core facility of Institut Curie</a >. His main interests are cancer biology/epigenetics/AI. Member of the <a href="https://nf-co.re/" target="_blank" >nf-core</a > community, he promotes Nextflow for good programming practices.',
  },
  {
    name: "Nishat Fatima Mohammad",
    img: "nishat.jpeg",
    country: "us",
    github: "NishatMohammad",
    linkedin: "nishat-fatima-mohammad",
    children:
      "Nishat is a Medical Doctor and Bioinformatician using ML/AI techniques to find solutions to Biomedical problems. She is passionate about scalable and reproducible workflows that streamline pipelines for the future of science.",
  },
  {
    name: "Nour Mahfel",
    img: "nour.jpeg",
    country: "gb",
    github: "nourmahfel",
    linkedin: "nour-mahfel-010568182",
    children:
      "Nour Mahfel is a trainee clinical scientist on the Scientist Training Programme, specialising in genetic diagnostics at Birmingham Women's Hospital NHS.",
  },
  {
    name: "Omer Ali",
    img: "omer.png",
    country: "nr",
    github: "Omer0191",
    linkedin: "omer0191",
    children:
      "Omer Ali is a Bioinformatician at Oslo University Hospital, Oslo Norway. He works with TSO500 and whole genome data. Part of the team who provides advance molecular diagnostics to cancer patients.",
  },

  {
    name: "Phil Ewels",
    img: "phil.jpg",
    country: "se",
    github: "ewels",
    linkedin: "philewels",
    twitter: "tallphil",
    children:
      'Phil Ewels is Senior Product Manager for Open Source Software at <a href="https://www.seqera.io" target="_blank" >Seqera</a >, working on strengthening the Nextflow and nf-core communities. He holds a PhD in Molecular Biology from the <a href="https://www.cam.ac.uk" target="_blank" >University of Cambridge</a >, UK. He co-founded <a href="https://nf-co.re/" target="_blank" >nf-core</a> and is a member of the <a href="https://nf-co.re/governance#core-team" target="_blank" >core team.</a>',
  },
  {
    name: "Pritam Kumar Panda",
    img: "pritam_photo.png",
    country: "us",
    github: "pritampanda15",
    linkedin: "pritam-kumar-panda",
    twitter: "pritamkpanda",
    children:
      'Pritam is a Postdoctoral scholar in the <a href="https://profiles.stanford.edu/pritam-panda"> Department of Anesthesiology, Perioperative and Pain Medicine</a> at <a href="https://med.stanford.edu/profiles/pritam-panda"> Stanford University, School of Medicine</a>, California, designing novel anesthetics suitable for battlefield conditions.',
  },
  {
    name: "Rachel Griffard-Smith",
    img: "rachelgriffardsmith.jpg",
    country: "us",
    github: "rachelgriffard",
    linkedin: "rachelgriffard",
    children:
      'Rachel is a bioinformatician in the <a href="https://www.kumc.edu/school-of-medicine/academics/departments/biostatistics-and-data-science.html"> Department of Biostatistics & Data Science</a> at the University of Kansas Medical Center in Kansas City, Kansas analyzing and developing tools for microbiome data and other next generation sequencing data.',
  },
  {
    name: "Ramiro Barrantes Reynolds",
    img: "ramiro.jpg",
    country: "us",
    github: "ramirobarrantes",
    linkedin: "ramiro-barrantes-reynolds",
    children:
      'Bioinformatician/statistician/data scientist for the <a href="https://www.med.uvm.edu/vigr/home">Genomics/Data Science Core</a > and the <a href="http://www.med.uvm.edu/tgircobre/home">Infectious Disease Center</a> at the <a href="https://www.uvm.edu/"> University of Vermont</a >, USA. Love Nextflow!!',
  },
  {
    name: "Raquel Manzano",
    img: "raquelmanzano.jpg",
    country: "gb",
    github: "RaqManzano",
    linkedin: "raquel-manzano-bioinformatician",
    twitter: "appletreewoman",
    children:
      'Raquel is a senior bioinformatician and final year PhD student at <a href="https://www.cruk.cam.ac.uk/" target="_blank" >CRUK Cambridge Institute</a >. She is also a trainer at the <a href="https://bioinfotraining.bio.cam.ac.uk" target="_blank" >University of Cambridge Bioinformatics Training Unit</a >. For her PhD she developed <a href="https://nf-co.re/rnadnavar/dev" target="_blank" >#rnadnavar.</a>',
  },
  {
    name: "Rayan Hassaïne",
    img: "RayanHassaine.jpg",
    country: "nl",
    github: "rhassaine",
    linkedin: "rayan-hassaine",
    children:
      'Rayan is a Bioinformatician at <a href="https://www.hartwigmedicalfoundation.nl/en/">Hartwig Medical Foundation</a > in Amsterdam. He is one of the developpers of WiGiTS & among the maintainers of <a href="https://github.com/nf-core/oncoanalyser">nf-core/oncoanalyser</a> pipeline.',
  },
  {
    name: "Riley Grindle",
    img: "rgrindle.jpg",
    country: "us",
    github: "Riley-Grindle",
    linkedin: "riley-grindle",
    children:
      'Riley is a Bioinformatician at <a href="https://mdibl.org/" target="_blank" >MDIBL</a >. He has helped develop pipelines for his institutions in-house use and is an active member of the <a href="https://nf-co.re/" target="_blank" >nf-core</a> community.',
  },

  {
    name: "Saba Nafees",
    img: "saba.jpeg",
    country: "us",
    github: "snafees",
    linkedin: "saba-nafees",
    twitter: "saba_nafees314",
  },
  {
    name: "Sameesh Kher",
    img: "SameeshKher.jpg",
    country: "de",
    github: "khersameesh24",
    linkedin: "khersameesh24",
    twitter: "khersameesh24",
    children:
      'Sameesh Kher is a Bioinformatician at <a href="https://www.ghga.de/">GHGA</a> in Heidelberg. He is the developer of the <a href="https://github.com/nf-core/spatialxe">nf-core/spatialxe</a> pipeline.',
  },
  {
    name: "Samuel Ruiz-Pérez",
    img: "samuelruizperez.jpg",
    country: "mx",
    github: "samuelruizperez",
    linkedin: "samuelruizperez",
    twitter: "samuelruizperez",
    mastodon: "https://genomic.social/@samuelruizperez",
    bluesky: "https://bsky.app/profile/samuelruizperez.bsky.social",
    children:
      'Sam is a MSc in Bioinformatics student at the University of Copenhagen and a bioinformatician at the <a href="https://www.cancer.dk/danish-cancer-institute/research-groups/epigenome-replication-and-maintenance/" target="_blank"> Center for Epigenetic Cell Memory (EpiC)</a>, Danish Cancer Institute. He works with nascent DNA sequencing data and develops pipelines and machine learning models to study replication and epigenome maintenance.',
  },
  {
    name: "Sanzida Akhter",
    img: "sanzida.jpg",
    country: "bd",
    github: "sanzidaanee",
    linkedin: "sanzida-akhter-anee-47817b288",
    children:
      "Sanzida is a CSE and Bioinformatics master's student and is working on developing cancer biomarkers with machine learning algorithms.",
  },
  {
    name: "Sarai Varona",
    img: "sarai_varona.png",
    country: "es-pv",
    github: "svarona",
    linkedin: "sarai-varona-fernández-abb51013a",
    children:
      'Sarai Varona is a bioinformatician from the Basque Country, currently based in Madrid, where she works at the <a href="https://www.isciii.es/ub/unidad" target="_blank">Bioinformatics Unit of the Carlos III Health Institute</a>. She is the current maintainer of the <a href="https://nf-co.re/viralrecon/" target="_blank">nf-core/viralrecon</a> pipeline.'
  },
  {
    name: "Shiyi Yin",
    img: "shiyi.jpeg",
    country: "us",
    github: "yinshiyi",
    linkedin: "shiyi-yin",
    children:
      "Shiyi Yin is a genomics scientist in the SF bay area, working in compbio using nextflow for biotech R&D.",
  },
  {
    name: "Sebastián González Tirado",
    img: "SebastianGT.jpeg",
    country: "de",
    github: "sebgoti",
    linkedin: "sebastián-gonzález-tirado-77b7981a7",
    children:
      'Sebastián is doing his PhD at the <a href="https://www.medizinische-fakultaet-hd.uni-heidelberg.de/einrichtungen/institute/institute-for-computational-biomedicine" target="_blank" >Institute for Computational Biomedicine</a > in Heidelberg. He is interested on the democratization of bioimage analysis workflows and spatial omics.',
  },
  {
    name: "Simon Murray",
    img: "simon_murray.jpg",
    country: "gb",
    github: "SimonDMurray",
    linkedin: "simon-murray-856102156",
    children:
      'Simon Murray is a Bioinformatics Engineer at <a href="https://www.genomicsengland.co.uk" target="_blank" >Genomics England</a > in Research Data Products where he develops and maintains pipelines for the <a href="https://www.genomicsengland.co.uk/research/research-environment" target="_blank" >Research Environment</a >.',
  },
  {
    name: "Simon Pearce",
    img: "simon_pearce.jpeg",
    country: "gb",
    github: "SPPearce",
    linkedin: "simonppearce",
    children:
      'Simon Pearce is a Principal Bioinformatician at the <a href="https://www.cruk.manchester.ac.uk/Our-Research/Cancer-Biomarker-Centre/CEP-Home/">Cancer Research UK National Biomarker Centre</a >, where he works on detecting and subtyping cancer using circulating cell-free DNA methylation.',
  },
  {
    name: "Sofia Stamouli",
    img: "sofia.jpeg",
    country: "se",
    github: "sofstam",
    linkedin: "sofia-stamouli-a22b5477",
  },
  {
    name: "Susan Collins",
    img: "susancollins.jpg",
    country: "us",
    github: "susancollins",
    linkedin: "susan--collins",
    children:
      'Susan Collins is a Bioinformatics Engineer based in Vermont, building scalable, production-grade genomics workflows with Nextflow.',
  },
  {
    name: "Thai-Huy Tran",
    img: "thai-huytran.jpg",
    country: "vn",
    github: "huymtraan",
    linkedin: "thaihuy-tran",
    children:
      'Thai-Huy Tran is a final-year student at International University (VNU-HCM), Vietnam, with experience developing cfDNA data analysis pipelines and performing single-cell RNA analysis for cancer biology.',
  },
  {
    name: "Ziad Al-Bkhetan",
    img: "Ziad.jpg",
    country: "au",
    github: "ziadbkh",
    linkedin: "ziad-al-bkhetan",
    children:
      '<a href="https://www.biocommons.org.au/lb-ziad">Ziad</a> is the Product Manager, Bioinformatics Platforms at <a href="http://biocommons.org.au/">Australian BioCommons</a >. He leads several bioinformatics focused services such as <a href="https://www.biocommons.org.au/seqera-platform">the Australian Nextflow Seqera Service</a > and <a href="https://www.biocommons.org.au/ables">ABLeS.</a>',
  },
  {
    name: "Ze Yu",
    img: "zeyu.png",
    country: "us",
    github: "EZUY",
    linkedin: "ze-yu-9024b215b",
    children:
      'Ze Yu is a Computational biologist at UT Southwestern <a href="https://labs.utsouthwestern.edu/bioinformatics-lab" target="_blank">Bioinformatics Lab</a>, focusing on scalable genomics workflows and production-grade Nextflow DSL2 pipelines.',
  },
  {
    name: "Zohaib Anwar",
    img: "zohaibanwar.jpeg",
    country: "ca",
    github: "anwarMZ",
    linkedin: "mzohaibanwar",
    twitter: "zohaibanwar_",
    children:
      'Zohaib is a Senior Postdoc Researcher at the <a href="https://cidgoh.ca/" target="_blank" >Center for Infectious Disease Genomics and One Health</a >. He\'s the developer of <a href="https://virusmvp.org" target="_blank" >VIRUS-MVP</a >, powered by a genomics workflow developed in Nextflow with nf-core modules.',
  },
  // Ambassador Program Staff
  {
    name: "Marcel Ribeiro-Dantas",
    img: "marcel.jpg",
    country: "br",
    github: "mribeirodantas",
    linkedin: "mribeirodantas",
    twitter: "mribeirodantas",
    title: "Ambassador Program Staff",
    children:
      'Marcel holds a PhD in Bioinformatics from <a href="https://www.sorbonne-universite.fr" target="_blank" >Sorbonne Université</a >. He is currently a Senior Developer Advocate at <a href="https://seqera.io" target="_blank" >Seqera</a> and a member in the <a href="https://nf-co.re/" target="_blank" >nf-core</a > <a href="https://nf-co.re/governance#core-team" target="_blank" >core team.</a>',
  },
  {
    name: "Geraldine Van der Auwera",
    img: "geraldine.png",
    country: "us",
    github: "vdauwera",
    linkedin: "geraldine-van-der-auwera-5a5811",
    twitter: "VdaGeraldine",
    title: "Ambassador Program Staff",
    children:
      'Geraldine is Lead Developer Advocate at <a href="https://seqera.io" target="_blank" >Seqera</a >. She has a bit of history with GATK and is also the author of Genomics in the Cloud (O\'Reilly, 2020).',
  },
  {
    name: "Ken Brewer",
    img: "kbrewer.jpg",
    country: "us",
    github: "kenibrewer",
    linkedin: "kenibrewer",
    twitter: "kenibrewer",
    title: "Ambassador Program Staff",
    children:
      "Ken is a open-source developer and advocate active within the areas of genomics and high-content imaging microscopy. He is a Senior Developer Advocate at Seqera.",
  },
  // Program Alumni
  {
    name: "Alex Peltzer",
    img: "alex.jpeg",
    country: "de",
    github: "apeltzer",
    linkedin: "apeltzer",
    twitter: "alex_peltzer",
    title: "Program Alumni",
  },
  {
    name: "André M. Ribeiro-dos-Santos",
    img: "andre_santos.jpeg",
    country: "br",
    github: "andremrsantos",
    linkedin: "am-ribeiro-dos-santos",
    twitter: "andremrsantos",
    title: "Program Alumni",
    children:
      "André is a Bioinformatician with expertise in gene regulation, population genetics, and machine learning. Currently working as visiting professor in Universidade Federal do Pará (UFPA, Belem - PA, BR).",
  },
  {
    name: "Anders Jemt",
    img: "anders.jpeg",
    country: "se",
    github: "jemten",
    linkedin: "anders-jemt-16b73367",
    title: "Program Alumni",
  },
  {
    name: "Bhargava Morampalli",
    img: "bhargava.jpeg",
    country: "nz",
    github: "bhargava-morampalli",
    linkedin: "bhargavamorampalli",
    twitter: "iambhargava",
    title: "Program Alumni",
  },
  {
    name: "Carson Miller",
    img: "CarsonMiller.png",
    country: "us",
    github: "CarsonJM",
    linkedin: "carson-miller-a94911173",
    twitter: "carsonjmiIIer",
    title: "Program Alumni",
    children:
      'Carson Miller is a Microbiology Ph.D. Candidate in Lucas Hoffman\'s lab at the <a href="https://microbiology.washington.edu/" target="_blank" >University of Washington</a >. He is also a <a href="https://nf-co.re/governance#maintainers" target="_blank" >maintainer</a > in the <a href="https://nf-co.re" target="_blank" >nf-core</a> project involved in developing/updating pipelines related to microbial metagenomics.',
  },
  {
    name: "Christopher Mohr",
    img: "christopher.jpeg",
    country: "de",
    github: "christopher-mohr",
    linkedin: "christopher-mohr-8b6a52197",
    twitter: "cmohr_tue",
    title: "Program Alumni",
  },
  {
    name: "Colby T. Ford",
    img: "Colby_Ford.jpg",
    country: "us",
    github: "colbyford",
    linkedin: "colbyford",
    title: "Program Alumni",
    children:
      '<a href="https://colbyford.com" target="_blank">Colby T. Ford, Ph.D.</a> is a Principal Consultant at <a href="https://tuple.xyz" target="_blank" >Tuple</a >. He is the author of <a href="https://www.oreilly.com/library/view/genomics-in-the/9781098139032/" target="_blank" >Genomics in the Azure Cloud</a > (O\'Reilly, 2022) and specializes in cloud architecture for large bioinformatics workloads.',
  },
  {
    name: "Edwin Simjaya",
    img: "edwin.jpg",
    country: "id",
    linkedin: "edwin-simjaya",
    github: "KalbeDigitalLab",
    title: "Program Alumni",
    children:
      'Edwin Simjaya is the Head of AI & Software Center at <a href="https://www.kalbe.co.id" target="_blank" >Kalbe Farma</a >. His focus is on AI-augmented nutrigenomics using LLMs and is the driving force behind projects at <a href="https://kdu.kalbe.co.id" target="_blank" >Kalbe Digital University</a >.',
  },
  {
    name: "Florian Wünnemann",
    img: "fwuennemann.jpg",
    country: "de",
    github: "flowuenne",
    linkedin: "florian-wünnemann-322357189",
    twitter: "flowuenne",
    title: "Program Alumni",
    children:
      'Florian is a Bioinformatics Engineer at <a href="https://seqera.io/" target="_blank" >Seqera</a >, living in Montreal, with experience in single-cell and spatial omics technologies. He is a developer for <a href="https://nf-co.re/molkart/dev" target="_blank" >nf-core/molkart</a > and in the nf-core maintainers team.',
  },
  {
    name: "Gisela Gabernet",
    img: "gisela.jpeg",
    country: "us",
    github: "ggabernet",
    linkedin: "gisela-gabernet-89a9779b",
    twitter: "GGabernet",
    title: "Program Alumni",
  },
  {
    name: "Jasmin Frangenberg",
    img: "jasmin.jpeg",
    country: "de",
    github: "jasmezz",
    linkedin: "jasmin-frangenberg",
    title: "Program Alumni",
  },
  {
    name: "Jaykishan Solanki",
    img: "Jaykishan.jpg",
    country: "in",
    github: "jaykishan91",
    linkedin: "jaykishan-solanki-7615b2146",
    twitter: "JaykishanSolan9",
    title: "Program Alumni",
    children:
      "Jaykishan Solanki is an Bioinformatics Developer at NCGM Global working on omics data analysis and machine learning.",
  },
  {
    name: "Juan A. Ugalde",
    img: "jugalde.jpg",
    country: "cl",
    github: "juanu",
    linkedin: "juanugalde",
    twitter: "juanu",
    title: "Program Alumni",
    children:
      'Juan is an <a href="https://www.ugaldelab.org" target="_blank" >assistant professor</a> at the Center for Bioinformatics and Integrative Biology, at Universidad Andres Bello, Santiago, Chile.',
  },
  {
    name: "Juan Martinez Villalobos",
    img: "martinezvbs.jpeg",
    country: "mx",
    github: "martinezvbs",
    title: "Program Alumni",
    children:
      'Juan is a Postgraduate Associate at the Augert Lab at <a href="https://medicine.yale.edu/profile/juan-martinezvillalobos/" target="_blank" >Yale School of Medicine</a > His research focuses on CRISPR screening and multi-omics integration to characterize the epigenetic landscape of Small Cell Lung Cancer (SCLC).',
  },
  {
    name: "Krešimir Beštak",
    img: "kbestak.jpg",
    country: "de",
    github: "kbestak",
    linkedin: "krešimir-beštak-a027ab211",
    twitter: "KresBestak",
    title: "Program Alumni",
    children:
      'Krešimir is a PhD student at the <a href="https://www.medizinische-fakultaet-hd.uni-heidelberg.de/einrichtungen/institute/institute-for-computational-biomedicine" target="_blank" >Institute for Computational Biomedicine</a >, University Hospital Heidelberg. He\'s the developer of <a href="https://nf-co.re/molkart/dev" target="_blank" >nf-core/molkart</a >, and is working on spatial omics nf-core pipelines.',
  },
  {
    name: "Melak Getu",
    img: "melak.jpg",
    country: "et",
    github: "MelakG13",
    title: "Program Alumni",
    linkedin: "melak-getu-531095193",
    bluesky: "https://bsky.app/profile/melak13.bsky.social",
    children:
      'Melak Getu is a Bioinformatics Researcher at <a href="https://www.ephi.gov.et/" target="_blank" >EPHI</a >. He focuses specifically in viral and bacterial genome analysis, and works at the Integrated Pathogen Genomic Sequencing and Bioinformatics Facility of EPHI.',
  },
  {
    name: "Mitali Nimkar",
    img: "mitali_nimkar.png",
    country: "in",
    github: "turningbrain",
    linkedin: "mitali-nimkar-40714438",
    title: "Program Alumni",
  },
  {
    name: "Moritz E. Beber",
    img: "moritz.jpeg",
    country: "dk",
    github: "Midnighter",
    linkedin: "moritz-beber-b597a55a",
    twitter: "me_beber",
    title: "Program Alumni",
    children:
      'Moritz is a Senior IT Solution Architect at <a href="https://www.novonordisk.com" target="_blank" referrerpolicy="no-referrer" >Novo Nordisk A/S</a >. He holds a PhD in bioinformatics and worked in computational systems biology. He strongly supports open knowledge and open source initiatives such as <a href="https://nf-co.re/" target="_blank" referrerpolicy="no-referrer" >nf-core</a >.',
  },
  {
    name: "Natalia Coutouné",
    img: "natalia.jpeg",
    country: "br",
    github: "nat-coutoune",
    linkedin: "natalia-coutouné-a421b91b5",
    title: "Program Alumni",
  },
  {
    name: "Nicolas A. Schcolnicov",
    img: "nschcolnicov.jpg",
    country: "ar",
    github: "nschcolnicov",
    linkedin: "naschcolnicov",
    title: "Program Alumni",
    children:
      'Nicolas Schcolnicov is a Bioinformatics Consultant focused on the development and maintenance of Nextflow pipelines. He contributes to the <a href="https://nf-co.re/" target="_blank" >nf-core</a> workflows and his interests are scientific software development and data science.',
  },
  {
    name: "Nico Trummer",
    img: "NicoTrummer.JPG",
    country: "de",
    github: "nictru",
    title: "Program Alumni",
    linkedin: "nico-trummer-401860199",
    children:
      'Nico is a student research assistant at the <a href="https://www.mls.ls.tum.de/en/daisybio/home/" target="_blank" >DaiSyBio group</a> at TUM in Munich, Germany. He is a member of the <a href="https://nf-co.re/governance#maintainers" target="_blank" >nf-core maintainers team</a > and leads the development of several nf-core pipelines.',
  },
  {
    name: "Noah A. Legall",
    img: "noahlegall.JPG",
    country: "us",
    github: "noahaus",
    linkedin: "noah-legall-phd-aa424210a",
    title: "Program Alumni",
    children:
      "Noah Legall (Ph. D.) is a Bioinformatic Scientist who uses Nextflow to make comparative genomics of zoonotic and antimicrobial resistant bacterial pathogens seamless.",
  },
  {
    name: "Opeyemi Alayande",
    img: "opeyemi_alayande.jpg",
    country: "ng",
    github: "ycode-sh",
    linkedin: "yemistats",
    twitter: "@ycode_sh",
    title: "Program Alumni",
    children:
      'Opeyemi creates and uses Nextflow-powered tools for infectious disease management at the <a href="https://dfg-mrtc.org/" target="_blank" >Damien Foundation Genomics and Mycobacteria Research and Training centre</a >, located in the <a href="http://www.ui.edu.ng/" target="_blank" >University of Ibadan</a >.',
  },
  {
    name: "Priyanka Surana",
    img: "priyanka.jpeg",
    country: "gb",
    github: "priyanka-surana",
    linkedin: "priyanka-surana",
    twitter: "psurana",
    title: "Program Alumni",
  },
  {
    name: "Ramprasad Neethiraj",
    img: "ramprasadneethiraj.jpg",
    country: "se",
    github: "ramprasadn",
    linkedin: "ramprasad-neethiraj",
    title: "Program Alumni",
    children:
      'Ram is a bioinformatician in the <a href="https://www.scilifelab.se/units/clinical-genomics-stockholm/">Clinical Genomics Unit</a> at <a href="https://www.scilifelab.se/">SciLifeLab</a>, Stockholm, developing diagnostic pipelines for clinical settings. He has a Ph.D. in Population Genetics from <a href="https://www.su.se/english/">Stockholm University.</a>',
  },
  {
    name: "Robert Petit",
    img: "robertpetit.png",
    country: "us",
    github: "rpetit3",
    linkedin: "rpetit3",
    twitter: "rpetit3",
    title: "Program Alumni",
    children:
      'Robert is a Senior Bioinformatics Scientist at the <a href="https://health.wyo.gov/publichealth/lab/" target="_blank" >Wyoming Public Health Laboratory</a >. He\'s the developer of <a href="https://bactopia.github.io/" target="_blank" >Bactopia</a >, <a href="https://bactopia.github.io/v3.0.0/impact-and-outreach/enhancements/#nf-coremodules-contributions" target="_blank" >a regular contributor to nf-core</a >, and a member of the Bioconda Core team.',
  },
  {
    name: "Sateesh Peri",
    img: "SateeshPeri.png",
    country: "in",
    github: "sateeshperi",
    linkedin: "sateesh-peri-b0130476",
    title: "Program Alumni",
  },
  {
    name: "Solenne Correard",
    img: "Solenne.png",
    country: "fr",
    github: "scorreard",
    linkedin: "solenne-correard-2a4675114",
    twitter: "CorreardSolenne",
    title: "Program Alumni",
    children:
      'Solenne Correard is the scientific director of <a href="https://www.wiseancestors.org/">Wise Ancestors.</a> She is a molecular biologist and bioinformatician using nextflow for projects focusing on biodiversity conservation and equity in health-care.',
  },
  {
    name: "Sombuddha Roy Bhowmick",
    img: "Sombuddha_Roy_Bhowmick_Image.jpg",
    country: "in",
    github: "Sombuddha-Roy-Bhowmick",
    linkedin: "sombuddha-roy-bhowmick-b76b1a13a",
    title: "Program Alumni",
    children:
      "Sombuddha Roy Bhowmick is a Bioinformatician/Computational Biologist trying to make a difference in a patient's life, one at a time.",
  },
];

const AmbassadorsList: React.FC = () => {
  const [selectedCountries, setSelectedCountries] = useState<string[]>([]);

  const ambassadorsByCategory = useMemo(() => {
    const regular = ambassadors.filter((amb) => !amb.title || amb.title === "Nextflow Ambassador");
    const staff = ambassadors.filter((amb) => amb.title === "Ambassador Program Staff");
    const alumni = ambassadors.filter((amb) => amb.title === "Program Alumni");

    return { regular, staff, alumni };
  }, []);

  const filteredCategories = useMemo(() => {
    if (selectedCountries.length === 0) {
      return ambassadorsByCategory;
    }

    // Function to check if an ambassador should be included based on selected countries
    const shouldIncludeAmbassador = (ambassador: Ambassador) => {
      if (selectedCountries.includes("es") && ambassador.country === "es-ct") {
        return true;
      }
      return selectedCountries.includes(ambassador.country);
    };

    return {
      regular: ambassadorsByCategory.regular.filter(shouldIncludeAmbassador),
      staff: ambassadorsByCategory.staff.filter(shouldIncludeAmbassador),
      alumni: ambassadorsByCategory.alumni.filter(shouldIncludeAmbassador),
    };
  }, [selectedCountries, ambassadorsByCategory]);

  const handleFilterChange = (countries: string[]) => {
    setSelectedCountries(countries);
  };

  const renderAmbassadorGrid = (ambassadorsArray: Ambassador[], defaultTitle: string = "Nextflow Ambassador") => (
    <div className="grid grid-cols-1 sm:grid-cols-2 md:grid-cols-3 lg:grid-cols-4 xl:grid-cols-5 gap-4 w-full mx-auto">
      {ambassadorsArray.map((ambassador, index) => (
        <AmbassadorCard
          key={`${ambassador.name}-${index}`}
          name={ambassador.name}
          img={ambassador.img}
          country={ambassador.country}
          github={ambassador.github || ""}
          linkedin={ambassador.linkedin}
          twitter={ambassador.twitter}
          mastodon={ambassador.mastodon}
          bluesky={ambassador.bluesky}
          title={ambassador.title || defaultTitle}
        >
          {ambassador.children && <div dangerouslySetInnerHTML={{ __html: ambassador.children }} />}
        </AmbassadorCard>
      ))}
    </div>
  );

  const totalFiltered =
    filteredCategories.regular.length + filteredCategories.staff.length + filteredCategories.alumni.length;

  return (
    <div>
      <AmbassadorFilter onFilterChange={handleFilterChange} />

      {filteredCategories.regular.length > 0 && (
        <div className="container container-lg mb-16">{renderAmbassadorGrid(filteredCategories.regular)}</div>
      )}

      {filteredCategories.staff.length > 0 && (
        <div>
          <h1 id="program-staff" className="text-center text-4xl font-display my-12">
            Ambassador Program Staff
          </h1>
          <div className="container container-lg mb-16">
            {renderAmbassadorGrid(filteredCategories.staff, "Ambassador Program Staff")}
          </div>
        </div>
      )}

      {filteredCategories.alumni.length > 0 && (
        <div>
          <h1 id="alumni" className="text-center text-4xl font-display my-12">
            Program Alumni
          </h1>
          <p className="text-center mb-8">
            Individuals who are no longer active Nextflow Ambassadors but have participated in the program in the past.
          </p>
          <div className="container container-lg mb-16">
            {renderAmbassadorGrid(filteredCategories.alumni, "Program Alumni")}
          </div>
        </div>
      )}

      {totalFiltered === 0 && selectedCountries.length > 0 && (
        <div className="text-center mt-8">
          <p className="">No ambassadors found in the selected countries.</p>
        </div>
      )}
    </div>
  );
};

export default AmbassadorsList;
