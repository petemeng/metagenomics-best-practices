# Benchmarking second and third-generation sequencing platforms for microbial metagenomics #

### Purpose ###
We compared seven platforms, encompassing second generation short read sequencers (Illumina Hiseq300, DNBSEQ-G400, DNBSEQ-T7, Ion S5 and Ion Proton P1) and third generation long read sequencers (Oxford Nanopore Technologies MinION and Pacific Biosciences PacBio) on three synthetic microbial communities composed of up to 87 microbial species per mock and spanning 29 phyla. 


### Requirements ###
* [fastp](https://github.com/OpenGene/fastp) v0.20.0
* [AlienTrimmer](https://github.com/aghozlane/masque/blob/master/AlienTrimmer_0.4.0/src/AlienTrimmer.java) v2.0
* [bowtie2](https://github.com/BenLangmead/bowtie2) v2.3.5.1
* [minimap2](https://github.com/lh3/minimap2) v2.15-r915-dirty
* [SPAdes](https://github.com/ablab/spades) v3.14.1
* metaFlye v2.8.1-b1688
* [metaquast](https://github.com/ablab/quast/blob/master/metaquast.py) v4.6.3
* Python3

The following R librairies are required:
```
library(ggpubr)
library(reshape2)
library(dplyr)
library(stringr)
```

### Authors ###
* Mathieu Almeida: mathieu.almeida@inrae.fr
* Victoria Meslier: victoria.meslier@inrae.fr
* Kevin Da Silva: kevin.dasilva.fr@gmail.com
* Florian Plaza Onate: florian.plaza-onate@inrae.fr


### Funding ###
This work was supported by the European FP7 Marie Sk?odowska-Curie actions AgreenSkillsPlus PCOFUND-GA-2013-609398 grant to MA. Additional funding was from the Metagenopolis grant ANR-11-DPBS-0001. MP was supported by Genomic Science Program, U.S. Department of Energy, Office of Science, Biological, and Environmental Research as part of the Plant Microbe Interfaces Scientific Focus Area, and by grant R01DE024463 from the National Institute of Dental and Craniofacial Research of the US National Institutes of Health. Oak Ridge National Laboratory is managed by UT-Battelle, LLC, for the U.S. Department of Energy under contract DE-AC05-00OR22725. 


### Acknowlegements ###
We thank Dr. Cynthia Gilmour (Smithsonian Research Institute, Edgewater, Maryland, USA) for providing some of the purified gDNA used in this study. We also thank David Stucki, Deborah Moine and Jules Bourgon (Pacific Biosciences, Menlo Park, California, USA) for PacBio sequencing, Jihua Sun and Yong Hou (BGI, Shenzhen, Guangdong, China) for DNBseq sequencing and Clemence Genthon (Plateforme Génomique GeT INRAE Transfert, Toulouse, France) for Illumina sequencing. 
