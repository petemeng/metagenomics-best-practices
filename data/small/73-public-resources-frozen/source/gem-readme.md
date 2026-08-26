## Supporting data for: A Genomic Catalogue of Earth's Microbiomes


The reconstruction of bacterial and archaeal genomes from shotgun metagenomes has enabled unprecedented insights into the ecology and evolution of environmental and host-associated microbiomes. Here we applied this powerful approach to over 10,000 metagenomes collected from diverse habitats covering all of Earth's continents and oceans, human- and animal-host associated microbiomes, engineered environments, and natural and agricultural soils to capture extant microbial metabolic and functional potential. We present a comprehensive catalogue of 52,515 metagenome-assembled genomes representing 12,556 novel candidate species-level operational taxonomic units (OTUs), spanning 135 phyla, which expand the known phylogenetic diversity of Bacteria and Archaea by 44%. We also demonstrate the utility of this collection for secondary metabolite biosynthetic potential and predicting host-virus linkage, which can provide a view into the global distribution of lysogenic viruses. This resource underscores the value of leveraging genome-centric approaches to reveal genomic properties of uncultivated microbes that impact on ecosystem processes.

### Data usage policy
The MAGs contained in the GEM data bundle are free to use. Underlying metagenomes are protected by JGI data release and utilization policies found here: https://jgi.doe.gov/user-programs/pmo-overview/policies/#data-util, and we encourage contacting PIs for any planned analyses or publications that may overlap with existing project goals.


#### --genomes/
* 52,515 MAGs from the current study
* CDS predicted using Prodigal (genomes/faa, genomes/ffn)
* Quality statistics, environmental & geographic metadata, OTU assignments, and other info found in genome_metadata.tsv


#### --otus/
* 52,515 MAGs and 524,046 reference genomes were clustered into 45,599 species-level OTUs
* All genomes were evaluated using CheckM and were required to have a quality score > 50
* All OTUs were taxonomically annotated using GTDB-tk
* Representive genomes were selected with the highest CheckM quality scores with isolate genomes prioritized over MAGs and SAGs -- see the /genomes directory for detailed information

#### --tree/
* A phylogenetic tree was constructed for a subset of 43,979 species-level OTUs
* Based on a concatenated alignment of 30 universal, single-copy marker genes
* Tree constructed using FastTree


#### --bgcs/ 
* Biosynthetic gene clusters identified in 52,515 MAGs using antiSMASH v5.1
* Any contig <5 kb was not considered.


#### --prophages/
* 24,313 putative proviruses were identified in the 52,515 MAGs using VirSorter v1.0.3
* All predictions of categories 4 and 5 were retained


#### --protclusts/  
* CDS from the 52,515 MAGs were grouped into 5,794,145 protein clusters using MMseqs2

