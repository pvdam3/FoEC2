# Fusarium effectors

This pipeline can identify putative effectors in a provided set of *Fusarium oxysporum* genomes and show their presence/absence variation across all input genomes.

## Concept
The following steps are executed in this pipeline:

**1. Putative effector identification (per genome)**
  * Miniature impala(mimp) terminal inverted repeat (TIR) identification based on TIR consensus sequence.
  * Determination of genome regions downstream of mimps/TIRs.
  * Finding possible Open Reading Frames (ORFs) within these downstream regions by translating the sequence in three or six frames (default three frames) and finding the first Methionine (M) residue followed by a sequence of threshold length (default 10 amino acids) and a STOP codon/end of contig
  * Update gene models (when possible) using [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) gene prediction. 
  * Filter for secreted proteins with SignalP.
  * Filter for size (default between 10 and 600 amino acids) and cysteine content (default 0 cysteines).

**2. Duplicate effector candidates are removed:**
  * A BLAST database is created from all the putative effectors (from all input genomes) and each of the putative effectors are BLASTed against this database using [DIAMOND](https://github.com/bbuchfink/diamond). The similarity among putative effectors is shown in their BLAST hits. Clusters are created with this information using [MCL](https://github.com/micans/mcl), thereby essentially marking redundancy.
  * A Hidden Markov Model (HMM) profile is created per putative effector cluster using [HMMER](http://hmmer.org/). This profile captures the variability present in each cluster and is used to represent the cluster.

**3. Identifying presence-absence patterns in the genomes:**
  * A search for each putative effector HMM profile across all input genomes is conducted with HMMER's nHMMER. These hits, found in (`output/03.presenceabsence/00_genome_effector_hits.out`), are filtered by a minimum E-value and query coverage (default thresholds 10E-10 and 0.8 respectively). The total number of valid hits per genome is recorded. 
  * These values are stored in a table with the effectors on one axis and the genomes on the other.
  * Putative effectors with excessive hits are filtered out (default threshold 20 hits for a single genome and 5 hits on average across all genomes). The final table can be seen at `output/03.presenceabsence/01_presence_absence.tsv`.

**4. Visualization:**
  * To discover which genomes are most alike in terms of effector pallette, the table is imported in an R script, which applies hierarchical clustering on the rows and columns. The clustering method can be selected in the RShiny App.
  * Dendrograms representing this clustering can be visualized in the app, and their newick files can be downloaded.
  * Multiple sequence alignments (MSAs) found in `output/03.presenceabsence/putative_effector_msas` can be visualized in the app.

## Usage
Please first make sure to have all the dependencies installed (see below).
This pipeline requires that the paths to the input genomes be in a specific configuration file. In order to facilitate this, the `foec_setup.sh` script automatically creates this config file. Usage, where `infolder` is the path to a directory which only contains the input genomes in FASTA format: 
```bash
./foec_setup.sh -g [infolder]
```

FoEC2 can also be run with an existing list of effectors. To do this, FASTA files representing these effectors will need to be placed into a separate directory `effectors`. The pipeline will then skip the putative effector prediction steps. Usage:
```bash
./foec_setup.sh -g [infolder] -e [effectors]
```

Type `./foec_setup.sh -h` for a detailed help page including options.

More configuration options for the pipeline are available in the config file [config/config.yaml](config/config.yaml). 

Once all the genome configuration files are ready, the pipeline can be run. Usage:
```bash
snakemake --use-conda --cores [N]
```

Type `snakemake -h` for a detailed help page for more Snakemake related options.

Once the pipeline has finished, the output files can be visualized in the RShiny app (`app.R`). Usage (R):
```R
library(shiny)
runApp("scripts/app.R")
```

The app will open at the 'Data' tab, where the output from the pipeline can be uploaded:
* PAV table: the presence/absence variation table (`output/03.presenceabsence/01_presence_absence.tsv`)
* Genome metadata table: extra information concerning genomes to include in the visualization (i.e. *formae speciales*, sample location, etc.). A template can be found at `config/visualization_config.csv` after running the pipeline.
* Putative effector metadata table: extra information concerning putative effectors to include in the visualization (i.e. suspected SIX genes, other genes of interest). A template can be found at `config/visualization_config_effectors.csv` after running the pipeline.

The plots can be seen in the 'Plots' tab:
* Heatmap: this plot depicts the presence / absence patterns of putative effectors across genomes.
* Dendrograms: these dendrograms reflect those which are shown on each axis of the heatmap.
* MSAs (Multiple Sequence Alignments): per putative effector cluster.

## Dependencies
The pipeline relies a number of different 3rd party programs and libraries:
* [Conda](https://docs.conda.io/en/latest/)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
* [Singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps) (only required to run FoEC2 with a container - recommended for users experiencing problems with installing conda environments)
* [SignalP](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp) (only versions 4 and 5 are currently supported)
* R with the following libraries installed:
  * shiny
  * shinythemes
  * dendextend
  * RColorBrewer
  * pals
  * pheatmap
  * phylogram
  * DT
  * rhandsontable
  * [msaR](https://zachcp.github.io/msaR/)
  ```R
  install.packages(c("shiny", "shinythemes", "dendextend", "RColorBrewer", "pals", "pheatmap", "phylogram", "DT", "rhandsontable", "msaR"))
  ```
* Python with the following package installed:
  * [BioPython](http://biopython.org/wiki/Download)

## References
van Dam, P., Fokkens, L., Schmidt, S. M., Linmans, J. H. J., Kistler, H. C., Ma, L.-J., & Rep, M. (2016). Effector profiles distinguish <I>formae speciales </I>of <I>Fusarium oxysporum</I>. Environmental Microbiology. http://doi.org/10.1111/1462-2920.13445