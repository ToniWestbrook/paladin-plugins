# paladin-plugins
Pipeline plugins for PALADIN, providing HPC support, abundance (taxonomy, go terms), automation, etc.

PALADIN-plugins can be chained in a pipeline configuration using a single or multiple execution runs to process PALADIN related data through multiple stages of analysis.  Included plugins are as follows:

- Aggregation: Combine multiple PALADIN outputs (SAM output and UniProt report) into a single output
- Automate: Automate PALADIN executions across multiple sets of reads (recursively searches directories using regex patterns)
- Decluster: Generate a reference database containing each protein sequence used by UniProt to generate each UniRef cluster (during the CD-HIT clustering process) for those clusters detected in a PALADIN UniProt report.  This may be used for a second PALADIN pass, refining alignments from representative hits to member hits
- Difference: Analyze relative differences between two PALADIN taxonomy reports.  Can also report specific contributing reads for each difference
- GO: Perform gene ontology term grouping and abundance reporting
- HPC: Distribute PALADIN execution across cluster nodes
- Plotting: Generate plots in PNG format from pipeline generated data
- Taxonomy: Perform taxonomic grouping and abundance reporting
- Uniprot: Download custom UniProt reports for a PALADIN prepared SAM alignment
- Internal Commands: Flush, Write (screen and file output)

Plugins in the pipeline are prepended with "@@", and are then given plugin specific arguments (which may be enumerated for any plugin with the "-h" argument).  Multiple plugin calls may be specified in a single execution (see examples below).

INSTALLATION
--
**Dependencies**

- PALADIN
- Python 3
- Python libraries: matplotlib, mpi4pyi, numpy, requests

```
git clone https://github.com/twestbrookunh/paladin-plugins.git
```

SAMPLE COMMANDS
--

List available plugins.
```
paladin-plugins -l
```

Run PALADIN on all fastq (fq) files in a subdirectory that have "R1" in the filename using 16 threads, placing the output in each subdirectory:
```
paladin-plugins @@automate reference.fasta.gz /path/subdir .*R1.*\.fq.gz -t 16 -o outfile
```
Combine the SAM and TSV outputs of the multiple PALADIN runs above into a single output
```
paladin-plugins @@aggregation -r /path/subdir -s outfile\.sam -o combined
```

Group GO terms and write abundances to a file, filtering for a mapping quality of 20
```
paladin-plugins @@go -i input.tsv -q 20 @@write output.txt
```

Group all bacterial species and write abundances to a file, then plot data to a pie chart, filtering for a mapping quality of 30 and limited number of values shown on graph to 10

```
paladin-plugins @@taxonomy -i input.tsv -q 30 -t species -r Bacteria @@write taxonomy.txt @@plotting -i taxonomy.txt -o chart.png -l 10 -L "My Chart" "My X-Axis", "My Y-Axis" -p -s 12 12
```

Group flattened kingdoms (child ranks of domains [level 1]) to one file, group all species to second file, then show both charts side-by-side in a single PNG:
```
paladin-plugins @@taxonomy -i input.tsv -q 30 -t children -l 1 @@write first.txt @@taxonomy -i input.tsv -q 30 -t species -l 0 @@write second.txt @@plotting -s 20 15 -g 1 2 @@plotting -i first.txt -l 10 -c 0 0 @@plotting -i second.txt -l 10 -c 0 1 -o chart.png
```

HPC Plugin
--
The HPC plugin may be used to distribute a PALADIN execution across multiple cluster nodes.  It is used in conjunction with mpirun and the appropriate batch submission script. Below is an example Slurm batch script for splitting a PALADIN execution across 10 nodes, using 24 threads per node, and routing any communication with uniprot.org through a central proxy server on port 3128:

```
#!/bin/bash

#SBATCH -p target_partition
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=10
#SBATCH --job-name="PALADIN"
#SBATCH --output=running.log

mpirun paladin-plugins.py @@hpc reference.fasta.gz reads.fq.gz output -t 24 -P http://proxy:3128
```

Note, usage of the HPC plugin may display an error if not running in an MPI aware environment.  
