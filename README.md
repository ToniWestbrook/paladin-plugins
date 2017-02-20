# paladin-plugins
Pipeline plugins for PALADIN, providing HPC support, abundance (taxonomy, go terms), automation, etc.

PALADIN-plugins can be chained in a pipeline configuration using a single or multiple execution runs to process PALADIN related data through multiple stages of analysis.  Included plugins are as follows:

- Aggregation: Combined multiple PALADIN outputs (SAM output and UniProt report) into a single output
- Automate: Automate PALADIN executions across multiple sets of reads (recursively searches directories using regex patterns)
- Decluster: Generate a reference database containing each protein sequence used by UniProt to generate each UniRef cluster (during the CD-HIT clustering process) for those clusters detected in a PALADIN UniProt report.  This may be used for a second PALADIN pass, refining alignments from representative hits to member hits
- Difference: Analyze relative differences between two PALADIN taxonomy reports.  Can also report specific contributing reads for each difference
- GO: Perform gene ontology term grouping and abundance reporting
- HPC: Distribute PALADIN execution across cluster nodes
- Plotting: Generate plots in PNG format from pipeline generated data
- Taxonomy: Perform taxonomic grouping and abundance reporting
- Uniprot: Download custom UniProt reports for a PALADIN prepared SAM alignment

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
