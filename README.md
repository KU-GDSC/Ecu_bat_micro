# Bioinformatics Pipeline for Bat Pathogens

*Identification from Host-Depleted Nanopore Sequencing Data (Optimized for bacterial detection of 
known human pathogens in bat tissues).*

The pipeline described in this tutorial requires that you are working 
from a Unix-based command-line terminal. If your computer is running
a modern Mac or Linux operating system then you should be able to follow
all of these instructions without additional setup. If your computer is
running Windows 10 or Windows 11, you can follow [these instructions from
Microsoft to install the Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install), which will install a Unix-based
command-line terminal for Windows.

## Step 0: Installations

### Dorado

A basecalling algorithm developed by ONT for decoding the raw electrical signals produced from sequencing, and translating them into nucleotide sequences. It also outputs a quality score for each base, which will be used later. [Follow these instructions to download the latest version of Dorado](https://software-docs.nanoporetech.com/dorado/1.4.0/#installation) for your operating system.

### Mamba / Conda

[Mamba](https://mamba.readthedocs.io/en/latest/) is a reimplementation of
the popular package management system [Conda](https://anaconda.org/channels/anaconda/packages/conda/overview), which is much faster and works with
exactly the same code syntax.

[Follow these directions to install mamba](https://github.com/conda-forge/miniforge#unix-like-platforms-macos-linux--wsl).

In short, first download the installation script:

```{bash}
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```

And then run the installation script, which will begin an interactive session
to install mamba on your computer.

```{bash}
bash Miniforge3-$(uname)-$(uname -m).sh
```

Then, to set up the mamba environment and download all of the remaining
programs we'll need for this pipeline, use the following command:

```{bash}
mamba create -n bat-micro -f setup/bat-micro.yml
```

This will start an interactive series of prompts that will walk you
through downloading all the necessary programs.

The programs installed can now be accessed by activating the mamba environment

```{bash}
mamba activate bat-micro
```

Note: you will need to activate the environment each time you open a new
terminal window.

## Step 1: Basecalling & Quality Control

### Basecalling

This step can be very computationally expensive and time-consuming. 

```{bash}
dorado basecaller sup@v5.2.0 pod5_pass/ --kit-name SQK-NBD114-24 --min-qscore 10 --recursive > mepa_pathogens_v5.2.0.bam
```

This process can be run substantially faster on a computer with GPU
resources, in which case add the parameter `-x cuda:all` to the above
command before `> mepa_pathogens_v5.2.0.bam`.

For comparison, basecalling on these data completed in 2.5 hours using
GPUs on a high-performance computing cluster with 128G of memory and 10
CPUs.

### Demultiplexing

The basecalled BAM file contains data from all samples, this command
will demultiplex the base calls into FASTQ files separately by barcode
(note: the `--no-classify` parameter is used because the barcodes were
identified in the previous step, and trimmed)

```{bash}
dorado demux mepa_pathogens_v5.2.0.bam -o mepa_pathogens_v5.2.0 --emit-fastq --emit-summary --no-classify
```

The resulting directory structure is a bit cumbersome, so we can make
some links to the FASTQ data for convenience:

```{bash}
mkdir -p reads/qc
pushd reads
for i in `find ../mepa_pathogens_v5.2.0 -name "*.fastq"`; do ln -s ${i} .; done
popd
```

### Host Sequence Filtering

Even though our data was generated using host-depleted sequencing,
the host genome can still make it through. Because of this, it is
necessary to remove these sequences before proceeding. For this we
will be using the latest reference genome of *Desmodus rotundus*.

1. Download the reference genome from NCBI

```{bash}
datasets download genome accession GCF_022682495.2 --include gff3,rna,cds,protein,genome,seq-report
```

2. Decompress the reference data file

```{bash}
unzip ncbi_data.zip
```

3. Filter the read data [following a method described at this link](https://linsalrob.github.io/ComputationalGenomicsManual/Deconseq/).

First we align the reads for each barcode to the *Desmodus rotundus* genome

```{bash}
for i in $(seq -w 01 24)
  do
    barcode="barcode${i}"
    minimap2 --split-prefix=tmp$$ -a -xsr ncbi_dataset/data/GCF_022682495.2/GCF_022682495.2_HLdesRot8A.1_genomic.fna reads/FAZ94206_pass_${barcode}_1608bcd0_00000000_0.fastq | samtools view -bh | samtools sort -o reads/qc/${barcode}_host_aligned.bam
done
```

Then we use the [samtools flags](https://broadinstitute.github.io/picard/explain-flags.html) to filter the reads that mapped to the host:

```{bash}
for i in $(seq -w 01 24)
  do
    barcode="barcode${i}"
    samtools fastq -F 3588 reads/qc/${barcode}_host_aligned.bam > reads/qc/${barcode}_host.fastq
done
```

as well as the reads that did not map to the host:

```{bash}
for i in $(seq -w 01 24)
  do
    barcode="barcode${i}"
    samtools fastq -F 3584 -f 4 reads/qc/${barcode}_host_aligned.bam > reads/qc/${barcode}_nonhost.fastq
done
```

**UNREVISED BELOW**

---

Each individual base-call is accompanied by another character which indicate the error probability for that base-call, called the Phred Quality score, which is calculated with the following formula:

Q = −10 × log10 p

Where p represents the estimated error probability. For example, a base-call with a Q score of 20 would have a probability of 1/100 of being incorrect. For a more in-depth explanation on this, read Ewing & Green (1998) (<http://genome.cshlp.org/cgi/pmidlookup?view=long&pmid=9521922>)

2.  To check any of the generated BAM files, run:

``` bash
samtools head -n 100 barcodexx.bam
```

Remember to replace the `barcodexx.bam` with your file-name.

add the bash script to run in the cluster. Add notes


### Data Visualization

3.  In order to clean up the data, first look into the quality of the reads, as well as their length, this can be done with NanoPlot. This package prefers fastq files, so first transform them using samtools:

``` bash
for i in $(seq -w 01 24)
  do
    barcode="barcode${i}"
    samtools fastq dorado_sup_out/barcode${i}.bam > samtools_fastq_out/barcode${i}.fastq
done
```

4.  Now, use NanoPlot to generate a report for each barcode:

``` python
for i in $(seq -w 01 24)
  do
    NanoPlot -t 8 --fastq samtools_fastq_out/barcode${i}.fastq --loglength --plots kde --title barcode${i} -o nanoplot_out/barcode${i}
done
```

`-t 8` makes the package run in eight threads, this can speed up the process with more powerful gpus.

`--fastq samtools_fastq_out/barcode${i}.fastq` gives the file type and the file name.

`--loglength` puts read lengths on a logarithmic scale.

`kde` creates a kernel density estimate plot (kde plot), which is a method for visualizing the distribution of the data, similar to a histogram (<https://seaborn.pydata.org/generated/seaborn.kdeplot.html>).

`--title barcode${i}` Puts the title on all plots.

`-o nanoplot_out/barcode${i}` puts all generated graphs in this directory.

### Data Cleaning

5.  Remove short and low-quality reads using chopper:

``` python
for i in $(seq -w 01 24)
  do
chopper --threads 8 -q 20 -l 20 -i samtools_filter_out/barcode${i}.fastq > chopper_out/barcode${i}_filtered.fastq
done
```

`--threads 8` makes the package run in eight threads, this can speed up the process with more powerful gpus.

`-q 20` sets the minimum Phred Quality Score to 20. Any sequence with an average score below this will get removed.

`-l 20` sets the minimum read length to 20bp.

`-i samtools_fastq_out/barcode${i}.fastq` tells the package the name of the file to filter.

## Step 2: Taxonomic Classification (UTI Pathogen Focus)

### De-novo Assembly

1.  Assemble contigs using metaFlye (WARNING: this step is relatively taxing on computer resources, it might take a while to run the command):

``` python
for i in $(seq -w 01 24)
  do
    barcode="barcode${i}"
flye --meta --read-error 0.03 --nano-hq chopper_out/barcode${i}_filtered.fastq --out-dir flye_out/barcode${i} --threads 8
done
```

`--meta` activates metaFlye mode.

`--read-error 0.03` specifies that data has already been filtered to Q20 (in step 1.5.).

`--nano-hq` specifies that the data was basecalled using Dorado's super accurate mode (in step 1.1.).

`chopper_out/barcode${i}_filtered.fastq` tells the package the name of the file to assemble.

`--out-dirflye_out/barcode${i}` specifies the output directory

`--threads 8` makes the package run in eight threads, this can speed up the process with more powerful gpus.

2.  Extract the assemblies from their enclosing folders:

``` python
for i in $(seq -w 01 24)
  do
    barcode="barcode${i}"
    assembly="assembly{i}"
cp flye_out/barcode${i}/assembly.fasta flye_assembly/assembly${i}.fasta
cp flye_out/barcode${i}/assembly_info.txt flye_assembly/assembly${i}_info.txt
done
```

### UTI Identification

Use Kraken2, Centrifuge, or Kaiju for taxonomic assignment:

``` console
bash
Copy
Edit
kraken2 --db kraken_db --threads 8 --report bacteria_report.txt --output bacteria_classified.txt --use-names filtered.fastq
```

#### Key UTI Pathogens to Check For:

-   *Escherichia coli* (UPEC - Uropathogenic *E. coli*)
-   *Klebsiella pneumoniae*
-   *Proteus mirabilis*
-   *Enterococcus faecalis*
-   *Staphylococcus saprophyticus*
-   *Pseudomonas aeruginosa*
-   *Morganella morganii*

### Visualization of UTI Pathogens

Generate interactive taxonomic plots using KronaTools.

``` console
bash
Copy
Edit
cut -f2,3 bacteria_classified.txt | ktImportTaxonomy -o bacteria_krona.html
```

## Step 3: UTI Pathogen Confirmation & Functional Annotation

### Targeted UTI Pathogen Screening

Use MetaPhlAn or PathoScope to refine UTI pathogen detection.

``` console
bash
Copy
Edit
metaphlan filtered.fastq --input_type fastq -o pathogen_abundance.txt
```

Alternative: Check for pathogenic UTI genes using MASH screen.

``` console
bash
Copy
Edit
mash screen -w -i 0.9 UTI_reference_db.msh filtered.fastq > mash_results.txt
```

### Antimicrobial Resistance (AMR) Screening for UTI Bacteria

Use ResFinder or CARD (Comprehensive Antibiotic Resistance Database).

``` console
bash
Copy
Edit
rgi main --input_sequence filtered.fastq --output rgi_output.txt --aligner DIAMOND --local
```

#### Look for AMR genes common in UTI pathogens, such as:

-   **Beta-lactam resistance**: *blaCTX-M, blaTEM, blaSHV*
-   **Fluoroquinolone resistance**: *gyrA, parC*
-   **Aminoglycoside resistance**: *aac(6')-Ib*
-   **Sulfonamide resistance**: *sul1, sul2*

### Virulence Factor Detection for UTI Pathogens

Use ABRICATE with the VFDB (Virulence Factor Database).

``` console
bash
Copy
Edit
abricate --db vfdb filtered.fastq > virulence_report.txt
```

#### Key virulence genes in UTI bacteria to look for:

-   *Escherichia coli* (*fimH, papG, sfa, iroN*)
-   *Proteus mirabilis* (*hpmA, mrpA*)
-   *Klebsiella pneumoniae* (*rmpA, yersiniabactin*)

## Step 4: Bacterial Genome Assembly & UTI Pathogen Strain Analysis

### De Novo Assembly of UTI Pathogen Genomes

Use Flye for assembling long-read bacterial genomes.

``` console
bash
Copy
Edit
flye --nano-raw filtered.fastq --out-dir assembly_output --genome-size 5m
```

### Polishing Assembly to Improve Accuracy

Use Medaka for error correction.

``` console
bash
Copy
Edit
medaka_consensus -i filtered.fastq -d assembly_output/assembly.fasta -o polished_assembly
```

### Confirm UTI Pathogen Identity via BLAST

Align assembled genomes to NCBI’s bacterial reference database.

``` console
bash
Copy
Edit
blastn -query polished_assembly/consensus.fasta -db nt -out blast_results.txt -outfmt 6
```

## Step 5: Phylogenetics & Comparative Genomics for UTI Pathogens

### Phylogenetic Placement of UTI Bacteria

Use Mashtree to analyze evolutionary relationships.

``` console
bash
Copy
Edit
mashtree --numcpus 4 polished_assembly/*.fasta > phylogeny_tree.nwk
```

### Comparative Genomics for UTI Pathogens

Use Roary or Panaroo for pan-genome analysis.

``` console
bash
Copy
Edit
roary -e -n -v *.gff
```

## Key Additions for UTI Pathogen Analysis

✅ Focus on key UTI bacteria (*E. coli, Klebsiella, Proteus*, etc.)\
✅ Detect antimicrobial resistance genes relevant to UTI treatment\
✅ Identify virulence genes involved in UTI pathogenesis\
✅ Ensure accurate pathogen strain identification

## Final Outputs

-   List of UTI pathogens present in the sample\
-   AMR profile of identified bacteria\
-   Virulence factor annotations\
-   Assembled bacterial genomes\
-   Phylogenetic relationships of UTI-causing bacteria
