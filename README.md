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

We will first define three variables, `${POD5DIR}` is the directory
that contains the `.pod5` files containing the Oxford Nanopore
sequence data, `${KITNAME}` is the barcoding kit used to prepare
the libraries, and `${OUTPUT}` will be used to name output files
and directories for the basecalls and demultiplexed data.

**These will differ between projects, and the examples below should
be modified to reflect your own data.**

```{bash}
POD5DIR=pod5_pass
KITNAME=SQK-NBD114-24
OUTPUT=mepa_pathogens
```

The following command will use the "super accurate" basecalling
model v5.2.0

```{bash}
dorado basecaller sup@v5.2.0 pod5_pass/ --kit-name ${KITNAME} --min-qscore 10 --recursive > ${OUTPUT}_v5.2.0.bam
```

This process can be run substantially faster on a computer with GPU
resources, in which case add the parameter `-x cuda:all` to the above
command before `> ${OUTPUT}_v5.2.0.bam`.

For comparison, basecalling on these data completed in 2.5 hours using
GPUs on a high-performance computing cluster with 128G of memory and 10
CPUs.

### Demultiplexing

The basecalled BAM file contains data from all samples, this command
will demultiplex the base calls into FASTQ files separately by barcode
(note: the `--no-classify` parameter is used because the barcodes were
identified in the previous step, and trimmed)

```{bash}
dorado demux ${OUTPUT}_v5.2.0.bam -o ${OUTPUT}_v5.2.0 --emit-fastq --emit-summary --no-classify
```

### Host Sequence Filtering

Even though our data was generated using host-depleted sequencing,
the host genome can still make it through. Because of this, it is
necessary to remove these sequences before proceeding. We use the
list of accession numbers from the sample sheet to download the
reference genomes for read alignment.

1. Prepare list of accession numbers to download

```{bash}
cut -f 4 -d ',' sample_sheet.csv  | tail -n +2 | sort | uniq > ref_accessions.txt
```

2. Download the reference genomes from NCBI

```{bash}
datasets download genome accession --inputfile ref_accessions.txt
```

3. Decompress the reference data files

```{bash}
unzip ncbi_data.zip
```

4. Filter the read data [following a method described at this link](https://linsalrob.github.io/ComputationalGenomicsManual/Deconseq/).

First create a new directory and then create some arrays to store the values for `BARCODE` and `ACCESSION` from `sample_sheet.csv`, and find the corresponding `FASTQ` and `REFERENCE` files n the reads for each barcode to the host genome listed in `sample_sheet.csv`.   

```{bash}
mkdir -p reads/qc
index=0
while read -r LINE;
    do BARCODE[index]=$(echo $LINE | cut -f 1 -d ',');
    FASTQ[index]=$(find . -name "*_${BARCODE[index]}_*.fastq");
    ACCESSION[index]=$(echo $LINE | cut -f 4 -d ',' | tr -d '\r\n');
    REFERENCE[index]=$(find . -name "${ACCESSION[index]}*.fna");
    let "index++"
done < sample_sheet.csv
```

Next, iterate through these arrays and map the reads to the host genome

```{bash}
for i in ${!BARCODE[@]};
    do minimap2 --split-prefix=tmp$$ -a -xsr ${REFERENCE[i]} ${FASTQ[i]} | samtools view -bh | samtools sort -o reads/qc/${BARCODE[i]}_host_aligned.bam
done
```
Following alignment we use the [samtools flags](https://broadinstitute.github.io/picard/explain-flags.html) to filter the reads that mapped to the host (`reads/qc/${BARCODE}_host.fastq`) and those that did 
not map to the host (`reads/${BARCODE}_nonhost.fastq`).

```{bash}
for i in ${!BARCODE[@]};
    do samtools fastq -F 3588 reads/qc/${BARCODE[i]}_host_aligned.bam > reads/qc/${BARCODE[i]}_host.fastq
    samtools fastq -F 3584 -f 4 reads/qc/${BARCODE[i]}_host_aligned.bam > reads/${BARCODE[i]}_nonhost.fastq
done
```

### Read classification using kraken2

[kraken2](https://github.com/DerrickWood/kraken2) is a program that provides taxonomic
classification for FASTQ read data. It requires a database in order to perform this classification,
and there are a [large number of options available here](https://benlangmead.github.io/aws-indexes/k2).

For the purposes of this tutorial we will download the Standard database, which includes sequences
for archaea, bacteria, viruses, plasmids, and human. The full database is quite large, [but there is
a smaller version (Standard-8)](https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260226.tar.gz) that is capped at 8 GB and should be sufficient for our purposes.


Download and extract database files

```{bash}
mkdir -p kraken_db

wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260226.tar.gz -O kraken_db/k2_standard_08_GB_20260226.tar.gz

tar -C kraken_db/ -xvf kraken_db/k2_standard_08_GB_20260226.tar.gz
```

Now iterate through the FASTQ files and classify them with kraken2

```{bash}
mkdir -p kraken_results

for i in ${!BARCODE[@]};
    do k2 classify --db kraken_db --report kraken_results/${BARCODE[i]}_kraken2_report.txt --output kraken_results/${BARCODE[i]}_kraken2_classified.txt --use-names reads/${BARCODE[i]}_nonhost.fastq
done
```
