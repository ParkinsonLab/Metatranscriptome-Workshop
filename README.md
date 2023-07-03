# Metatranscriptomics Practical Lab - CBW 2023

**This work is licensed under a [Creative Commons Attribution-ShareAlike 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

**This tutorial was produced by Mobolaji Adeolu (adeolum@mcmaster.ca), John Parkinson (john.parkinson@utoronto.ca) & Xuejian Xiong (xuejian@sickkids.ca)**

## Overview

This tutorial will take you through a pipeline for processing metatranscriptomic data. The pipeline, developed by the Parkinson lab, consists of various steps which are as follows:

1.  Remove adapter sequences, which are added during library preparation and sequencing steps, and trim low quality bases and sequencing reads.
2.  Remove duplicate reads to reduce processing time for following steps.
3.  Remove vector contamination (reads derived from cloning vectors, spike-ins, and primers).
4.  Remove host reads (if exploring a microbiome in which the host is an issue).
5.  Remove abundant rRNA sequences which typically dominate metatranscriptomic datasets despite the use of rRNA removal kits.
6.  Add duplicated reads, removed in step 2, back to the data set to improve quality of assemblies.
7.  Classify reads to known taxonomic groups and visualize the taxonomic composition of your dataset.
8.  Assemble the reads into contigs to improve annotation quality.
9.  Annotate reads to known genes.
10.  Map identified genes to the swiss-prot database to identify enzyme function
11.  Generate normalized expression values associated with each gene.
12. Visualize the results using KEGG metabolic pathways as scaffolds in Cytoscape.

The whole metatranscriptomic pipeline includes existing bioinformatic tools and a series of Python scripts that handle file format conversion and output parsing. We will go through these steps to illustrate the complexity of the process and the underlying tools and scripts.

New, faster, and/or more accurate tools are being developed all the time, and it is worth bearing in mind that any pipelines need to be flexible to incorporate these tools as they get adopted as standards by the community. For example, over the past two years, our lab has transitioned from cross\_match to Trimmomatic and from BLAST to DIAMOND.
Note:  This workshop was designed for use with DIAMOND v0.826.  Newer versions of DIAMOND will be incompatible with the pre-compiled database files we have made as part of this exercise.  
To illustrate the process we are going to use sequence reads generated from the contents of the colon of a mouse. These are 150 bp single-end reads. Paired-end reads can also be used, and are often preferred because they can improve annotation quality when there is enough overlap in the read pairs to improve the effective average read length. Working with paired-end data involves an additional data processing step (merging of overlapping reads) produces more files during data processing (files for merged/singleton reads, forward reads, and reverse reads), but the structure of a pipeline for paired-end data is similar to the pipeline described here and can be readily adapted.

Rather than use the entire set of 25 million read, which might take several days to process on a desktop, the tutorial will take you through processing a subset of 100,000 reads.

## Preliminaries

### Work directory

Create a new directory that will store all of the files created in this lab.

```
mkdir metatranscriptomics
cd metatranscriptomics/
```

### Python Scripts

We have written a number of scripts to extract and analyze data from the tools you will be using. Download our package for the metatranscriptomics workshop and extract our python scripts.

```
wget https://github.com/ParkinsonLab/2017-Microbiome-Workshop/releases/download/CBW2023/precomputed_files_cbw_2023.tar.gz -O precomputed_files.tar.gz
tar --wildcards -xvf precomputed_files.tar.gz *.py
```

### Input files

Our data set consists of 150 bp single-end Illumina reads generated from mouse colon contents. To inspect its contents:

```
tar -xvf precomputed_files.tar.gz mouse1.fastq
less mouse1.fastq
```

**Notes**:

-   Type `q` to exit `less`.

### Checking read quality with FastQC

```
fastqc mouse1.fastq
```

The FastQC report is generated in a HTML file, `mouse1_fastqc.html`. You'll also find a zip file which includes data files used to generate the report.

To open the HTML report file use the following command `firefox mouse1_fastqc.html` then you can go through the report and find the following information:

-   Basic Statistics: Basic information of the mouse RNA-seq data, e.g. the total number of reads, read length, GC content.
-   Per base sequence quality: An overview of the range of quality values across all bases at each position.
-   Per Base Sequence Content: A plot showing nucleotide bias across sequence length.
-   Adapter Content: Provides information on the level of adapter contamination in your sequence sample.

## Processing the Reads

### Step 1. Remove adapter sequences and trim low quality sequences. 

Trimmomatic can rapidly identify and trim adaptor sequences, as well as identify and remove low quality sequence data - It is already installed on the PCs

```
tar -xvf precomputed_files.tar.gz TruSeq3-SE.fa
ln -s TruSeq3-SE.fa Adapters
java -jar /usr/local/trimmomatic-0.36.jar SE mouse1.fastq mouse1_trim.fastq ILLUMINACLIP:Adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
```

**Notes**:

-   `ln -s TruSeq3-SE.fa Adapters` is used to create a symbolic link to the Trimmomatic supplied single-end adapter sequence files suitable for use with sequences produced by HiSeq and MiSeq machines. However, this file should be replaced with known adapter files from your own sequencing project if possible.
-   The command line parameters are:
    -   `SE`: The input data are single-end reads.
    -   `ILLUMINACLIP:Adapters:2:30:10`: remove the adapters.
    -   `LEADING:3`: Trims bases at the beginning of a read if they are below quality score of 3.
    -   `TRAILING:3`: Trims bases at the end of a read if they are below quality score of 3.
    -   `SLIDINGWINDOW:4:15`: Scan with a window of size 4 for reads with local quality below a score of 15, and trim if found.
    -   `MINLEN:50`: Delete a sequence with a length less than 50.

<!-- ***Question 1: How many low quality sequences have been removed?*** -->

Checking read quality with FastQC:

```
fastqc mouse1_trim.fastq
firefox mouse1_trim_fastqc.html
```

Compare with the previous report to see changes in the following sections:

-   Basic Statistics
-   Per base sequence quality

**Optional: Paired-end read merging**

If you were working with a paired-end dataset, we could identify pairs of sequence reads that overlap and can therefore be merged into a single sequence. For this we use the tool VSEARCH which can be found at this [website](https://github.com/torognes/vsearch):

```
Example only, do not run!
vsearch --fastq_mergepairs mouse1_trim.fastq --reverse mouse2_trim.fastq --fastqout mouse_merged_trim.fastq --fastqout_notmerged_fwd mouse1_merged_trim.fastq --fastqout_notmerged_rev mouse2_merged_trim.fastq
```

**Notes**:

-   The command line parameters are:
    -   `--fastq_mergepairs` Instructs VSEARCH to use the read merging algorithm to merge overlapping paired-end reads
    -   `--reverse` Indicates the name of the file with the 3' to 5' (reverse) paired-end reads
    -   `--fastqout` Indicates the output file contain the overlapping paired-end reads
    -   `--fastqout_notmerged_fwd` and `--fastqout_notmerged_rev` Indicates the output files containing the non-overlapping paired-end reads

If you want to see the distribution of merged read lengths you can use fastqc to examine the read properties:

```
Also example only!
fastqc mouse_merged_trim.fastq
firefox mouse_merged_trim_fastqc.html
```

**Read quality filtering**

Trimmomatic, which was used to remove the adapters and trim low quality bases in the reads, uses a sliding window method to remove contigous regions of low quality bases in reads. However, it is worthwhile to impose an overall read quality threshold to ensure that all reads being used in our analyses are of sufficiently error-free. For this we use the tool VSEARCH which can be found at this [website](https://github.com/torognes/vsearch) (when processing paired-end data, this step should come **after** the read merging step):

```
vsearch --fastq_filter mouse1_trim.fastq --fastq_maxee 2.0 --fastqout mouse1_qual.fastq
```

**Notes**:

-   The command line parameters are:
    -   `--fastq_filter ` Instructs VSEARCH to use the quality filtering algorithm to remove low quality reads
    -   `--fastq_maxee 2.0` The expected error threshold. Set at 1. Any reads with quality scores that suggest that the average expected number of errors in the read are greater than 1 will be filtered.
    -   `--fastqout` Indicates the output file contain the quality filtered reads

Checking read quality with FastQC:

```
fastqc mouse1_qual.fastq
firefox mouse1_qual_fastqc.html
```

Compare with the previous reports to see changes in the following sections:

-   Basic Statistics
-   Per base sequence quality
-   Per sequence quality

<!-- ***Question 2: How has the per read sequence quality curve changed?*** -->

### Step 2. Remove duplicate reads

To significantly reduce the amount of computating time required for identification and filtering of rRNA reads, we perform a dereplication step to remove duplicated reads using the software tool CD-HIT which can be obtained from this [website](https://github.com/weizhongli/cdhit).

```
/usr/local/cd-hit-v4.8.1-2019-0228/cd-hit-auxtools/cd-hit-dup -i mouse1_qual.fastq -o mouse1_unique.fastq
```

**Notes**:

-   The command line parameters are:
    -   `-i`: The input fasta or fastq file.
    -   `-o`: The output file containing dereplicated sequences, where a unique representative sequence is used to represent each set of sequences with multiple replicates.
-   A second output file `mouse1_unique.fastq.clstr` is created which shows exactly which replicated sequences are represented by each unique sequence in the dereplicated file and a third, empty, output file, `mouse1_unique.fastq2.clstr` is also created which is only used for paired-end reads.

<!-- ***Question 3: Can you find how many unique reads there are?*** -->

While the number of replicated reads in this small dataset is relatively low, with larger datasets, this step can reduce file size by as much as 50-80%

### Step 3. Remove vector contamination

To identify and filter reads from sources of vector, adapter, linker, and primer contamination we the Burrows Wheeler aligner (BWA) and the BLAST-like alignment tool (BLAT) to search against a database of cow sequences. As a reference database for identifying contaminating vector and adapter sequences we rely on the UniVec\_Core dataset which is a fasta file of known vectors and common sequencing adapters, linkers, and PCR Primers derived from the NCBI Univec Database. Please download it into your working directory first.

```
wget ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core
```

Now we must generate an index for these sequences for BWA and BLAT using the following commands:

```
bwa index -a bwtsw UniVec_Core
samtools faidx UniVec_Core
makeblastdb -in UniVec_Core -dbtype nucl
```

Next we can perform alignments for the reads with BWA and filter out any reads that align to our vector database with Samtools using the following commands:

```
bwa mem -t 4 UniVec_Core mouse1_unique.fastq > mouse1_univec_bwa.sam
samtools view -bS mouse1_univec_bwa.sam > mouse1_univec_bwa.bam
samtools fastq -n -F 4 -0 mouse1_univec_bwa_contaminats.fastq mouse1_univec_bwa.bam
samtools fastq -n -f 4 -0 mouse1_univec_bwa.fastq mouse1_univec_bwa.bam
```

**Notes**:

-   The commands to the following tasks:
    -   `bwa mem`: Generates alignments of reads to the vector contaminant database
    -   `samtools view`: Converts the .sam output of bwa into .bam for the following steps
    -   `samtools fastq`: Generates fastq outputs for all reads that mapped to the vector contaminant database (`-F 4`) and all reads that did not map to the vector contaminant database (`-f 4`)

<!-- ***Question 4: Can you find how many reads BWA mapped to the vector database?*** -->

Now we want to perform additional alignments for the reads with BLAT to filter out any remaining reads that align to our vector contamination database. However, BLAT only accepts fasta files so we have to convert our reads from fastq to fasta. This can be done using VSEARCH.

```
vsearch --fastq_filter mouse1_univec_bwa.fastq --fastaout mouse1_univec_bwa.fasta
```

**Notes**:

-   The VSEARCH command used, `--fastq_filter`, is the same as the command used to filter low quality reads in Step 1. However, here we give no filter criteria so all input reads are passed to the output fasta file.

Now we can use BLAT to perform additional alignments for the reads against our vector contamination database.

```
blat -noHead -minIdentity=90 -minScore=65  UniVec_Core mouse1_univec_bwa.fasta -fine -q=rna -t=dna -out=blast8 mouse1_univec.blatout
```

**Notes**:

-   The command line parameters are:
    -   `-noHead`: Suppresses .psl header (so it's just a tab-separated file).
    -   `-minIdentity`: Sets minimum sequence identity is 90%.
    -   `-minScore`: Sets minimum score is 65. This is the matches minus the mismatches minus some sort of gap penalty.
    -   `-fine`: For high-quality mRNAs.
    -   `-q`: Query type is RNA sequence.
    -   `-t`: Database type is DNA sequence.

Lastly, we can run a small python script to filter the reads that BLAT does not confidently align to any sequences from our vector contamination database.

```
pip install biopython
python 1_BLAT_Filter.py mouse1_univec_bwa.fastq mouse1_univec.blatout mouse1_univec_blat.fastq mouse1_univec_blat_contaminats.fastq
```

**Notes**:

The argument structure for this script is:
`1_BLAT_Filter.py <Input_Reads.fq> <BLAT_Output_File> <Unmapped_Reads_Output> <Mapped_Reads_Output>`

Here, BLAT does not identify any additional sequences which align to the vector contaminant database. However, we have found that BLAT is often able find alignments not identified by BWA, particularly when searching against a database consisting of whole genomes.

some alignments to vector contaminants missed by BWA in large multi-million read datasets.

### Step 4. Remove host reads

To identify and filter host reads (here, reads of mouse origin) we repeat the steps above using a database of mouse DNA sequences. For our purposes we use a [mouse genome database](ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/Mus_musculus.GRCm39.cds.all.fa.gz) downloaded from Ensembl.

```
wget ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/Mus_musculus.GRCm39.cds.all.fa.gz
gzip -d Mus_musculus.GRCm39.cds.all.fa.gz
mv Mus_musculus.GRCm39.cds.all.fa mouse_cds.fa
```

Then we repeat the steps above used to generate an index for these sequences for BWA and BLAT:

```
bwa index -a bwtsw mouse_cds.fa
samtools faidx mouse_cds.fa
makeblastdb -in mouse_cds.fa -dbtype nucl
```

Now we align and filter out any reads that align to our host sequence database using BWA and Samtools:

```
bwa mem -t 4 mouse_cds.fa mouse1_univec_blat.fastq > mouse1_mouse_bwa.sam
samtools view -bS mouse1_mouse_bwa.sam > mouse1_mouse_bwa.bam
samtools fastq -n -F 4 -0 mouse1_mouse_bwa_contaminats.fastq mouse1_mouse_bwa.bam
samtools fastq -n -f 4 -0 mouse1_mouse_bwa.fastq mouse1_mouse_bwa.bam
```

Finally, we use BLAT to perform additional alignments for the reads against our host sequence database.

```
vsearch --fastq_filter mouse1_mouse_bwa.fastq --fastaout mouse1_mouse_bwa.fasta
blat -noHead -minIdentity=90 -minScore=65  mouse_cds.fa mouse1_mouse_bwa.fasta -fine -q=rna -t=dna -out=blast8 mouse1_mouse.blatout
python 1_BLAT_Filter.py mouse1_mouse_bwa.fastq mouse1_mouse.blatout mouse1_mouse_blat.fastq mouse1_mouse_blat_contaminats.fastq
```

<!-- ***Question 5: How many reads did BWA and BLAT align to the mouse host sequence database?*** -->

***Optional:*** In your own future analyses you can choose to complete steps 3 and 4 simultaneously by combining the vector contamination database and the host sequence database using `cat UniVec_Core mouse_cds.fa > contaminants.fa`. However, doing these steps together makes it difficult to tell how much of your reads came specifically from your host organism.

### Step 5. Remove abundant rRNA sequences

rRNA genes tend to be highly expressed in all samples and must therefore be screened out to avoid lengthy downstream processing times for the assembly and annotation steps. You could use sequence similarity tools such as BWA or BLAST for this step, but we find [Infernal] (http://infernal.janelia.org/), albeit slower, is more sensitive as it relies on a database of covariance models (CMs) describing rRNA sequence profiles based on the Rfam database. Due to the reliance on CMs, Infernal, can take as much as 4 hours for ~100,000 reads on a single core. So we will skip this step and use a precomputed file, `mouse1_rRNA.infernalout`, from the tar file `precomputed_files.tar.gz`.

``` 
tar -xzf precomputed_files.tar.gz mouse1_rRNA.infernalout
```

**Notes**:

-   The commands you would use to generate this output with infernal are given below:
	 -   `vsearch --fastq_filter mouse1_mouse_blat.fastq --fastaout mouse1_mouse_blat.fasta`
    -   `cmsearch -o mouse1_rRNA.log --tblout mouse1_rRNA.infernalout --anytrunc --rfam -E 0.001 Rfam.cm mouse1_mouse_blat.fasta`
-   The command line parameters are:
    -   `-o`: the infernal output log file.
    -   `--tblout`: the simple tabular output file.
    -   `--noali`: omit the alignment section from the main output. This can greatly reduce the output volume.
    -   `--anytrunc`: relaxes the thresholds for truncated alignments
    -   `--rfam`: use a strict filtering strategy devised for large database. This will speed the search at a potential cost to sensitivity.
    -   `-E`: report target sequences with an E-value of 0.001.

From this output file we need to use a script to filter out the rRNA reads:

```
python 2_Infernal_Filter.py mouse1_mouse_blat.fastq mouse1_rRNA.infernalout mouse1_unique_mRNA.fastq mouse1_unique_rRNA.fastq
```

**Notes**:

The argument structure for this script is:
`2_Infernal_Filter.py <Input_Reads.fq> <Infernal_Output_File> <mRNA_Reads_Output> <rRNA_Reads_Output>`

Here, we only remove a few thousand reads than map to rRNA, but in some datasets rRNA may represent up to 80% of the sequenced reads.

<!-- ***Question 6: How many rRNA sequences were identified? How many reads are now remaining?*** -->


### Step 6. Rereplication

After removing contaminants, host sequences, and rRNA, we need to replace the previously removed replicate reads back in our data set.

```
python 3_Reduplicate.py mouse1_qual.fastq mouse1_unique_mRNA.fastq mouse1_unique.fastq.clstr mouse1_mRNA.fastq
```

**Notes**:

The argument structure for this script is:
`3_Reduplicate.py <Duplicated_Reference_File> <Deduplicated_File> <CDHIT_Cluster_File> <Reduplicated_Output>`

<!-- ***Question 7: How many putative mRNA sequences were identified? How many unique mRNA sequences?*** -->

Now that we have filtered vectors, adapters, linkers, primers, host sequences, and rRNA, check read quality with FastQC:

```
fastqc mouse1_mRNA.fastq
firefox mouse1_mRNA_fastqc.html
```

<!-- ***Question 8: How many total contaminant, host, and rRNA reads were filtered out?*** -->

### Step 7. Taxonomic Classification

Now that we have putative mRNA transcripts, we can begin to infer the origins of our mRNA reads. Firstly, we will attempt to use a reference based short read classifier to infer the taxonomic orgin of our reads. Here we will use [Kaiju] (https://github.com/bioinformatics-centre/kaiju) to generate taxonomic classifications for our reads based on a reference database. Kaiju can classify prokaryotic reads at speeds of millions of reads per minute using the proGenomes database on a system with less than 16GB of RAM (~13GB). Using the entire NCBI nr database as a reference takes ~43GB. Similarly fast classification tools require >100GB of RAM to classify reads against large databases. However, Kaiju still takes too much memory for the systems in the workshop so we have precompiled the classifications, `mouse1_classification.tsv`, in the tar file `precomputed_files.tar.gz`.

```
tar --wildcards -xzf precomputed_files.tar.gz kaiju*
chmod +x kaiju*
tar -xzf precomputed_files.tar.gz mouse1_classification.tsv nodes.dmp names.dmp
```

**Notes**:

-   The kaiju command you would use is given below:
    -   `kaiju -t nodes.dmp -f kaiju_db.fmi -i mouse1_mRNA.fastq -z 4 -o mouse1_classification.tsv`
-   The command line parameters are:
    -   `-t`: The hierarchal representation of the taxonomy IDs
    -   `-f`: The precomputed index for kaiju
    -   `-i`: The input reads
    -   `-z`: The number of threads supported on your system
    -   `-o`: The output file for the kaiju taxonomic classifications

We can then take the classified reads and perform supplemental analyses. Firstly, we'll restrict the specificity of the classifications to Genus-level taxa which limits the number of spurious classifications.

```
python 4_Constrain_Classification.py genus mouse1_classification.tsv nodes.dmp names.dmp mouse1_classification_genus.tsv
```

**Notes**:

The argument structure for this script is:
`4_Constrain_Classification.py <Minimum_Taxonomic_Rank> <kaiju_Classification> <nodes_file> <names_file> <Output_Classifications>`

Then we generate a human readable summary of the classification using Kaiju.

```
./kaijuReport -t nodes.dmp -n names.dmp -i mouse1_classification_genus.tsv -o mouse1_classification_Summary.txt -r genus
```

**Notes**:

-   The command line parameters are:
    -   `-t`: The hierarchal representation of the taxonomy IDs
    -   `-n`: The taxonomic names corresponding to each taxonomy ID
    -   `-i`: The kaiju taxonomic classifications
    -   `-o`: The summary report output file
    -   `-r`: The taxonomic rank for which the summary will be produced

<!-- ***Question 9: How many reads did kaiju classify?*** -->

Lastly, we will use [Krona] (https://github.com/marbl/Krona/wiki) to generate a hierarchical multi-layered pie chart summary of the taxonomic composition of our dataset.

```
./kaiju2krona -t nodes.dmp -n names.dmp -i mouse1_classification_genus.tsv -o mouse1_classification_Krona.txt
tar -xzf precomputed_files.tar.gz KronaTools
KronaTools/scripts/ImportText.pl -o mouse1_classification.html mouse1_classification_Krona.txt
```

We can then view this pie chart representation of our dataset using a web browser:

```
firefox mouse1_classification.html
```

<!-- 
***Question 10: What is the most abundant family in our dataset? What is the most abundant phylum?  
Hint: Try decreasing the `Max depth` value on the top left of the screen and/or double clicking on spcific taxa.***
-->

### Step 8. Assembling reads

Previous studies have shown that assembling reads into larger contigs significantly increases our ability to annotate them to known genes through sequence similarity searches. Here we will apply the SPAdes genome assemblers' transcript assembly algorithm to our set of putative mRNA reads.

```
/usr/local/bin/spades.py --rna -s mouse1_mRNA.fastq -o mouse1_spades
mv mouse1_spades/transcripts.fasta mouse1_contigs.fasta
```

**Notes**:

-   The command line parameters are:
    -   `--rna`: Uses the mRNA transcript assembly algorithm
    -   `-s`: The single-end input reads
    -   `-o`: The output directory
-   SPAdes assembles reads into contigs which are placed into a file named `mouse1_spades/transcripts.fasta`

<!--
***Question 11: How many assemblies did SPAdes produce?  
Hint: try using the command`tail mouse1_contigs.fasta`***
-->
In order to extract unassembled reads we need to map all putative mRNA reads to our set of assembled contigs by BWA.

First, we need to build an index to allow BWA to search against our set of contigs:

```
bwa index -a bwtsw mouse1_contigs.fasta
```

Next we attempt to map the entire set of putative mRNA reads to this contig database:

```
bwa mem -t 4 mouse1_contigs.fasta mouse1_mRNA.fastq > mouse1_contigs.sam
```

We then extract unmapped reads into a fastq format file for subsequent processing and generate a mapping table in which each contig is associated with the number of reads used to assemble that contig. This table is useful for determining how many reads map to a contig and is used for determining relative expression (see Steps 6 and 8).

```
python 5_Contig_Map.py mouse1_mRNA.fastq mouse1_contigs.sam mouse1_unassembled.fastq mouse1_contigs_map.tsv
```

**Notes**:

The argument structure for this script is:
`5_Contig_Map.py <Reads_Used_In_Alignment> <Output_SAM_From_BWA> <Output_File_For_Unassembed_Reads> <Output_File_For_Contig_Map>`

<!--
***Question 12: How many reads were not used in contig assembly? How many reads were used in contig assembly? How many contigs did we generate?***
-->
### Step 9. Annotate reads to known genes/proteins

Here we will attempt to infer the specific genes our putative mRNA reads originated from. In our pipeline we rely on a tiered set of sequence similarity searches of decreasing accuracy - BWA and DIAMOND. While BWA provides high stringency, sequence diversity that occurs at the nucleotide level results in few matches observed for these processes. Nonetheless it is quick. To avoid the problems of diversity that occur at the level of nucleotide, particularly in the absence of reference microbial genomes, we use DIAMOND searches to provide more sensitive peptide-based searches, which are less prone to sequence changes between strains.

Since BWA utilizes nucleotide searches, we rely on a [microbial genome database] (ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/all.ffn.tar.gz) that we obtained from the NCBI which contains 5231 ffn files. We then merge all 5231 ffn files into one fasta file `microbial_all_cds.fasta` and build indexes for this database to allow searching via BWA. For DIAMOND searches we use the [Non-Redundant (NR) protein database] (ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz) also from NCBI.

**Notes**:

-   The systems used in the workshop do not have enough memory to handle indexing or searching large databases like `microbial_all_cds.fasta` (9GB) and `nr` (>60GB). The descriptions in this section are purely for your information. Please use our precomputed gene, protein, and read mapping files from the tar file `tar -xzf precomputed_files.tar.gz mouse1_genes_map.tsv mouse1_genes.fasta mouse1_proteins.fasta`
-   While we only utilize BWA here, it is possible to use BWA followed BLAT for a more thorough search of the `microbial_all_cds.fasta` like we described in Steps 3 and 4. This limits the number of searches that need to be completed against the much larger NCBI nr database by Diamond/BLAST.
-   the commands used to build the indexed databases are as follows (You don't need to do these!)
    -   `bwa index -a bwtsw microbial_all_cds.fasta`
    -   `samtools faidx microbial_all_cds.fasta`
    -   `diamond makedb -p 8 --in nr -d nr`

**BWA searches against microbial genome database**

-  If you were to run BWA yourself, you would use the following commands to search the `microbial_all_cds.fasta` database:
   -  `bwa mem -t 4 microbial_all_cds.fasta mouse1_contigs.fasta > mouse1_contigs_annotation_bwa.sam`
   -  `bwa mem -t 4 microbial_all_cds.fasta mouse1_unassembled.fasta > mouse1_unassembled_annotation_bwa.sam`

Then you would run the following python script to extract high confidence alignments to the `microbial_all_cds.fasta` database and generate a read to gene mapping table. Here we are only taking one gene per contig, but it is possible that contigs may have more than one genes (e.g. co-transcribed genes).

-  `python 6_BWA_Gene_Map.py microbial_all_cds.fasta mouse1_contigs_map.tsv mouse1_genes_map.tsv mouse1_genes.fasta mouse1_contigs.fasta mouse1_contigs_annotation_bwa.sam mouse1_contigs_unmapped.fasta mouse1_unassembled.fastq mouse1_unassembled_annotation_bwa.sam mouse1_unassembled_unmapped.fasta`

The argument structure for this script is:

-  `6_BWA_Gene_Map.py <Gene_database> <Contig_Map> <Output_File_For_Gene_Map> <Output_File_For_Gene_sequences> <Contigs_File> <Contig_BWA_SAM> <Unmapped_Contigs> <Unassembled_Reads_File> <Unassembled_Reads_BWA_SAM> <Unmapped_Unassembled_Reads>`

**DIAMOND against the non-redundant (NR) protein DB**

DIAMOND is a BLAST-like local aligner for mapping translated DNA query sequences against a protein reference database (BLASTX alignment mode). The speedup over BLAST is up to 20,000 on short reads at a typical sensitivity of 90-99% relative to BLAST depending on the data and settings. However, searching time for the nr database is still long (timing scales primarily by size of reference database for small numbers of reads).

-  If you were to run DIAMOND yourself, you would use the following commands:
   -   `mkdir -p dmnd_tmp`
   -   `diamond blastx -p 4 -d nr -q mouse1_contigs_unmapped.fasta -o mouse1_contigs.dmdout -f 6 -t dmnd_tmp -k 10 --id 85 --query-cover 65 --min-score 60`
   -   `diamond blastx -p 4 -d nr -q mouse1_unassembled_unmapped.fasta -o mouse1_unassembled.diamondout -f 6 -t dmnd_tmp -k 10 --id 85 --query-cover 65 --min-score 60`
-   The command line parameters are:
    -   `-p`: Number of threads to use in the search is 4.
    -   `-q`: Input file name.
    -   `-d`: Database name.
    -   `-e`: Expectation value (E) threshold for saving hits.
    -   `-k`: Maximum number of aligned sequences to keep is 10.
    -   `-t`: Temporary folder.
    -   `-o`: Output file name.
    -   `-f`: Output file is in a tabular format.

From the output of these searches, you would need to extract the top matched proteins using the script below. Here we consider a match if 85% sequence identity over 65% of the read length - this can result in very poor e-values (E = 3!) but the matches nonetheless appear reasonable.

-  `python 7_Diamond_Protein_Map.py nr mouse1_contigs_map.tsv mouse1_genes_map.tsv mouse1_proteins.fasta mouse1_contigs_unmapped.fasta mouse1_contigs.dmdout mouse1_contigs_unannotated.fasta mouse1_unassembled_unmapped.fasta mouse1_unassembled.dmdout mouse1_unassembled_unannotated.fasta`

The argument structure for this script is:

-  `7_Diamond_Protein_Map.py <Protein_database> <Contig_Map> <Gene_Map> <Output_File_For_Protein_sequences> <Unmapped_Contigs_File> <Contig_Diamond_Output> <Output_For_Unannotated_Contigs> <Unmapped_ Unassembled_Reads_File> <Unassembled_Reads_Diamond_Output> <Output_For_Unannotated_Unassembled_Reads>`

Because the non-redundant protein database contains entries from many species, including eukaryotes, we often find that sequence reads can match multiple protein with the same score. From these multiple matches, we currently select the first (i.e. 'top hit'). As mentioned in the metagenomics lecture, more sophisticated algorithms could be applied, however our current philosophy is that proteins sharing the same sequence match are likely to possess similar functions in any event; taxonomy is a separate issue however!

<!--
***Question 13: How many reads were mapped in each step? How many genes were the reads mapped to? How many proteins were the genes mapped to?***
-->
-   Total number of mapped-reads with BWA = 3356 reads
-   Total number of mapped genes (BWA) = 1234
-   Total number of mapped-reads with DIAMOND = 51936 reads
-   Total number of mapped proteins (DIAMOND) = 21699

Thus of ~83000 reads of putative microbial mRNA origin, we can annotate ~55000 of them to almost ~23000 genes!

Remember, to extract the precomputed output files for this step:

```
Run this one!  
tar -xzf precomputed_files.tar.gz mouse1_genes_map.tsv mouse1_genes.fasta mouse1_proteins.fasta
```

### Step 10. Enzyme Function Annotation

To help interpret our metatranscriptomic datasets from a functional perspective, we rely on mapping our data to functional networks such as metabolic pathways and maps of protein complexes. Here we will use the KEGG carbohydrate metabolism pathway.

To begin, we need to first match our annotated genes the enzymes in the KEGG pathway. To do this, we will use Diamond to identify homologs of our genes/proteins from the SWISS-PROT database that have assigned enzyme functions. Diamond is a relatively coarse and straight forward way to annotate enzyme function by homology. We have chosen to use it here in order to avoid having to introduce additional tools. However, more robust methods for enzymatic function annotation exist in literature, such as our own probability density based enzyme function annotation tool, [Detect] (https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btq266).

```
mkdir -p dmnd_tmp
tar -xzf precomputed_files.tar.gz swiss_db.dmnd swiss_map.tsv
```

For microbial *genes* identified through our BWA searches:

```
diamond blastx -p 4 -d swiss_db -q mouse1_genes.fasta -o mouse1_genes.diamondout -f 6 -t dmnd_tmp -e 10 -k 1
```

For *proteins* identified through our DIAMOND searches:

```
diamond blastp -p 4 -d swiss_db -q mouse1_proteins.fasta -o mouse1_proteins.diamondout -f 6 -t dmnd_tmp -e 10 -k 1
```

We then need to generate a mapping file which lists our gene/protein and the enzyme commission (EC) number, describing enzymatic function, which corresponds to it:

```
python 8_Gene_EC_Map.py swiss_map.tsv mouse1_genes.diamondout mouse1_proteins.diamondout mouse1_EC_map.tsv
```

The argument structure for this script is:

`8_Gene_EC_Map.py <SWISS-PROT_EC_Mappings> <Diamond_Output_For_Genes> <Diamond_Output_For_Proteins> <Output_EC_Mapping_File>`

<!--
***Question 14: How many unique enzyme functions were identified in our dataset?***
-->
### Step 11. Generate normalized expression values associated with each gene

We have removed low quality bases/reads, vectors, adapters, linkers, primers, host sequences, and rRNA sequences and annotated reads to the best of our ability - now lets summarize our findings. We do this by looking at the relative expression of each of our genes in our microbiome.

```
python 9_RPKM.py nodes.dmp mouse1_classification.tsv mouse1_genes_map.tsv mouse1_EC_map.tsv mouse1_RPKM.txt mouse1_cytoscape.txt
```

**Notes**:

-   The argument structure for this script is:

    -   `9_RPKM.py <nodes_file> <kaiju_Classification> <Gene_Mapping_File> <EC_Mapping_File> <RPKM_Output> <Cytoscope_Attributes_Output>`

-   The structure of the output file `mouse1_RPKM.txt` is:
    -   `[geneID/proteinID, length, #reads, EC#, Total RPKM, RPKM per phylum]`
    -   `gi|110832861|ref|NC_008260.1|:414014-415204 1191 1 3.9.3.5 106.98 0 0 45.89 6.86 20.77 7.35 2.3 0 4.63`

<!--
***Question 15: have a look at the `mouse1_RPKM.txt` file. What are the most highly expressed genes? Which phylum appears most active?***
-->
### Step 12. Visualize the results using a KEGG Pathway as a scaffold in Cytoscape.

To visualize our processed microbiome dataset in the context of the carbohydrate metabolism pathways, we use the network visualization tool - Cytoscape together with the enhancedGraphics and KEGGscape plugins. Some useful commands for loading in networks, node attributes and changing visual properties are provided below (there are many cytoscape tutorials available online).


**Download the metabolic pathway**

First, download the carbohydrate metabolism pathways from KEGG using the following commands:

```
wget https://github.com/ParkinsonLab/2017-Microbiome-Workshop/releases/download/EC/ec00010.xml
wget https://github.com/ParkinsonLab/2017-Microbiome-Workshop/releases/download/EC/ec00500.xml
```

You can find other [pathways on KEGG] (http://www.genome.jp/kegg-bin/get_htext?htext=br08901.keg) which can also be imported into Cytoscape by selecting the `Download KGML` option on the top of the page for each pathway.

**Install the Cytoscape plugins**

-   Select `Apps` -> `App Manager`
-   Search for `enhancedGraphics`
-   Select `enhancedGraphics` in the middle column then click `Install` in the bottom right
-   Search for `KEGGScape`
-   Select `KEGGScape` in the middle column then click `Install` in the bottom right

**Import an XML from KEGG into Cytoscape**

-   Select `File` -> `Import` -> `Network` -> `File...`
-   Select the XML file, `ec00010.xml` or `ec00500.xml` and click `Open`
-   Check `Import pathway details from KEGG Database` box then select `OK`

**Loading a node attribute text file (.txt) - this will map attributes to nodes in your network which you can subsequently visualize**

-   Select `File` -> `Import` -> `Table` -> `File...`
-   Select the `mouse1_cytoscape.txt` file and click `Open`
-   Change the `Key Column for network` from `shared name` to `KEGG_NODE_LABEL`
-   Click OK

**Visualizing your node attributes**

-   In the left `Control Panel` select the `Style` tab
-   Check the `Lock node width and height` box
-   Click the left-most box by the `Size` panel and change the default node size to 20.0
-   Click the blank box immediately to the right of the box you clicked to change the default size, change the `Column` field to `RPKM` and the `Mapping Type` field to `Continuous Mapping`
-   Click the left-most box by the `Image/Chart 1` panel, switch to the `Charts` tab, Click the doughnut ring icon, and press the `>>` "add all" button between the two column fields before clicking apply (make sure to remove overall RPKM from the fields that are added to the doughnut ring)
-   If you do not see the `Image/Chart 1` panel, select `Properties` -> `Paint` -> `Custom Paint 1` -> `Image/Chart 1` from the to left corner of the control panel
-   To improve the visualization you can modify colour properties under `Image/Chart 1` -> `Charts` -> `Options`, or modify other properties such as Label Font Size, Label Position, Fill Color, Node location, and edge properties

**Notes**:

-   A cytoscape file with node attributes precalculated is provided for your convenience, `tar -xzf precomputed_files.tar.gz Example.cys`, feel free to open it and play with different visualizations and different layouts - compare the circular layouts with the spring embedded layouts for example. If you want to go back to the original layout then you will have to reload the file.
-   Cytoscape can be temperamental. If you don't see pie charts for the nodes, they appear as blank circles, you can show these manually. Under the 'properties' panel on the left, there is an entry labeled 'Custom Graphics 1'. Double click the empty box on the left (this is for default behavior) - this will pop up a new window with a choice of 'Images' 'Charts' and 'Gradients' - select 'Charts', choose the chart type you want (pie chart or donut for example) and select the different bacterial taxa by moving them from "Available Columns" to "Selected Columns". Finally click on 'Apply' in bottom right of window (may not be visible until you move the window).

**Visualization Questions:**

- Which genes are most highly expressed in these two systems?
- Which taxa are responsible for most gene expression?
- Can you identify sub-systems (groups of interacting genes) that display anomalous taxonomic profiles?
- Think about how you might interpret these findings; for example are certain taxa responsible for a specific set of genes that operate together to fulfill a key function?
- Can you use the gene annotations to identify the functions of these genes through online searches?
- Think about the implications of sequence homology searches, what may be some caveats associated with interpreting these datasets?
