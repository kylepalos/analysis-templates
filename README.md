# Repo to host code, scripts, and data for [DOI here]

## The HAMR algorithm was run using [HAMRLINC](https://github.com/harrlol/HAMRLINC).

## ModTect was run using raw bam files from STAR mapping. 

### Make a STAR index

```
STAR --runMode genomeGenerate \
--genomeDir genome \ # prefix for index files
--genomeFastaFiles reference_genome.fa \ 
--sjdbGTFfile genome_ann.gff3 \
--sjdbGTFtagExonParentTranscript Parent # GFF3 specific option

```

### Trim reads using [fastp](https://github.com/OpenGene/fastp) in automatic mode

```
for i in *_1.fastq.gz
> do
> name=$(basename ${i} _1.fastq.gz);
> fastp -i ${name}_1.fastq.gz -I ${name}_2.fastq.gz \
> -o trimmed/${name}_1.fq.gz -O trimmed/${name}_2.fq.gz \
> -w 8 -h ${name}.html -j ${name}.json
> done
```

### Map reads

```
for i in *_1.fq.gz
> do
> name=$(basename ${i} _1.fq.gz);
> STAR --runMode alignReads --readFilesIn ${name}_1.fq.gz ${name}_2.fq.gz \
--runThreadN 16 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${name} \
--alignIntronMax 10000 --outFilterIntronMotifs RemoveNoncanonical \ # 10KB maximum intron length
--genomeDir /home/kpalos/sequencing_reads/trimmed/sb_star_index/ --outSAMstrandField intronMotif
> done
```

### Run [ModTect](https://github.com/ktan8/ModTect)


```
for i in *.bam
> do
> name=$(basename ${i} Aligned.sortedByCoord.out.bam); # STAR aligned BAM file
> python2 /home/kpalos/ModTect/ModTect_1_7_5_1.py \
${name}Aligned.sortedByCoord.out.bam \
reference_genome.fa 1 1 2 \ # the "1 1 2" is a ModTect place-holder to specify chromosome, start and end to search for mods - we're actually using a bed file
--regionFile genes.bed --readlength ${} \ # insert the correct read length
--minBaseQual 30 --label ${name} # SAM/BAM minimum alignment quality threshold
> done
```

### Intersect ModTect output with gene bed file
To make the ModTect output work with [BedTools](https://bedtools.readthedocs.io/en/latest/), you need to duplicate the second column so it's in bed format, e.g.:

**chromosome** **start** **end**

**I used the folloinwg bash code to duplicate that column:**

```
for i in *.txt # each ModTect file is a .txt file
do
name=$(basename ${i} .combined.txt);
awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' ${name}.combined.txt > ${name}.txt # duplicate the second column - enforce tab delimited format
done
```

**Bedtools intersect:**

```
for i in *.txt
do
name=$(basename ${i} .txt);
bedtools intersect -wao -a ${name}.txt -b genes.bed > ${name}_intersected.txt
done
```

**I'm going to combine all the PRM files together, add the filename to each file so experiment can be linked to PRMs:**

```
for i in *.txt
do
name=$(basename ${i} _intersected.txt);
awk 'BEGIN{OFS="\t"}{print $0, FILENAME}' ${name}_intersected.txt > ${name}.txt
done
```




This repository is meant to host analysis templates for the various research projects that I have worked on.

I will try to provide the necessary input data where I can, but some projects may be too sensitive/far from publishing that I cannot make the data freely available.

Additionally, as projects become published, I will post all the necessary analyses here.

## Model RNA mods XGBoost

I work on many RNA modification projects. RNA modidifications can be predicted from regular RNA-seq with algorithms like [HAMR](https://github.com/GregoryLab/HAMR) and [ModTect](https://github.com/ktan8/ModTect). 

ModTect detects substantially more modifications relative to HAMR, but HAMR was build to classify RNA modifications using a reverse transcriptase error profile modeled on Yeast tRNA modifications. I wanted to the scale of RNA modifications predicted by ModTect but have the RNA modification class that HAMR reports. This script uses XGBoost to predict ModTect modification classes from HAMR output. Briefly, the important bits of data to model are the reference nucleotide of the position, the distribution of observed nucleotides at that position, and the relative mutation ratio of the position. Input data is present in this repo. The analysis HTML can be accessed in this repo.
