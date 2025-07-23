# Heamoproteus majoris lineage WW2: Pipeline of the assembly and annotations

## ASSEMBLY
### STEP 1: Assembly of Delurb4 with FLYE (raw data, fastq) and decontamination with Centrifuge/Recentrifuge
#### Filtering with NanoFilt (v2.8.0)
gunzip -c nanopore_long_reads_file.fastq.gz | NanoFilt -l 300 --readtype 1D | gzip > NANOFILT_file.fastq.gz
#### Assembly with FLYE (v2.8.3)
flye --nano-raw NANOFILT_file.fastq.gz --out-dir $FLYE_folder -t 30
#### Decontamination/taxonomic classification with centrifuge/recentrifuge. NOTE: The decontamination has been done at this stage due to large data. It was not possible to run PILON.
#### Download all the bank necessary to the decontamination phase. I had to add folder/file path for the decontamination part. Command lines of centrifuge/recentrifuge from ILRA software. Options are: Same options used in ILRA software.
#### Centrifuge (v1.0.4_beta). Command line from ILRA software.
centrifuge -f -x .../databases/nt -U FLYE_assembly.fasta -p 30 --report-file report.txt -S classification.txt --min-hitlen 100 1> centrifuge_log_out.txt 2> centrifuge_log_out_warnings_errors.txt
cat report.txt
#### Recentrifuge (v1.3.2) 
#### Extract contigs classified as different organisms. Command line from ILRA software. rcf command from RECENTRIFUGE. Perl script: fasta_to_fastq.pl from ILRA. 
rcf -n .../databases/taxdump -f classification.txt -o recentrifuge_contamination_report.html -e CSV &> /rcf_log_out.txt
perl -S fasta_to_fastq.pl FLYE_assembly.fasta ? > FLYEdecontaminated_assembly.fa.fq
#### Excluded contigs. rextract command from RECENTRIFUGE. Command line from ILRA software.
rextract -f classification.txt -i 77521 -i 5820 -x 2 -x 10239 -x 40674 -x 81077 -n .../databases/taxdump -q FLYEdecontaminated_assembly.fa.fq &> rextract_log_out.txt
sed -n '1~4s/^@/>/p;2~4p' *.fastq > FLYEdecontaminated_assembly.fasta
#### Save contigs.Command line from ILRA software.
echo -e "\n" >> ../ExcludedWW2.contigs.fofn; echo -e "
comm -23 <(cat FLYE_assembly.fasta | grep ">" | sort) <(cat FLYEdecontaminated_assembly.fasta | grep ">" | sort) >> ../ExcludedWW2.contigs.fofn
#### WARNING: CENTRIFUGE and RECENTRIFUGE command lines from ILRA software remove sequences defined as unclassified. The unclassified sequences are trimmed by hand if there are(blast). The sequences not considered as contamination are added back for the next step.

###  STEP 2: Correction of the assembly with racon and medaka to improve the consensus accuracy of an assembly of ONT sequencing reads 
#### First correction with Racon: BWA (v0.7.17) and Racon (v1.4.20)
bwa index FLYEkeepunclassified_assembly.fasta
bwa mem -t 30 -x ont2d FLYEkeepunclassified_assembly.fasta NANOFILT_file.fastq.gz > FLYEkeepunclassified_assembly.sam
racon -t 30 NANOFILT_file.fastq.gz FLYEkeepunclassified_assembly.sam FLYEkeepunclassified_assembly.fasta > FLYE_racon.fasta
#### Second correction with Medaka (v1.6.0)
medaka_consensus -i NANOFILT_file.fastq.gz -d FLYE_racon.fasta -o medaka_folder -t 30 -m r941_min_high_g360

### STEP 3: Postprocessing of the assemblies with ILRA (only step 1, 2, 3). It is a version not published in github https://github.com/ThomasDOtto/ILRA. The version 1 of the article have been posted in biorxiv in August 1, 2021-19:01; https://www.biorxiv.org/content/10.1101/2021.07.30.454413v3.article-info
cd /home/amelie/raid8/ILRA
source /home/amelie/raid8/ILRA/external_software/ILRA/path_to_source
bash ILRA.sh -a medaka.fasta -o Output_folder -c NANOFILT_file.fastq.gz -n name_of_folder -r reference.fasta -I illumina_shortreads_folder -t 30 -T 77521 -g reference.gff -L ont &> Output_folder/output_ILRA.txt

### STEP 4: Subsampling short reads and Correction with Pilon. There is to much reads to run pilon on our cluster. we add to downsampling our high coverage.
#### PILON round1
#### Mapping short reads with bowtie2, transform the sam in bam and then sorting/indexing the bam with samtools and use Pilon to correct the assembly.For the first step we downsampling the data to a mean of 1000 reads per position. 
#### Align the short reads
bowtie2-build --threads 30 output_step3ILRA_folder output_step3ILRA.fasta
bowtie2 --threads 30 -x output_step3ILRA.fasta -1 illumina_shortreads_1.fastq.gz -2 illumina_shortreads_2.fastq.gz -S output_step3ILRA.sam
#### Transform the output .sam to a .bam
samtools view -S -b output_step3ILRA.sam -@ 30 > output_step3ILRA.bam
samtools sort output_step3ILRA.bam -@ 30 -o output_step3ILRAsorted.bam
samtools index output_step3ILRAsorted.bam
#### variant bam for downsampling. For heterogenous genome, treshold: above 1000 reads 
variant $output_step3ILRAsorted.bam -m 1000 -o output_step3ILRAsortedsubsamp.bam -b
#### Sort the .bam
samtools sort output_step3ILRAsortedsubsamp.bam -@ 30 -o output_step3ILRAsortedsubsampsorted.bam
samtools index output_step3ILRAsortedsubsampsorted.bam
#### pilon_command line
pilon --genome output_step3ILRA.fasta --frags output_step3ILRAsortedsubsampsorted.bam --output Pilon_round1_folder --changes --fix all --threads 30 --verbose | tee Pilon_round1.pilon
#### PILON round2
bowtie2-build --threads 30 Pilon_round1_folder Pilon_round1.fasta
bowtie2 --threads 30 -x Pilon_round1.fasta -1 illumina_shortreads_1.fastq.gz -2 illumina_shortreads_2.fastq.gz -S Pilon_round1.sam
samtools view -S -b Pilon_round1.sam -@ 30 > Pilon_round1.bam
samtools sort Pilon_round1.bam -@ 30 -o ILRA_round1_sorted.bam
samtools index Pilon_round1_sorted.bam
variant Pilon_round1_sorted.bam -m 1000 -o Pilon_round1_sortedsubsamp.bam -b
samtools sort Pilon_round1_sortedsubsamp.bam -@ 30 -o Pilon_round1_sortedsubsampsorted.bam
samtools index Pilon_round1_sortedsubsampsorted.bam
pilon --genome Pilon_round1.fasta --frags Pilon_round1_sortedsubsampsorted.bam --output Pilon_round2_folder --changes --fix all --threads 30 --verbose | tee Pilon_round2.pilon
#### PILON round3
bowtie2-build --threads 30 Pilon_round2_folder Pilon_round2.fasta
bowtie2 --threads 30 -x Pilon_round2.fasta -1 illumina_shortreads_1.fastq.gz -2 illumina_shortreads_2.fastq.gz -S Pilon_round2.sam
samtools view -S -b Pilon_round2.sam -@ 30 > Pilon_round2.bam
samtools sort Pilon_round2.bam -@ 30 -o Pilon_round2_sorted.bam
samtools index Pilon_round2_sorted.bam
variant Pilon_round2_sorted.bam -m 1000 -o Pilon_round2_sortedsubsamp.bam -b
samtools sort Pilon_round2_sortedsubsamp.bam -@ 30 -o Pilon_round2_sortedsubsampsorted.bam
samtools index Pilon_round2_sortedsubsampsorted.bam
pilon --genome Pilon_round2.fasta --frags Pilon_round2_sortedsubsampsorted.bam --output Pilon_round3_folder --changes --fix all --threads 30 --verbose | tee Pilon_round3.pilon

### STEP 5: Postprocessing of the assemblies with ILRA only step 5 and 7. It is a version not published in github https://github.com/ThomasDOtto/ILRA. The version 1 of the article have been posted in biorxiv in August 1, 2021-19:01; https://www.biorxiv.org/content/10.1101/2021.07.30.454413v3.article-info
cd /home/amelie/raid8/ILRA
source /home/amelie/raid8/ILRA/external_software/ILRA/path_to_source
#### ILRA command line
bash ILRA.sh -a Pilon_round3.fasta -o Output_folder_Step5 -c NANOFILT_file.fastq.gz -n name_of_folder_Step5 -r reference.fasta -I illumina_shortreads_folder -t 30  -T 77521 -g reference.gff -L ont -d step5 &> Output_folder/output_ILRA_Step5.txt

### STEP6: assembly quality
#### QUAST
#### BUSCO

## ANNOTATION

### STEP 1: Gene annotation using the online version of Companion (v2.2.11). We used P. relictum SGS1-like as reference. The reference is inside Companion (Release 63). 
### Option used: Use of reference protein evidence, Pseudogene detection, Structural annotation using AUGUSTUS with a score of 0.2 and the RATT transfer tool with the Assembly transfert type, Taxon ID: 77521. OrthoFinder based on Diamond hits (e-value cutoff: 10-3). 
https://companion.gla.ac.uk/

### STEP 2: TE annotation using EDTA (v2.2.2).
### Option used: --force 1, --sensitive 1, --anno 1 and --evaluate 1
perl /home/amelie/EDTA/EDTA_raw.pl --genome file.fasta --type ltr --threads 30
perl /home/amelie/EDTA/EDTA.pl --genome file.fasta --overwrite 1 --force 1 --sensitive 1 --anno 1 --evaluate 1 --threads 30

### STEP 3: Gene re-annotation using the online version of Companion (v2.2.11). We used P. relictum SGS1-like as reference. The reference is inside Companion (Release 63). 
### Option used: Use of reference protein evidence, Pseudogene detection, Structural annotation using AUGUSTUS with a score of 0.2 and the RATT transfer tool with the Assembly transfert type, Taxon ID: 77521. OrthoFinder based on Diamond hits (e-value cutoff: 10-3). 
https://companion.gla.ac.uk/

#### BUSCO protein
