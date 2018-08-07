# ExportBowtie
export Bowtie=/hpc/grid/shared/ngsapp/bowtie2-2.2.4
export PATH=$Bowtie:$PATH
Indexing a reference genome
To create an index for the Lambda phage reference genome included with Bowtie 2, create a new temporary directory (it doesn't matter where), change into that directory, and run:
$BT2_HOME/bowtie2-build $BT2_HOME/example/reference/lambda_virus.fa lambda_virus

bowtie2-build GAPDH.fa GAPDH
Stay in the directory created in the previous step, which now contains the lambda_virus index files. Next, run:

$ bowtie2 -x GAPDH -U reads_1.fq -S eg1.sam
bowtie2 -x GAPDH -U ../../../fastq/170707_NS500482_0311_AHCM7NBGX2_Andre_Franco/A2pos_S36_2.fastq.gz -S eg1.sam

This runs the Bowtie 2 aligner, which aligns a set of unpaired reads to the Lambda phage reference genome using the index generated in the previous step. The alignment results in SAM format are written to the file eg1.sam, and a short alignment summary is written to the console. (Actually, the summary is written to the "standard error" or "stderr" filehandle, which is typically printed to the console.)
To see the first few lines of the SAM output, run:
head eg1.sam


[sjelinsk@amrndhl1315]$ bowtie2 -x GAPDH -U ../../../fastq/170707_NS500482_0311_AHCM7NBGX2_Andre_Franco/A11pos_S59_2.fastq.gz -S eg1.sam
#index reference genome
Samtools fa

####Sort and index and convert to BAM file
samtools import eg1.fai eg1.sam eg1.bam

2.  Run bowtie2-build command
ml GCCcore/.5.4.0
ml Bowtie2
bowtie2-build FN1.fa FN1

3. To use bowtie2 to map this data, run the following command:

/path-to-bowtie-programs/bowtie2 -p <# cpu> -x <genome index prefix> <fastq file>  > <output filename>
i.e.
/programs/bowtie2 -p 8 -x hg19 Experiment1.fastq > Experiment1.sam

bowtie2 -x FN1 -U /hpc/grid/iandi/AMP/AMP_RA_SLE.Phase1/RA/RNASeq_lowinputSS2_Broad/fastq/S45.2.fastq.gz -S Exp1a.sam   
to mapped paired samples
bowtie2 -x FN1 -1 /hpc/grid/iandi/AMP/AMP_RA_SLE.Phase1/RA/RNASeq_lowinputSS2_Broad/fastq/S45.1.fastq.gz -2 /hpc/grid/iandi/AMP/AMP_RA_SLE.Phase1/RA/RNASeq_lowinputSS2_Broad/fastq/S45.2.fastq.gz -S Exp1paired.sam 

test filtered sam
samtools  view -q 50 -f 0x2 unfiltered.sam > filtered_by_samtools.sam

3.  Convert to BAM
samtools view -b -S file.sam > file.bam
samtools view –b –S Exp1a.sam > Exp1a.bam
4.  Sort a BAM file
samtools sort -m 1000000000 file.bam outputPrefix
samtools sort –m 1000000000 Exp1a.bam Exp1a.bai
	samtools index Exp1a.bai.bam

3.5 Filter to remove non mapped sequences (-F 4) and low quality (-q 42)
filter bam files that map to genome and high high quality read score
samtools view -F 4 -q 42 Exp1a.bam> Exp1a.bam.f

5 View in IGV
ml IGV
 igv.sh
 

Load custom Genome 
Load 



Converting BAM to SAM and vice versa
A common thing you may want to do is view the contents of a BAM file.  Since BAM files are a binary format, it is hard to read them with a text editor.  To convert a bam file to a sam file:
samtools view -h file.bam > file.sam
To convert back to a bam file:
samtools view -b -S file.sam > file.bam
Sorting a BAM file
Many of the downstream analysis programs that use BAM files actually require a sorted BAM file.  This allows access to reads to be done more efficiently.  To sort a BAM file:
samtools sort -m 1000000000 file.bam outputPrefix
The number of -m specifies the maximum memory to use, and can be changed to fit your system.  The output will be placed in file named "outputPrefix.bam".
Creating a BAM index
Some tools require a BAM Index file to more efficiently access reads in a BAM file.  An example of this is viewing alignments using the UCSC Genome Browser.  To create a BAM index, you must first sort the BAM file.  After that, run the following command to make a BAM index:
samtools index sorted.bam
This will create a file named sorted.bam.bai which contains the index.



  940  ml RHEL6-apps
  941  ml bowtie
  942  bowtie2-build
  943  ml bowtie2
  944  module spider bowtie2
  945  ml Bowtie2
  946  ml GCCcore/.5.4.0
  947  ml Bowtie2
  948  bowtie2-build
  949  ls
  950  bowtie2-build FN1.fa FN1
  951  ls
  952  ls -l
  953  bo FN1
  954  bowtie2 -build FN1.fa FN1
  955  bowtie2 -x FN1 -U
  956  ls
  957  cd ..
  958  ls
  959  cd Reference/
  960  ls
  961  bowtie2 -x FN1 -U /hpc/grid/iandi/AMP/AMP_RA_SLE.Phase1/RA/RNASeq_lowinputSS2_Broad/fastq/S45.1.fastq.gz -S eg1.sam
  962  ls
  963  head eg1.sam
  964  igv
  965  IGV
  966  bowtie2 -x FN1 -U   -S eg2.sam
  967  head eg2.sam
  968  ml igv
  969  module spider igv
  970  ml IGV
  971  run igv.sh
  972  igv.sh
  973  ls
  974  samtools view -bh -o eg2.bam eg2.sam
  975  samtools view -bh -o eg1.bam eg1.sam
  976  samtools index eg2.sam eg2.sort
  977  ls
  978  samtools import eg2.sam eg2.bam
  979  ls
  980  samtools import eg2.fai eg2.sam eg2.bam
  981  samtools faidx FN1.fa
  982  ls
  983  samtools import FN1.fa.fai eg2.sam eg2.bam
  984  ls
  985  samtools faidx eg2.sam
  986  samtools faidx eg2.bam
  987  ls -l
  988  samtools index eg2.bam
  989  samtools sort eg2.bam
  990  samtools sort eg2.bam .bai
  991  ls -l
  992  ls -al
  993  mv .bai.bam eg2.bai.bam
  994  less eg2.sam
  995  ls
  996  samtools sort -m eg2.bam egsort
  997  samtools sort eg2.bam egsort
  998  samtools index egsort.bam
  999  ls -l
 1000  more egsort.bam.bai
 1001  ls
 1002  history		
