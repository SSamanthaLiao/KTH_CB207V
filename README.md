1. Download Fastq files from the ENCODE and NCBI GEO datasets by SRAtoolkits.(Linux) \
   JunD: `fasterq-dump SRR502543 && mv SRR502543 jund.fastq` \
   Input: `fasterq-dump SRR650337 && mv SRR650337 input.fastq` \
   IgG: `fasterq-dump SRR650333  && mv SRR650333 IgG.fastq`
2. Fastqc quality check (Linux)
  `for i in jund input IgG; do fastqc ${i}.fastq; done`
3. Mapping data to the reference genome(Linux)
   ```bash
   bowtie2-build chm13v2.0.fa chm13v2_bt2
   for i in jund input IgG; do
     bowtie2 -x chm13v2_bt2 -U ${i}.fastq | awk '$2 != 4 {print}' | samTools view -S -b '-' > ${i}.bam;
   done
   ```
4. Peak calling (Linux)
   ```bash
   for i in jund input IgG; do
     samtools sort ${i}.bam -o ${i}_sorted.bam 
     samtools index ${i}_sorted.bam;
   done
   
   macs3 callpeak -t jund.bam -c input.bam -f BAM -g hs -B --call-summits -n jund
   macs3 callpeak -t IgG.bam -c input.bam -f BAM -g hs -B --call-summits -n IgG
   ```
5. Visualising reads and peaks in IGV (haven't done yet) (Linux)
   ```bash
   for bdg in *_treat_pileup.bdg do
    name=$(echo $bdg | awk -F"/" '{print $NF}' | awk -F"_treat_pileup" '{print $1}') echo $name
    bedGraphToBigWig ${name}_treat_pileup.bdg CHM13v2chromSizes.txt ${name}.bigWig done
   bedGraphToBigWig HS30_HSF1_control_lambda.bdg CHM13v2chromSizes.txt Input.bigWig
   ```
6. Removing unspecific peaks (Linux)
   ```bash
   bedtools intersect -v -a jund_peaks.narrowPeak -b igg_peaks.narrowPeak > jund.narrowpeak
   ```
7. Identifying binding sites across the genome (R)
   ```R
   ## manage the referece genome files
   library(dplyr)
   names = c("chr","peakStart","peakEnd","name","score","strand","fold_change","pvalue","qvalue","summit")
   ref = read.csv("~/Documents/Courses/Computational Analyses of sequence/Chipseq/chm13v2_geneList.txt",header = T,sep = "\t")
   mrna = read.csv("~/Documents/Courses/Computational Analyses of sequence/mRNA/K562_hs60_against_hs0_mRNAseq_DESeq2output.txt",
                header = T,sep = "\t") %>%
   rownames_to_column(., var = "refGeneName") %>%
   right_join(ref, ., by = "refGeneName")

   pl = subset(mrna, strand=="+"
   pl$TSS = pl$txStart
   pl$CPS = pl$txEnd

   mn = subset(mrna, strand=="-")
   mn$TSS = mn$txEnd
   mn$CPS = mn$txStart

   mrna = rbind(pl,mn)
   mrna <- mrna[with(mrna, order(chr, TSS)),] ## reorder the chromosome position
   mrna <- mrna %>% mutate(TSS.1 = TSS) %>% select(1,14,16,2:13,15)
   #write.table(mrna, "chm_mRNA_TSS.txt",col.names = T,row.names = F,sep = "\t",quote = F)

   ## manage the JunD narrowpeak file
   jund.igg <- read.csv("3 jund.IgGfiltered.narrowpeak",header = F,sep = "\t")
   data = jund.igg
   colnames(data) <- names
   data$summitCoordinate = data$peakStart + data$summit
   data <-  data %>%
   mutate(summitCoordinate.1 = summitCoordinate) %>%
   select(chr, summitCoordinate,summitCoordinate.1, peakStart, peakEnd, name, score,strand,
         fold_change, pvalue, qvalue, summit) %>%
   arrange(chr,summitCoordinate)
   write.table(data,"jund.igg.summitcoord.txt",quote = F,sep = "\t",row.names = F,col.names = F) 
   ```
   8. Identifying the nearest genes for JunD binding sites (Linux) \
      `bedtools closest -D b -a jund.igg.summitCoord.txt -b chm_mRNA_TSS.txt > jund_distance_to_TSS.bed`
   10.  
