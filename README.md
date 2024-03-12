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
5. Visualising reads and peaks in IGV (Linux)
   ```bash
   for bdg in *_treat_pileup.bdg; do
    name=$(echo $bdg | awk -F"/" '{print $NF}' | awk -F"_treat_pileup" '{print $1}') echo $name
    bedGraphToBigWig ${name}_treat_pileup.bdg CHM13v2chromSizes.txt ${name}.bigWig;
   done
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
9. Identifying the JunD-bound promoters (R and Linux)
    ```R
    ## Generating the promoter region in reference genome 
    mrna$qStart = mrna$TSS-2500
    mrna$qEnd = mrna$TSS+2500
    write.table(mrna[,c("chr", "qStart", "qEnd", "geneName")], file="chm13_promoter_as_TSSpm2500.txt", col.names=F, row.names=F, quote=F, sep="\t")
   ```
   ```bash
   bedtools intersect -wa -a jund_distance_to_TSS.bed -b chm13_promoter_as_TSSpm2500.txt > jund_bound_promoters.bed
   awk '{print $NF}' jund_bound_promoters.bed > jund_bound_promoters.txt
   ```
   ```R
   ## extract promoter genes to submit to DAVID website
   go = read.csv("jund_bound_promoters.bed",header = F,sep = "\t")
   genes  = go$V19 %>% unique()
   write.table(genes, "promoter.genes.txt",row.names = F,col.names = F,quote = F)

   ## GO terms visualization  
   goterm = read.csv("GOterm.txt",header = T,sep = "\t") %>%
      separate(Term,into = c("col1","Term"),sep = "~") %>% select(-col1) %>%
      separate(Category, into = c("c1","Category","c2"),sep = "_") %>% select(-c2,-c1)

   ## CC, BP, and MF were conducted by the following codes separately
   cc = goterm %>% filter(Category == "BP") %>% slice(1:10) %>% mutate(log = -log10(FDR))
   KEGG = cc
   KEGG$Term[7] = "proteasome-mediated ubiquitin-dependent\n protein catabolic process"
   KEGG$Term <- fct_reorder(KEGG$Term,KEGG$PValue,.desc = T)
   p = ggplot(data = KEGG, aes(x =log,y = Term,color = FDR))+
     theme_bw(base_size = 10,base_line_size = 1,base_rect_size = 2)+
     theme(aspect.ratio = 1)+
     geom_point(stat = 'Identity',size = 5,alpha = 0.8)+
     scale_color_gradient(low = 'red',high = 'blue')+
     ylab(label = 'GO_Biological Process')+
     xlab(label = expression(paste("-Log"[10], "(padj)")))+
     labs(color = "padj")+
     ggtitle("GO_Biological Process")+
     theme(axis.title.x = element_text(size = 15,colour = 'black',face = 'bold',family = 'Times'),
        axis.text.x = element_text(size = 13,colour = 'black',face = 'bold',family = 'Times'),
        axis.text.y = element_text(size = 13,colour = 'black',face = 'bold',family = 'Times'),
        axis.title.y = element_text(size = 13,colour = 'black',face = 'bold',family = 'Times'),
        plot.title = element_text(size = 20,colour = 'black',face = 'bold',family = 'Times'))+
     theme(legend.key.size = unit(0.6,'cm'),
        #legend.key.height = unit(0.3,'cm'),
        legend.text = element_text(family = 'Times',size = 8,face = 'bold'),
        legend.title = element_text(family = 'Times',size = 10,face = 'bold',vjust = 0.95),
        legend.position = 'right')+
     scale_x_continuous(breaks = seq(round(min(KEGG$log),digits = 0)-1,round(max(KEGG$log),digits = 0),10))  
   
   ```
10. Functional genomic regions analysis (Linux and R)
    ```bash
    bedtools intersect -loj -wa -a jund_distance_to_TSS.bed -b chm13v2_functionalGenomicRegions_withRegionName.bed > jund_at_functional_genomic_regions.bed
    ```
    ```R
    # Functional genomic region plots
      
      jund_fGRs = read.table("6_jund_at_functional_genomic_regions.bed", sep="\t", header = F)
      
      jund_fGRs$regionOrder=7
      jund_fGRs$regionOrder<-ifelse(jund_fGRs$V39=="promoter", 1, jund_fGRs$regionOrder)
      jund_fGRs$regionOrder<-ifelse(jund_fGRs$V39=="enhancer", 2, jund_fGRs$regionOrder)
      jund_fGRs$regionOrder<-ifelse(jund_fGRs$V39=="divergent", 3, jund_fGRs$regionOrder)
      jund_fGRs$regionOrder<-ifelse(jund_fGRs$V39=="geneBody", 4, jund_fGRs$regionOrder)
      jund_fGRs$regionOrder<-ifelse(jund_fGRs$V39=="CPS", 5, jund_fGRs$regionOrder)
      jund_fGRs$regionOrder<-ifelse(jund_fGRs$V39=="TW", 6, jund_fGRs$regionOrder)
      jund_fGRs$V39 <-ifelse(jund_fGRs$regionOrder==7, "untranscribed", jund_fGRs$V39)
      
      jund_fGRs = jund_fGRs[order(jund_fGRs$regionOrder),]
      jund_fGRs = jund_fGRs[!duplicated(jund_fGRs$V6),]
      
      jund_fGRs$rgbCol ="grey"
      jund_fGRs$rgbCol <- ifelse(jund_fGRs$V39=="promoter","#f38400", jund_fGRs$rgbCol)
      jund_fGRs$rgbCol <- ifelse(jund_fGRs$V39=="enhancer","#73d47a", jund_fGRs$rgbCol)
      jund_fGRs$rgbCol <- ifelse(jund_fGRs$V39=="divergent","#b23bd4", jund_fGRs$rgbCol)
      jund_fGRs$rgbCol <- ifelse(jund_fGRs$V39=="geneBody","yellow", jund_fGRs$rgbCol)
      jund_fGRs$rgbCol <- ifelse(jund_fGRs$V39=="CPS","#67c8f9", jund_fGRs$rgbCol)
      jund_fGRs$rgbCol <- ifelse(jund_fGRs$V39=="TW","#ff3662", jund_fGRs$rgbCol)
      
      jund_fGRs$region = factor(jund_fGRs$V39, c("promoter", "enhancer", "divergent", "geneBody", "CPS", "TW", "untranscribed"))
      
      a = as.data.frame(table(jund_fGRs$V39)/dim(jund_fGRs)[1]*100)
      names(a) =c("region", "pct")
      
      fGRs_score_plot = boxplot(jund_fGRs$V7~jund_fGRs$V39, col=unique(jund_fGRs$rgbCol), log = "y", xlab = "JUND region", ylab ="score")
      
      fGRs_percentage_plot = ggplot(a, aes(x = "", y = pct, fill = region)) +
        geom_bar(stat = "identity", position = position_dodge(width = 1)) +
        labs(x = "Genomic region", y = "Percentage") + 
        scale_fill_manual(values = unique(jund_fGRs$rgbCol)) +
        theme_minimal() +
        ylim(0, 50)
      
      # Selecting top 1000 scores
      
      top_scores_fGRs <- jund_fGRs %>%
        top_n(1000, V7) %>% arrange(desc(V7))
      
      a_top_1000 = as.data.frame(table(top_scores_fGRs$V39)/dim(top_scores_fGRs)[1]*100)
      names(a_top_1000) =c("region", "pct")
      
      fGRs_score_plot_top = boxplot(top_scores_fGRs$V7~top_scores_fGRs$V39, col=unique(jund_fGRs$rgbCol), log = "y", xlab = "JUND region", ylab ="score")
      
      fGRs_percentage_plot_top = ggplot(a_top_1000, aes(x = "", y = pct, fill = region)) +
        geom_bar(stat = "identity", position = position_dodge(width = 1)) +
        labs(x = "Genomic region", y = "Percentage") + 
        scale_fill_manual(values = unique(jund_fGRs$rgbCol)) +
        theme_minimal() +
        ylim(0, 50)
      
      library(gridExtra)
      grid.arrange(fGRs_percentage_plot, fGRs_percentage_plot_top, ncol = 2)
      
      # Top 5 scored peaks per functional genomic region
      
      top_scores <- jund_fGRs %>%
        group_by(V39) %>%
        slice(which.max(V7)) %>%
        select(V39, top_score = V7, V19)
      
      top_5_scores <- jund_fGRs %>%
        group_by(V39) %>%
        top_n(5, V7) %>%
        arrange(V39, desc(V7)) %>%
        select(V39, V7, V19, V6, V13)
      print(top_5_scores, n=35)
      
      # MEME analysis
      
      summitPM50 = read.csv("jund.igg.summitcoord.txt",header = F, sep = "\t")
      names(summitPM50) = c("chr","summitCoordinate","summitCoordinate.1","peakStart","peakEnd","name",
                            "score","strand","fold_change","pvalue","qvalue","summit")
      
      top_scores <- summitPM50 %>%
        top_n(1000, score) %>% arrange(desc(score))
      
      meme = top_scores %>% select(chr,summitCoordinate) %>%
        mutate(qStart = summitCoordinate -50) %>%
        mutate(qEnd = summitCoordinate+ 50) %>% select(-summitCoordinate)
        
      write.table(meme, "jund_summitPM50_top_1000.txt", col.names=F, row.names=F, quote=F, sep="\t")
    ```
11. De novo motif analyses with MEME-CHIP (R and Linux) 
    ```R
    summitPM50 = read.csv("jund.igg.summitcoord.txt",header = F, sep = "\t")
    names(summitPM50) = c("chr","summitCoordinate","summitCoordinate.1","peakStart","peakEnd","name",
                      "score","strand","fold_change","pvalue","qvalue","summit")
    meme = summitPM50 %>% select(chr,summitCoordinate) %>%
       mutate(qStart = summitCoordinate -50) %>%
       mutate(qEnd = summitCoordinate+ 50) %>% select(-summitCoordinate)
    write.table(meme, "jund_summitPM50.txt", col.names=F, row.names=F, quote=F, sep="\t")
    ```
    ```bash
    bedtools getfasta -fi chm13v2.0.fa -bed jund_summitPM50.txt > jund_summitPM50.fasta
    ```
