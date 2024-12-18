mkdir r1r2_separate
#######Scan the chimeric reads supporting circRNA-target RNA interactions from read1########
perl step1.print_read_blocks.pl read1.toGenome.bwa.sam 1 > ./r1r2_separate/out1.read1.toGenome.bwa.blocks
perl step2.find_circ_Link.pl ./r1r2_separate/out1.read1.toGenome.bwa.blocks circAtlas.circRNA_SpliceSite_10nt.bedpair > ./r1r2_separate/out2.read1.circRNA_linked.blocks
perl step3.remove_potential_intraPairing.pl ./r1r2_separate/out2.read1.circRNA_linked.blocks > ./r1r2_separate/out3.read1.circRNA_linked.targets 2> ./r1r2_separate/out4.read1.circRNA_linked.bed
perl step4.reMap_target_loci.pl ./r1r2_separate/out4.read1.circRNA_linked.bed  read1.fq > ./r1r2_separate/out5.read1.target_block.readSeq.fq
#Remap the candidate target reads to the genome by Bowtie2 and HISAT2
bowtie2 --sensitive --end-to-end -p 12 -N 1 -L 17 -x genome.fa -U ./r1r2_separate/out5.read1.target_block.readSeq.fq -S ./r1r2_separate/out6.read1.target_block.readSeq.reMap.bowtie2.sam 2> ./r1r2_separate/out6.read1.target_block.readSeq.reMap.bowtie2.log
hisat2 -p 12 -x genome.fa -U ./r1r2_separate/out5.read1.target_block.readSeq.fq -S ./r1r2_separate/out6.read1.target_block.readSeq.reMap.hisat2.sam 2> ./r1r2_separate/out6.read1.target_block.readSeq.reMap.hisat2.log
#Scan for the reads that mapped to rRNA
perl fq2fa.pl ./r1r2_separate/out5.read1.target_block.readSeq.fq > ./r1r2_separate/out5.read1.target_block.readSeq.fa
blastn -query ./r1r2_separate/out5.read1.target_block.readSeq.fa -task megablast -db human_rRNA.fa -out ./r1r2_separate/out6.read1.target_block.readSeq.blast_rRNA.out -evalue 1 -word_size 5 -outfmt 7 -num_threads 16
#Remove the reads thay multiply mapped to the genome or mapped to rRNA
perl step5.filter_low_MAPQ.pl ./r1r2_separate/out3.read1.circRNA_linked.targets ./r1r2_separate/out4.read1.circRNA_linked.bed ./r1r2_separate/out6.read1.target_block.readSeq.reMap.bowtie2.sam ./r1r2_separate/out6.read1.target_block.readSeq.reMap.hisat2.sam ./r1r2_separate/out6.read1.target_block.readSeq.blast_rRNA.out read1.circRNA_linked 2> log.filtered.read1.blocks

#######Scan the chimeric reads supporting circRNA-target RNA interactions from read1########
perl step1.print_read_blocks.pl read2.toGenome.bwa.sam 2 > ./r1r2_separate/out1.read2.toGenome.bwa.blocks
perl step2.find_circ_Link.pl ./r1r2_separate/out1.read2.toGenome.bwa.blocks  circAtlas.circRNA_SpliceSite_10nt.bedpair >./r1r2_separate/out2.read2.circRNA_linked.blocks
perl step3.remove_potential_intraPairing.pl ./r1r2_separate/out2.read2.circRNA_linked.blocks > ./r1r2_separate/out3.read2.circRNA_linked.targets 2> ./r1r2_separate/out4.read2.circRNA_linked.bed
perl step4.reMap_target_loci.pl ./r1r2_separate/out4.read2.circRNA_linked.bed  read2.fq > ./r1r2_separate/out5.read2.target_block.readSeq.fq
#Remap the candidate target reads to the genome by Bowtie2 and HISAT2
bowtie2 --sensitive --end-to-end -p 12 -N 1 -L 17 -x genome.fa -U ./r1r2_separate/out5.read2.target_block.readSeq.fq -S ./r1r2_separate/out6.read2.target_block.readSeq.reMap.bowtie2.sam 2> ./r1r2_separate/out6.read2.target_block.readSeq.reMap.bowtie2.log
hisat2 -p 12 -x genome.fa -U ./r1r2_separate/out5.read2.target_block.readSeq.fq -S ./r1r2_separate/out6.read2.target_block.readSeq.reMap.hisat2.sam 2> ./r1r2_separate/out6.read2.target_block.readSeq.reMap.hisat2.log
#Scan for the reads that mapped to rRNA
perl fq2fa.pl ./r1r2_separate/out5.read2.target_block.readSeq.fq > ./r1r2_separate/out5.read2.target_block.readSeq.fa
blastn -query ./r1r2_separate/out5.read2.target_block.readSeq.fa -task megablast -db human_rRNA.fa -out ./r1r2_separate/out6.read2.target_block.readSeq.blast_rRNA.out -evalue 1 -word_size 5 -outfmt 7 -num_threads 16
#Remove the reads thay multiply mapped to the genome or mapped to rRNA
perl step5.filter_low_MAPQ.pl ./r1r2_separate/out3.read2.circRNA_linked.targets ./r1r2_separate/out4.read2.circRNA_linked.bed ./r1r2_separate/out6.read2.target_block.readSeq.reMap.bowtie2.sam ./r1r2_separate/out6.read2.target_block.readSeq.reMap.hisat2.sam ./r1r2_separate/out6.read2.target_block.readSeq.blast_rRNA.out read2.circRNA_linked 2> log.filtered.read2.blocks
mv out7* ./r1r2_separate
