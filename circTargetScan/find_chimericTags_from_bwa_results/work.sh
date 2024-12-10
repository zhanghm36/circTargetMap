mkdir r1r2_separate
perl step1.print_read_blocks.pl read1.toGenome.bwa.sam 1 > ./r1r2_separate/out1.read1.toGenome.bwa.blocks
perl step2.find_circ_Link.v2.pl ./r1r2_separate/out1.read1.toGenome.bwa.blocks circAtlas.spliceSite_circRNA.good.bedpair > ./r1r2_separate/out2.read1.circRNA_linked.blocks
perl step3.remove_potential_intraPairing.v2.pl ./r1r2_separate/out2.read1.circRNA_linked.blocks > ./r1r2_separate/out3.read1.circRNA_linked.targets 2> ./r1r2_separate/out4.read1.circRNA_linked.bed
perl step4.reMap_target_loci.pl ./r1r2_separate/out4.read1.circRNA_linked.bed *read1_torRNA_Unmapped_really.fq.gz > ./r1r2_separate/out5.read1.target_block.readSeq.fq
bowtie2 --sensitive --end-to-end -p 12 -N 1 -L 17 -x genome.fa -U ./r1r2_separate/out5.read1.target_block.readSeq.fq -S ./r1r2_separate/out6.read1.target_block.readSeq.reMap.bowtie2.sam 2> ./r1r2_separate/out6.read1.target_block.readSeq.reMap.bowtie2.log
hisat2 -p 12 -x genome.fa -U ./r1r2_separate/out5.read1.target_block.readSeq.fq -S ./r1r2_separate/out6.read1.target_block.readSeq.reMap.hisat2.sam 2> ./r1r2_separate/out6.read1.target_block.readSeq.reMap.hisat2.log
perl fq2fa.pl ./r1r2_separate/out5.read1.target_block.readSeq.fq > ./r1r2_separate/out5.read1.target_block.readSeq.fa
blastn -query ./r1r2_separate/out5.read1.target_block.readSeq.fa -task megablast -db human_rRNA.fa -out ./r1r2_separate/out6.read1.target_block.readSeq.blast_rRNA.out -evalue 1 -word_size 5 -outfmt 7 -num_threads 16
perl step5.filter_low_MAPQ.pl ./r1r2_separate/out3.read1.circRNA_linked.targets ./r1r2_separate/out4.read1.circRNA_linked.bed ./r1r2_separate/out6.read1.target_block.readSeq.reMap.bowtie2.sam ./r1r2_separate/out6.read1.target_block.readSeq.reMap.hisat2.sam ./r1r2_separate/out6.read1.target_block.readSeq.blast_rRNA.out read1.circRNA_linked 2> log.filtered.read1.blocks

perl step1.print_read_blocks.pl read2.toGenome.bwa.sam 2 > ./r1r2_separate/out1.read2.toGenome.bwa.blocks
perl step2.find_circ_Link.v2.pl ./r1r2_separate/out1.read2.toGenome.bwa.blocks  circAtlas.spliceSite_circRNA.good.bedpair > ./r1r2_separate/out2.read2.circRNA_linked.blocks
perl step3.remove_potential_intraPairing.v2.pl ./r1r2_separate/out2.read2.circRNA_linked.blocks > ./r1r2_separate/out3.read2.circRNA_linked.targets 2> ./r1r2_separate/out4.read2.circRNA_linked.bed
perl step4.reMap_target_loci.pl ./r1r2_separate/out4.read2.circRNA_linked.bed *read2_torRNA_Unmapped_really.fq.gz > ./r1r2_separate/out5.read2.target_block.readSeq.fq
bowtie2 --sensitive --end-to-end -p 12 -N 1 -L 17 -x genome.fa -U ./r1r2_separate/out5.read2.target_block.readSeq.fq -S ./r1r2_separate/out6.read2.target_block.readSeq.reMap.bowtie2.sam 2> ./r1r2_separate/out6.read2.target_block.readSeq.reMap.bowtie2.log
hisat2 -p 12 -x genome.fa -U ./r1r2_separate/out5.read2.target_block.readSeq.fq -S ./r1r2_separate/out6.read2.target_block.readSeq.reMap.hisat2.sam 2> ./r1r2_separate/out6.read2.target_block.readSeq.reMap.hisat2.log
perl fq2fa.pl ./r1r2_separate/out5.read2.target_block.readSeq.fq > ./r1r2_separate/out5.read2.target_block.readSeq.fa
blastn -query ./r1r2_separate/out5.read2.target_block.readSeq.fa -task megablast -db human_rRNA.fa -out ./r1r2_separate/out6.read2.target_block.readSeq.blast_rRNA.out -evalue 1 -word_size 5 -outfmt 7 -num_threads 16
perl step5.filter_low_MAPQ.pl ./r1r2_separate/out3.read2.circRNA_linked.targets ./r1r2_separate/out4.read2.circRNA_linked.bed ./r1r2_separate/out6.read2.target_block.readSeq.reMap.bowtie2.sam ./r1r2_separate/out6.read2.target_block.readSeq.reMap.hisat2.sam ./r1r2_separate/out6.read2.target_block.readSeq.blast_rRNA.out read2.circRNA_linked 2> log.filtered.read2.blocks
mv out7* ./r1r2_separate

