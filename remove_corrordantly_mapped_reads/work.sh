#Map reads to the genome using HISAT2
hisat2 -x $genome_ref_fasta -U /read1.fq -S /read1.toGenome.hisat2.sam --un  --no-softclip --phred33 -p 8 2> /read1.toGenome.hisat2.log
hisat2 -x $genome_ref_fasta -U /read2.fq -S /read2.toGenome.hisat2.sam --un  --no-softclip --phred33 -p 8 2> /read2.toGenome.hisat2.log

#Discard concordantly mapped Read 1 and Read 2 using HISAT2
perl remove_both_mapping_reads.pl  read1.fq  read2.fq  read1.toGenome.hisat2.sam  read2.toGenome.hisat2.sam  read1.toGenome.unmapped.any.byHisat2.fq  read2.toGenome.unmapped.any.byHisat2.fq

#Remap unmapped reads to the genome using BWA, which supports for chimeric alignment
bwa mem -t 16 -k 12 -T 15 -o /read1.toGenome.bwa.sam $genome_ref_fasta /read1.toGenome.unmapped.any.byHisat2.fq
bwa mem -t 16 -k 12 -T 15 -o /read2.toGenome.bwa.sam $genome_ref_fasta /read2.toGenome.unmapped.any.byHisat2.fq
