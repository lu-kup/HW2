# mkdir qc
# fastqc raw_data/*.fastq.gz -t 8 -o qc

trim_galore -o reads_trimmed -j 8 --quality 20 --length 20 --paired raw_data/SRR11647649_1.fastq.gz raw_data/SRR11647649_2.fastq.gz
rm raw_data/SRR11647649_1.fastq.gz
rm raw_data/SRR11647649_2.fastq.gz

trim_galore -o reads_trimmed -j 8 --quality 20 --length 20 --paired raw_data/SRR11647651_1.fastq.gz raw_data/SRR11647651_2.fastq.gz
rm raw_data/SRR11647651_1.fastq.gz
rm raw_data/SRR11647651_2.fastq.gz

trim_galore -o reads_trimmed -j 8 --quality 20 --length 20 --paired raw_data/SRR11647659_1.fastq.gz raw_data/SRR11647659_2.fastq.gz
rm raw_data/SRR11647659_1.fastq.gz
rm raw_data/SRR11647659_2.fastq.gz

trim_galore -o reads_trimmed -j 8 --quality 20 --length 20 --paired raw_data/SRR11647661_1.fastq.gz raw_data/SRR11647661_2.fastq.gz
rm raw_data/SRR11647661_1.fastq.gz
rm raw_data/SRR11647661_2.fastq.gz

mkdir qc_trimmed
fastqc reads_trimmed/*.fq.gz -t 8 -o qc_trimmed

mkdir outputs

# genomo indeksavimas
conda activate bwameth_env
bwameth.py index references/GRCh38.primary_assembly.genome.fa 

# mapping

# bwameth.py --threads 8 --reference references/GRCh38.primary_assembly.genome.fa reads_trimmed/SRR11647649_1_val_1.fq.gz reads_trimmed/SRR11647649_2_val_2.fq.gz | samtools sort -@ 8 -o outputs/SRR11647649.bam
# bwameth.py --threads 8 --reference references/GRCh38.primary_assembly.genome.fa reads_trimmed/SRR11647651_1_val_1.fq.gz reads_trimmed/SRR11647651_2_val_2.fq.gz | samtools sort -@ 8 -o outputs/SRR11647651.bam
# bwameth.py --threads 8 --reference references/GRCh38.primary_assembly.genome.fa reads_trimmed/SRR11647659_1_val_1.fq.gz reads_trimmed/SRR11647659_2_val_2.fq.gz | samtools sort -@ 8 -o outputs/SRR11647659.bam
# bwameth.py --threads 8 --reference references/GRCh38.primary_assembly.genome.fa reads_trimmed/SRR11647661_1_val_1.fq.gz reads_trimmed/SRR11647661_2_val_2.fq.gz | samtools sort -@ 8 -o outputs/SRR11647661.bam

bwameth.py --threads 8 --reference references/GRCh38.primary_assembly.genome.fa reads_trimmed/SRR11647649_1_val_1.fq.gz reads_trimmed/SRR11647649_2_val_2.fq.gz > outputs/SRR11647649.sam 2> outputs/SRR11647649.log
samtools sort -@ 8 -o outputs/SRR11647649.bam outputs/SRR11647649.sam 2> outputs/SRR11647649.bam.log

bwameth.py --threads 8 --reference references/GRCh38.primary_assembly.genome.fa reads_trimmed/SRR11647651_1_val_1.fq.gz reads_trimmed/SRR11647651_2_val_2.fq.gz > outputs/SRR11647651.sam 2> outputs/SRR11647651.log
samtools sort -@ 8 -o outputs/SRR11647651.bam outputs/SRR11647651.sam 2> outputs/SRR11647651.bam.log

bwameth.py --threads 8 --reference references/GRCh38.primary_assembly.genome.fa reads_trimmed/SRR11647659_1_val_1.fq.gz reads_trimmed/SRR11647659_2_val_2.fq.gz > outputs/SRR11647659.sam 2> outputs/SRR11647659.log
samtools sort -@ 8 -o outputs/SRR11647659.bam outputs/SRR11647659.sam 2> outputs/SRR11647659.bam.log

bwameth.py --threads 8 --reference references/GRCh38.primary_assembly.genome.fa reads_trimmed/SRR11647661_1_val_1.fq.gz reads_trimmed/SRR11647661_2_val_2.fq.gz > outputs/SRR11647661.sam 2> outputs/SRR11647661.log
samtools sort -@ 8 -o outputs/SRR11647661.bam outputs/SRR11647661.sam 2> outputs/SRR11647661.bam.log

conda activate picard

picard MarkDuplicates I=outputs/SRR11647649.bam O=outputs/SRR11647649_dedp.bam M=outputs/SRR11647649_metrics.txt REMOVE_DUPLICATES=false
picard MarkDuplicates I=outputs/SRR11647651.bam O=outputs/SRR11647651_dedp.bam M=outputs/SRR11647651_metrics.txt REMOVE_DUPLICATES=false
picard MarkDuplicates I=outputs/SRR11647659.bam O=outputs/SRR11647659_dedp.bam M=outputs/SRR11647659_metrics.txt REMOVE_DUPLICATES=false
picard MarkDuplicates I=outputs/SRR11647661.bam O=outputs/SRR11647661_dedp.bam M=outputs/SRR11647661_metrics.txt REMOVE_DUPLICATES=false

samtools view -b -F 1024 outputs/SRR11647649_dedp.bam > outputs/SRR11647649_deduplicated.bam
samtools view -b -F 1024 outputs/SRR11647651_dedp.bam > outputs/SRR11647651_deduplicated.bam
samtools view -b -F 1024 outputs/SRR11647659_dedp.bam > outputs/SRR11647659_deduplicated.bam
samtools view -b -F 1024 outputs/SRR11647661_dedp.bam > outputs/SRR11647661_deduplicated.bam

conda activate methyldackel_env
MethylDackel mbias references/GRCh38.primary_assembly.genome.fa outputs/SRR11647649_deduplicated.bam mbias_plot_SRR11647649
MethylDackel mbias references/GRCh38.primary_assembly.genome.fa outputs/SRR11647651_deduplicated.bam mbias_plot_SRR11647651
MethylDackel mbias references/GRCh38.primary_assembly.genome.fa outputs/SRR11647659_deduplicated.bam mbias_plot_SRR11647659
MethylDackel mbias references/GRCh38.primary_assembly.genome.fa outputs/SRR11647661_deduplicated.bam mbias_plot_SRR11647661

MethylDackel extract --OT 0,0,2,143 --OB 2,132,4,149 references/GRCh38.primary_assembly.genome.fa outputs/SRR11647649_deduplicated.bam sample_methylation_SRR11647649
MethylDackel extract --OT 0,0,2,143 --OB 2,132,4,149 references/GRCh38.primary_assembly.genome.fa outputs/SRR11647651_deduplicated.bam sample_methylation_SRR11647651
MethylDackel extract --OT 0,0,2,143 --OB 2,132,4,149 references/GRCh38.primary_assembly.genome.fa outputs/SRR11647659_deduplicated.bam sample_methylation_SRR11647659
MethylDackel extract --OT 0,0,2,143 --OB 2,132,4,149 references/GRCh38.primary_assembly.genome.fa outputs/SRR11647661_deduplicated.bam sample_methylation_SRR11647661
