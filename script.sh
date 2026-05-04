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

# samtools view -b -F 1024 outputs/SRR11647649_dedp.bam > outputs/SRR11647649_deduplicated.bam
# samtools view -b -F 1024 outputs/SRR11647651_dedp.bam > outputs/SRR11647651_deduplicated.bam
# samtools view -b -F 1024 outputs/SRR11647659_dedp.bam > outputs/SRR11647659_deduplicated.bam
# samtools view -b -F 1024 outputs/SRR11647661_dedp.bam > outputs/SRR11647661_deduplicated.bam

conda activate methyldackel_env
MethylDackel mbias references/GRCh38.primary_assembly.genome.fa outputs/SRR11647649_dedp.bam mbias_plot_SRR11647649
MethylDackel mbias references/GRCh38.primary_assembly.genome.fa outputs/SRR11647651_dedp.bam mbias_plot_SRR11647651
MethylDackel mbias references/GRCh38.primary_assembly.genome.fa outputs/SRR11647659_dedp.bam mbias_plot_SRR11647659
MethylDackel mbias references/GRCh38.primary_assembly.genome.fa outputs/SRR11647661_dedp.bam mbias_plot_SRR11647661

samtools view -@ 8 -b outputs/SRR11647649_dedp.bam chr20 > outputs/SRR11647649_chr20.bam
samtools view -@ 8 -b outputs/SRR11647651_dedp.bam chr20 > outputs/SRR11647651_chr20.bam
samtools view -@ 8 -b outputs/SRR11647659_dedp.bam chr20 > outputs/SRR11647659_chr20.bam
samtools view -@ 8 -b outputs/SRR11647661_dedp.bam chr20 > outputs/SRR11647661_chr20.bam

samtools index -@ 8 outputs/SRR11647649_chr20.bam
samtools index -@ 8 outputs/SRR11647651_chr20.bam
samtools index -@ 8 outputs/SRR11647659_chr20.bam
samtools index -@ 8 outputs/SRR11647661_chr20.bam

MethylDackel extract --OT 2,0,2,111 --OB 2,124,14,124 references/GRCh38.primary_assembly.genome.fa outputs/SRR11647649_chr20.bam
echo "done 1 chr only"
MethylDackel extract --OT 2,0,2,112 --OB 2,124,13,124 references/GRCh38.primary_assembly.genome.fa outputs/SRR11647651_chr20.bam 
echo "done 2 chr only"
MethylDackel extract --OT 3,0,2,120 --OB 2,123,5,124 references/GRCh38.primary_assembly.genome.fa outputs/SRR11647659_chr20.bam 
echo "done 3 chr only"
MethylDackel extract --OT 6,0,2,112 --OB 2,123,13,124 references/GRCh38.primary_assembly.genome.fa outputs/SRR11647661_chr20.bam
echo "done 4 chr only"

MethylDackel extract --OT 2,0,2,111 --OB 2,124,14,124 references/GRCh38.primary_assembly.genome.fa outputs/SRR11647649_dedp.bam
echo "done 1"
MethylDackel extract --OT 2,0,2,112 --OB 2,124,13,124 references/GRCh38.primary_assembly.genome.fa outputs/SRR11647651_dedp.bam 
echo "done 2"
MethylDackel extract --OT 3,0,2,120 --OB 2,123,5,124 references/GRCh38.primary_assembly.genome.fa outputs/SRR11647659_dedp.bam 
echo "done 3"
MethylDackel extract --OT 6,0,2,112 --OB 2,123,13,124 references/GRCh38.primary_assembly.genome.fa outputs/SRR11647661_dedp.bam
echo "done 4"

mv *.svg outputs/