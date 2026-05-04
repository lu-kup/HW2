# picard MarkDuplicates I=outputs/SRR11647661.bam O=outputs/SRR11647661_dedp.bam M=outputs/SRR11647661_metrics.txt REMOVE_DUPLICATES=false
# samtools view -b -F 1024 outputs/SRR11647649_dedp.bam > outputs/SRR11647649_deduplicated.bam
# samtools view -b -F 1024 outputs/SRR11647651_dedp.bam > outputs/SRR11647651_deduplicated.bam
# samtools view -b -F 1024 outputs/SRR11647659_dedp.bam > outputs/SRR11647659_deduplicated.bam
# samtools view -b -F 1024 outputs/SRR11647661_dedp.bam > outputs/SRR11647661_deduplicated.bam

# conda activate methyldackel_env
# MethylDackel mbias references/GRCh38.primary_assembly.genome.fa outputs/SRR11647649_dedp.bam mbias_plot_SRR11647649
# MethylDackel mbias references/GRCh38.primary_assembly.genome.fa outputs/SRR11647651_dedp.bam mbias_plot_SRR11647651
# MethylDackel mbias references/GRCh38.primary_assembly.genome.fa outputs/SRR11647659_dedp.bam mbias_plot_SRR11647659
# MethylDackel mbias references/GRCh38.primary_assembly.genome.fa outputs/SRR11647661_dedp.bam mbias_plot_SRR11647661

# samtools index -@ 8 outputs/SRR11647649_chr20.bam
# samtools index -@ 8 outputs/SRR11647651_chr20.bam
# samtools index -@ 8 outputs/SRR11647659_chr20.bam
# samtools index -@ 8 outputs/SRR11647661_chr20.bam

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