#!/bin/bash

if [ $# -lt 4 ]
then
    echo "sh $0 <species.fa/softmasked.fa> <species_dir> <softmasked.fa> <species_id>"
    exit 0
fi

fa=$1
dir=$2
Rfa=$3
id=$4

## protein sequences were mapped to DNA sequences. Blast and Genewise were used. 'genome.dna.fa' represents the genome assembly for a specific bird. Please note that this command line submitted parallel jobs to BGI Linux cluster (SGE).
###GALGA, Galgal4.85
perl ./1.bin/gene_finding/protein-map-genome/bin/protein_map_genome_b10k_phase2.pl --verbose --cpu 50 --resource vf=1G --run qsub --blast_eval 1e-5 --filter_rate 0.5 --extend_len 2000 --step 123468 --lines 500 --reqsub --outdir $dir/GALGA/ ./0.homolog_data/GALGA.pep $fa

###TAEGU, taeGut3.2.4.85
perl ./1.bin/gene_finding/protein-map-genome/bin/protein_map_genome_b10k_phase2.pl --verbose --cpu 50 --resource vf=1G --run qsub --blast_eval 1e-5 --filter_rate 0.5 --extend_len 2000 --step 123468 --lines 500 --reqsub --outdir $dir/TAEGU/ ./0.homolog_data/TAEGU.pep $fa

###HUMAN, GRCh38
perl ./1.bin/gene_finding/protein-map-genome/bin/protein_map_genome_b10k_phase2.pl --verbose --cpu 50 --resource vf=1G --run qsub --blast_eval 1e-5 --filter_rate 0.5 --extend_len 2000 --step 123468 --lines 500 --reqsub --outdir $dir/HUMAN/ ./0.homolog_data/HUMAN.pep $fa

###TRAN2, multiple homology transcriptome
perl ./1.bin/gene_finding/protein-map-genome/bin/protein_map_genome_b10k_phase2.pl --verbose --cpu 50 --resource vf=1G --run qsub --blast_eval 1e-5 --filter_rate 0.5 --extend_len 2000 --step 123468 --lines 500 --reqsub --outdir $dir/TRAN2/ ./0.homolog_data/TRAN2.pep $fa

###Newfilter & combine
sh ./1.bin/shell_bin/combine.pl $dir $Rfa $id
