export PATH="/hwfssz5/ST_DIVERSITY/B10K/PUB/USER/fangqi/local/bio-tool/hal/bin/:/hwfssz5/ST_DIVERSITY/B10K/PUB/local/basic_bio-tool/mafft-7.402/bin/:$PATH"
export SCRIPTS_DIR="/hwfssz5/ST_DIVERSITY/B10K/PUB/USER/fangqi/share/Cactus_Alignments_Tools/ortholog_B10K_synteny/"
ref=$1
tag=$2
export ALL_INPUT_FILES="/hwfssz5/ST_DIVERSITY/B10K/PUB/USER/fangqi/share/Cactus_Alignments_Tools/ortholog_B10K/test/input/"
export ALL_OUTPUT_FILES="/hwfssz5/ST_DIVERSITY/B10K/PUB/USER/fangqi/share/Cactus_Alignments_Tools/ortholog_B10K/test/"
hal="/hwfssz5/ST_DIVERSITY/B10K/USER/fangqi/03.alignments/03.363birds/00.data/birds-final.hal"

echo ==========start at : `date` ==========
echo step1: files preparation

[ -d ${ALL_OUTPUT_FILES}/${ref}.${tag} ] && rm -rf ${ALL_OUTPUT_FILES}/${ref}.${tag}
mkdir ${ALL_OUTPUT_FILES}/${ref}.${tag}

cd ${ALL_OUTPUT_FILES}/${ref}.${tag}

##reference: gff pep cds; target gff pep cds
[ -e ${ALL_INPUT_FILES}/${ref}.gff -a -e ${ALL_INPUT_FILES}/${tag}.gff ] || exit
awk '$3=="mRNA"{print $1"\t"$4-1"\t"$5"\t"$NF}' ${ALL_INPUT_FILES}/${ref}.gff | sed 's/;.*//g' | sed 's/ID=//g' | sort -k1,1 -k2,2n > ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.mrna.sort.bed
awk '$3=="CDS"{print $1"\t"$4-1"\t"$5"\t"$NF}' ${ALL_INPUT_FILES}/${ref}.gff | sed 's/;.*//g' | sed 's/Parent=//g' | sort -k1,1 -k2,2n > ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.cds.sort.bed
awk '$3=="CDS"{print $1"\t"$4-1"\t"$5"\t"$NF}' ${ALL_INPUT_FILES}/${tag}.gff | sed 's/;.*//g' | sed 's/Parent=//g' | sort -k1,1 -k2,2n > ${ALL_OUTPUT_FILES}/${ref}.${tag}/${tag}.cds.sort.bed

[ -e ${ALL_INPUT_FILES}/${ref}.cds.fa -a -e ${ALL_INPUT_FILES}/${tag}.cds.fa ] || exit

echo step2.1: find the intersection of putative orthologous regions from Cactus
[ -d ${ALL_OUTPUT_FILES}/${ref}.${tag}/aligned_pos/ ] && rm -rf ${ALL_OUTPUT_FILES}/${ref}.${tag}/aligned_pos/
mkdir ${ALL_OUTPUT_FILES}/${ref}.${tag}/aligned_pos/
perl ${SCRIPTS_DIR}/91.extract_aligned_pos.parallel.pl -i ${ref}.mrna.sort.bed -r ${ref} -t ${tag} -hal ${hal} -s ${ALL_OUTPUT_FILES}/${ref}.${tag}/aligned_pos/

echo PLEASE: perl /hwfssz5/ST_DIVERSITY/B10K/PUB/local/bin/qsub-sge.pl --resource vf=2G --num_proc 1 --convert no --lines 100 --maxjob 200 --jobprefix ex ${ALL_OUTPUT_FILES}/${ref}.${tag}/split_halLiftover_${ref}.sh
echo about 0.5~1.0 hour for each bird

cd -
echo ==========end at : `date` ==========
