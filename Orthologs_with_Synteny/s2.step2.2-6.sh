export PATH="/usr/hal/bin/:/usr/mafft-7.402/bin/:$PATH"
export SCRIPTS_DIR="/usr/ortholog_B10K/"
ref=$1
tag=$2
refId=$3
tagId=$4
export ALL_INPUT_FILES="/usr/ortholog_B10K/test/input/"
export ALL_OUTPUT_FILES="/usr/ortholog_B10K/test/"

echo ==========start at : `date` ==========

cd ${ALL_OUTPUT_FILES}/${ref}.${tag}

echo step2.2: find the intersection of putative orthologous regions from Cactus
awk '{print "cat "$NF}' ${ALL_OUTPUT_FILES}/${ref}.${tag}/split_halLiftover_${ref}.sh | sh | sort -k1,1 -k2,2n | bedtools intersect -a ${ALL_OUTPUT_FILES}/${ref}.${tag}/${tag}.cds.sort.bed -b - -wo > ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.${tag}.intersect
perl ${SCRIPTS_DIR}/92.get_overlap.pl ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.${tag}.intersect ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.cds.sort.bed ${ALL_OUTPUT_FILES}/${ref}.${tag}/${tag}.cds.sort.bed > ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.${tag}.intersect.ortholog.list

echo step4: gene order
perl ${SCRIPTS_DIR}/94.geneOrder_add_ortholog.pl ${ALL_INPUT_FILES}/${ref}.gff ${ALL_INPUT_FILES}/${tag}.gff ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.${tag}.intersect.ortholog.list > ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.${tag}.intersect.ortholog.list.order

perl ${SCRIPTS_DIR}/95.syntenic_ort.syn.pl ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.${tag}.intersect.ortholog.list.order ${ALL_INPUT_FILES}/${ref}.gff ${ALL_INPUT_FILES}/${tag}.gff > ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.${tag}.intersect.ortholog.list.order.ort

perl ${SCRIPTS_DIR}/96.ort2cor.syn.pl 3 ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.${tag}.intersect.ortholog.list.order.ort > ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.${tag}.intersect.ortholog.list.order.ort.cor

echo step5: RBH based on gene synteny
perl ${SCRIPTS_DIR}/97.ort_rbh.syn.pl ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.${tag}.intersect.ortholog.list.order.ort.cor ${refId} ${tagId} > ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.${tag}.intersect.ortholog.list.order.ort.cor.rbh

awk '{print $2"\t"$1"\t"$5"\t"$6}' ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.${tag}.intersect.ortholog.list.order.ort.cor.rbh | sort -k1,1 > ${ALL_OUTPUT_FILES}/${ref}.${tag}/${refId}_${tagId}.orthlog

echo step6: RBH based on gene synteny and retain all the Tandem
perl ${SCRIPTS_DIR}/97.ort_rbh.syn.all.pl ${ALL_OUTPUT_FILES}/${ref}.${tag}/${ref}.${tag}.intersect.ortholog.list.order.ort.cor ${refId} ${tagId} > ${ALL_OUTPUT_FILES}/${ref}.${tag}/${refId}_${tagId}.orthlog.all

cd -
echo ==========end at : `date` ==========
echo run 98.merge_table.pl list_of_all_.ortholog_files to get merged table, or 98.merge_table.v2.pl list_of_.ortholog.all_files to get merged table containing Tandem

