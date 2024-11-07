#!/bin/bash
#

set -e

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

CWD=`dirname ${CWD}`

help() {
  echo "usage: `basename $0` -o <output_dir> <input_file> ..."
  exit $1
}

output_dir="."
while getopts 'x:o:h' OPT; do
  case $OPT in
    o) output_dir="$OPTARG";;
    h) help 0;;
    ?) help 1;;
  esac
done
shift $((OPTIND - 1))

export PYTHONPATH=${CWD}/profold2

output_params="?attr_idx=attr.idx_all&mapping_idx=mapping.idx_all"


# convert csv to fasta files (TCR_pMHC data)
python ${CWD}/main.py csv_to_fasta -o "${output_dir}${output_params}" --target_uri "${output_dir}${output_params}" --pid_prefix tcr_pmhc_test_ --start_idx=0 -v $*

# make chain.idx
cat ${output_dir}/mapping.idx_all | \
  cut -f2 | \
  awk -F _ '{printf("%s",$1);for (i=2;i<NF;++i) printf("_%s", $i); printf(" %s\n", $NF);}' | \
  sort -T . | \
  awk -f ${CWD}/scripts/collapse.awk  > ${output_dir}/chain.idx_all

# build the dataset (test data included) mapping.idx and chain.idx
cat ${CWD}/data/tcr_pmhc_db/mapping.idx ${output_dir}/mapping.idx_all > ${output_dir}/mapping.idx
cat ${output_dir}/mapping.idx | \
  cut -f2 | \
  awk -F _ '{printf("%s",$1);for (i=2;i<NF;++i) printf("_%s", $i); printf(" %s\n", $NF);}' | \
  sort -T . | \
  awk -f ${CWD}/scripts/collapse.awk  > ${output_dir}/chain.idx


# build fasta for each chain
for c in "A" "B" "P" "M"; do
  echo ${c}
  find ${output_dir}/fasta -name "*_${c}.fasta" -exec awk '{print $0}' {} \; > ${output_dir}/tcr_pmhc_${c}.fa
  n=$(cat ${output_dir}/tcr_pmhc_${c}.fa | wc -l)
  if [ ${n} -eq 0 ]; then
    echo ">fake1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" > ${output_dir}/tcr_pmhc_${c}.fa
  fi
done

# align A B M with jackhmmer
for c in "A" "B" "M"; do
  find ${CWD}/data/tcr_pmhc_db/fasta -name "*_${c}.fasta" > ${output_dir}/tcr_pmhc_db_${c}
  cat ${output_dir}/tcr_pmhc_db_${c} | ${CWD}/bin/mapred -m "uniref90_db=${output_dir}/tcr_pmhc_${c}.fa mgnify_db=${CWD}/data/tcr_pmhc_db_${c}.fa PYTHONPATH=${CWD}/profold2 sh ${CWD}/scripts/run_jackhmmer.sh -o ${output_dir}/a3m" -c 10
  cat ${output_dir}/tcr_pmhc_db_${c} | ${CWD}/bin/mapred -m "PIPELINE_UNIREF_MAX_HITS=1000000 PIPELINE_MGNIFY_MAX_HITS=1000000 PYTHONPATH=${CWD}/profold2 sh ${CWD}/scripts/run_pipeline.sh -o ${output_dir}/a3m" -c 10
done

# align P with equal length
python ${CWD}/main.py align_peptide -o ${output_dir}/a3m --db ../db/tcr/tcr_pmhc_db_P.fa ${output_dir}/tcr_pmhc_P.fa -v ${CWD}/data/tcr_pmhc_db/fasta/*_P.fasta
for c in "P"; do
  find ${CWD}/data/tcr_pmhc_db/fasta -name "*_${c}.fasta" > ${output_dir}/tcr_pmhc_db_${c}
  cat ${output_dir}/tcr_pmhc_db_${c} | ${CWD}/bin/mapred -m "PIPELINE_UNIREF_MAX_HITS=1000000 PIPELINE_MGNIFY_MAX_HITS=1000000 PYTHONPATH=${CWD}/profold2 sh ${CWD}/scripts/run_pipeline.sh -o ${output_dir}/a3m" -c 10
done

# filter a3m with threshold=t
t="0.85"
if [ -d ${output_dir}/var ]; then
  rm -rf ${output_dir}/var
fi
echo "filter a3m: t=${t}"
cp -r ${output_dir}/a3m ${output_dir}/var
python a3m_filter.py --trim_gap -t ${t} -o ${output_dir}/var ${CWD}/data/tcr_pmhc_db/fasta/*_M.fasta

