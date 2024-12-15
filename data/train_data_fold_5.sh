#!/bin/bash
#

set -e

CWD=`readlink -f $0`
CWD=`dirname ${CWD}`

CWD=`dirname ${CWD}`

export PYTHONPATH=${CWD}/profold2

data_dir=${CWD}/data/tcr_pmhc_train_data_fold_5
output_dir=${data_dir}/output2

output_params="?attr_idx=attr.idx_all&mapping_idx=mapping.idx_all"


# convert csv to fasta files (pMHC data)
python ${CWD}/main.py csv_to_fasta -o "${output_dir}${output_params}" --target_uri "${output_dir}${output_params}" --pid_prefix pmhc_train_ --default_y=1.0 -v ${CWD}/data/HLA_Data2_with_negs_and_label.csv


# convert csv to fasta files (TCR_pMHC data)
for ((i=0;i<5;++i)); do
  start_idx=$((${i} * 1000000))
  echo -e "$i\t${start_idx}"
  python ${CWD}/main.py csv_to_fasta -o "${output_dir}${output_params}" --target_uri "${output_dir}${output_params}" --pid_prefix tcr_pmhc_train_ --start_idx=${start_idx} -v ${data_dir}/train${i}
  python ${CWD}/main.py csv_to_fasta -o "${output_dir}${output_params}" --target_uri "${output_dir}${output_params}" --pid_prefix tcr_pmhc_test_ --start_idx=${start_idx} -v ${data_dir}/test${i}
done

# make chain.idx
cat ${output_dir}/mapping.idx_all | \
  cut -f2 | \
  awk -F _ '{printf("%s",$1);for (i=2;i<NF;++i) printf("_%s", $i); printf(" %s\n", $NF);}' | \
  sort -T . | \
  awk -f ${CWD}/scripts/collapse.awk  > ${output_dir}/chain.idx_all

# filter out ones that has only one chain
#   1. load dict a (negotives or in test dataset) from attr.idx_all
#   2. filter out those that:
#      i.  has no peptide
#      ii. only have peptide & MHC and in dict a
#      iii.has only one chain
#
cat ${output_dir}/chain.idx_all | \
  awk -v attr_idx=${output_dir}/attr.idx_all 'BEGIN{
      while(getline<attr_idx) {
        split($3,x,"[: }]");
        if ($1~/^tcr_pmhc_train_[0-9]+/) {
          if (x[1]<1.0)
            a[$1]=1;
        } else if ($1~/^tcr_pmhc_test_[0-9]+/) {  // We DO NOT test pMHCs
          a[$1]=1;
        }
      }
    }{
      has_P = 0;
      has_M = 0;
      for (i=2; i<=NF; ++i) {
        if ($i == "P")
          has_P +=1;
        else if ($i=="M")
          has_M +=1;
      }
      if (has_P==0 || (NF==3 && has_P>0 && has_M>0 && ($1 in a)) || NF<=2)
        print $0;
    }' > ${output_dir}/chain.idx_all_blacklist

for ((i=0;i<5;++i)); do
  start_idx=$((${i} * 1000000))
  end_idx=$((${i} * 1000000 + 1000000))

  # make chain.idx for fold_i
  cat ${output_dir}/chain.idx_all | \
    awk  -v start_idx=${start_idx} -v end_idx=${end_idx} '{
        if ($1~/^pmhc_(train|test)_[0-9]+$/) {
          print $0;
        } else if ($1~/^tcr_pmhc_(train|test)_[0-9]+$/) {
          n=split($1, x, "_");
          if (start_idx <= x[n] && x[n] < end_idx) {
            print $0;
          }
        }
      }' > ${output_dir}/chain.idx_tcr_pmhc_${i}

  # make attr.idx for train fold_i
  cat ${output_dir}/attr.idx_all | \
    awk  -v start_idx=${start_idx} -v end_idx=${end_idx} -v blacklist=${output_dir}/chain.idx_all_blacklist 'BEGIN{
        a["xxxxxxxx"] = 1;
        while(getline<blacklist)
          a[$1] = 1;
      }{
        if (!($1 in a)) {
          if ($1~/^pmhc_train_[0-9]+$/) {
            print $0;
          } else if ($1~/^tcr_pmhc_train_[0-9]+$/) {
            n = split($1, x, "_");
            if (start_idx <= x[n] && x[n] < end_idx) {
              print $0;
            }
          }
        }
      }' > ${output_dir}/train_attr.idx_tcr_pmhc_db_fold${i}

  cat ${CWD}/data/tcr_pmhc_db/attr.idx >> ${output_dir}/train_attr.idx_tcr_pmhc_db_fold${i}

  # make attr.idx for test fold_i
  cat ${output_dir}/attr.idx_all | \
    awk  -v start_idx=${start_idx} -v end_idx=${end_idx} -v blacklist=${output_dir}/chain.idx_all_blacklist 'BEGIN{
        a["xxxxxxxx"] = 1;
        while(getline<blacklist)
          a[$1] = 1;
      }{
        if (!($1 in a)) {
          if ($1~/^tcr_pmhc_test_[0-9]+$/) {
            n = split($1, x, "_");
            if (start_idx <= x[n] && x[n] < end_idx) {
              print $0;
            }
          }
        }
      }' > ${output_dir}/test_attr.idx_tcr_pmhc_db_fold${i}

  cat ${CWD}/data/tcr_pmhc_db/attr.idx ${output_dir}/attr.idx_all > ${output_dir}/attr.idx

  # calculate weights
  python ${CWD}/main.py tcr_pmhc_to_pmhc \
      --target_uri "${output_dir}?mapping_idx=mapping.idx_all&chain_idx=chain.idx_tcr_pmhc_${i}&attr_idx=train_attr.idx_tcr_pmhc_db_fold${i}" \
      --output "${output_dir}?mapping_idx=mapping.idx_tcr_pmhc_${i}&attr_idx=train_attr.idx_tcr_pmhc_db_fold${i}" \
      --pid_topk=1000

  #cat ${output_dir}/mapping.idx_tcr_pmhc_${i} | \
  #  awk -v mapping_idx=${output_dir}/mapping.idx_all 'BEGIN{
  #      while(getline<mapping_idx)
  #        a[$0] = 1;
  #    }{
  #      if (!($0 in a)) {
  #        print $0;
  #      }
  #    }' > ${output_dir}/mapping.idx_tcr_pmhc_delta_${i}
done

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
python ${CWD}/main.py align_peptide -o ${output_dir}/a3m --db ../db/tcr/tcr_pmhc_db_P.fa ${output_dir}/tcr_pmhc_P.fa -v ../custom_data/tcr_pmhc_db/fasta/*_P.fasta
for c in "P"; do
  find ${CWD}/data/tcr_pmhc_db/fasta -name "*_${c}.fasta" > ${output_dir}/tcr_pmhc_db_${c}
  cat ${output_dir}/tcr_pmhc_db_${c} | ${CWD}/bin/mapred -m "PIPELINE_UNIREF_MAX_HITS=1000000 PIPELINE_MGNIFY_MAX_HITS=1000000 PYTHONPATH=${CWD}/profold2 sh ${CWD}/scripts/run_pipeline.sh -o ${output_dir}/a3m" -c 10
done

# filter a3m with threshold=t
for t in "0.85" "0.90" "0.95"; do
  if [ -d ${output_dir}/var${t} ]; then
    rm -rf ${output_dir}/var${t}
  fi
  echo "filter a3m: t=${t}"
  cp -r ${output_dir}/a3m ${output_dir}/var${t}
  python a3m_filter.py --trim_gap -t ${t} -o ${output_dir}/var${t} ${CWD}/data/tcr_pmhc_db/fasta/*_M.fasta
done


# make sampling weights
for ((i=0;i<5;++i)); do
  for t in "0.85" "0.90" "0.95"; do
    python ${CWD}/main.py align_complex \
      -o ${output_dir}/complex_${i}_${t} \
      --db_uri "${output_dir}?mapping_idx=mapping.idx&attr_idx=train_attr.idx_tcr_pmhc_db_fold${i}" \
      --target_uri "${CWD}/data/tcr_pmhc_db?a3m_dir=${output_dir}/var${t}"
    mkdir -p ${output_dir}/cache_${i}_${t}
    find ${output_dir}/complex_${i}_${t} -name "*.pkl" | awk -F "/" -v d=${output_dir}/cache_${i}_${t} '{
        printf("cp %s %s/\n", $0, d);
      }' | bash -
    find ${output_dir}/complex_${i}_${t} -name "*.a3m" | python ${CWD}/a3m_name_list.py > ${output_dir}/complex_${i}_${t}_a3m_name_list
    cat ${output_dir}/complex_${i}_${t}_a3m_name_list | awk '{
        n = split($1,x,"/");
        k = x[n];
        a[k] += $3;
        b[k] += 1;
        c += $3;
        d += 1;
      } END {
        n = length(a);
        for (k in a) {
          printf("%s\t%f\t%f\t%d\t%f\t%d\n", k, a[k] * n/c, a[k], b[k], c, d);
        }
      }' | sed "s/.a3m//g" > ${output_dir}/tcr_pmhc_db_size_idx_${i}_${t}
    cat ${CWD}/data/tcr_pmhc_db/name.idx | awk -v value_idx=${output_dir}/tcr_pmhc_db_size_idx_${i}_${t} 'BEGIN{
        while (getline<value_idx) {
          a[$1] = $2;
        }
      } {
        w = 0;
        for (i = 1; i <= NF; ++i) {
          n = split($i, x, "_");

          k = x[1];
          for (j = 2; j < n; ++j) {
            k = k"_"x[j];
          }

          if (k in a) {
            w += a[k];
          } else {
            w += 1;
          }
        }
        w = w / NF;
        if (w < 0.1) {
          w = 0.1;
        } else if (w > 10.0) {
          w = 10.0;
        }
        print w"\t"$0;
      }' > ${output_dir}/model_weights_2022111102_tcr_pmhc_db_fold${i}_${t}

  done
done
