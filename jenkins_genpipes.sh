#!/bin/bash
set -eu -o pipefail


usage() {
  echo "script usage: jenkins_genpipes.sh -h [-c cluster] [-p pipeline] [-t protocol] [-i input_file]"
  echo "Usage:"
  echo " -h                               Display this help message."
  echo " -c <cluster>                     Cluster name on which the analysis will be run."
  echo " -p <pipeline>                    Pipeline name to be used for the analysis."
  echo " -t <protocol>                    Protocol to be used for the analysis. (Optional)"
  echo " -i <input_file>                  Path to Input File to be used for the analysis. This file is a csv file with 1st
                                          column being to the readset file and the 2nd (optional) column being the pair file."
  exit 1
  }

while getopts 'hc:p::t:r:i:' OPTION; do
  case "$OPTION" in
    c)
      cluster="$OPTARG"
      if [[ $cluster == beluga ]]; then
        beluga_ini="$MUGQIC_PIPELINES_HOME/pipelines/common_ini/beluga.ini"
        path="/lustre03/project/6007512/C3G/projects/MOH_PROCESSING/MAIN"
        scheduler="slurm"
        export MUGQIC_INSTALL_HOME_DEV=/project/6007512/C3G/analyste_dev
        
      elif [[ $cluster == abacus ]]; then
        beluga_ini=""
        path="/lb/project/mugqic/projects/MOH/MAIN"
        scheduler="pbs"
        export MUGQIC_INSTALL_HOME_DEV=/lb/project/mugqic/analyste_dev
        export MUGQIC_INSTALL_HOME_PRIVATE=/lb/project/mugqic/analyste_private
      fi
      echo "cluster: $cluster"
      ;;
    p)
      pipeline="$OPTARG"
      echo "pipeline: $pipeline"
      ;;
    t)
      protocol="$OPTARG"
      echo "protocol: $protocol"
      ;;
    i)
      input_file="$OPTARG"
      echo "input_file: $input_file"
      ;;
    h)
      usage
      ;;
    ?)
      usage
      ;;
  esac
done

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
export PORTAL_OUTPUT_DIR=$MUGQIC_INSTALL_HOME_DEV/portal_out_dir
module use $MUGQIC_INSTALL_HOME/modulefiles $MUGQIC_INSTALL_HOME_DEV/modulefiles
export JOB_MAIL=c3g-processing@fakeemail.ca

##################################################
# Initialization
chunk_submit=false
module purge
module load mugqic/python/3.10.2


export MUGQIC_PIPELINES_HOME=${path}/genpipes_moh/genpipes

cd $path

while IFS=, read -r readset_file pair_file; do
  timestamp=$(date +%Y-%m-%dT%H.%M.%S)
  patient=$(awk 'NR==2, match($1, /^((MoHQ-(JG|HM|CM|GC|MU|MR|XX)-\w+)-\w+)/) {print substr($1, RSTART, RLENGTH)}' $readset_file)
  sample=$(awk 'NR>1{print $1}' $readset_file)
  # echo $sample
  echo '-> Running GenPipes...'
  # GenPipes call
  if test $pipeline == rnaseq_light; then
      # rnaseq_light
      echo "$MUGQIC_PIPELINES_HOME/pipelines/rnaseq_light/rnaseq_light.py \
-s 1-4 \
-c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq_light/rnaseq_light.base.ini ${beluga_ini}\
$MUGQIC_PIPELINES_HOME/pipelines/common_ini/Homo_sapiens.GRCh38.ini \
RNA_light.custom.ini \
-j $scheduler \
-r $readset_file \
-g rnaseq_light.sh \
--json-pt"
      chunk_submit=true
      genpipes_file=rnaseq_light_${patient}_${timestamp}.sh
  elif test $pipeline == rnaseq; then
      # rnaseq
      echo "$MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.py \
-t $protocol \
-s 1-8,11-28 \
-c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnase.base.ini ${beluga_ini}\
$MUGQIC_PIPELINES_HOME/pipelines/common_ini/Homo_sapiens.GRCh38.ini \
Custom_ini/tumor_rna.moh.ini \
-j $scheduler \
-r $readset_file \
-g rnaseq_cancer.sh \
--json-pt"
      chunk_submit=true
      genpipes_file=rnaseq_cancer_${patient}_${timestamp}.sh
  elif test $pipeline == tumor_pair; then
      # tumor_pair
      if test $protocol == ensemble; then
        steps="5-13,15-38"
        custom_ini="TP_ensemble.custom.ini"
        genpipes_file="tumor_pair_ensemble_${patient}_${timestamp}.sh"
      elif test $protocol == sv; then
        steps="12-16"
        custom_ini="TP_sv.custom.ini"
        genpipes_file="tumor_pair_sv_${patient}_${timestamp}.sh"
      fi
      $MUGQIC_PIPELINES_HOME/pipelines/tumor_pair/tumor_pair.py \
-t $protocol \
-s $steps \
-c $MUGQIC_PIPELINES_HOME/pipelines/tumor_pair/tumor_pair.base.ini \
$MUGQIC_PIPELINES_HOME/pipelines/tumor_pair/tumor_pair.extras.ini ${beluga_ini}\
$MUGQIC_PIPELINES_HOME/pipelines/common_ini/Homo_sapiens.GRCh38.ini \
$custom_ini \
-j $scheduler \
-r $readset_file \
-p $pair_file \
-g $genpipes_file \
--json-pt
      chunk_submit=true
  fi
  # Chunking & Submission
  if test $chunk_submit == true && test -f "$genpipes_file"; then
    today=$(date +%Y-%m-%dT)
    chmod 664 *.${protocol}.${today}*.config.trace.ini
    chmod 774 $genpipes_file
    echo '-> Chunking...'
    $MUGQIC_PIPELINES_HOME/utils/chunk_genpipes.sh $genpipes_file ${patient}_${timestamp}_chunks
    exit
    chmod 775 ${patient}_${timestamp}_chunks
    chmod 664 ${patient}_${timestamp}_chunks/*
    echo '-> Submitting...'
    $MUGQIC_PIPELINES_HOME/utils/submit_genpipes ${patient}_${timestamp}_chunks
    cat /dev/null > ${patient}_${timestamp}.txt
    (sleep 1 && submit_genpipes ${patient}_${timestamp}_chunks >> ${patient}_${timestamp}.txt 2>&1) & echo -n "PID: " >> ${patient}_${timestamp}.txt
    echo $! >> ${patient}_${timestamp}.txt
    echo "PATIENT: ${patient}" >> ${patient}_${timestamp}.txt
    echo "LOG:" >> ${patient}_${timestamp}.txt
  fi
done < ${path}/${input_file}