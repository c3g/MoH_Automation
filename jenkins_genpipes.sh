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

##################################################
# Initialization
chunk_submit=false
module purge
module load mugqic/python/3.10.2


export MUGQIC_PIPELINES_HOME=${path}/genpipes_moh/genpipes

cd $path

while IFS=, read -r readset_file pair_file; do
  timestamp=$(date +%Y-%m-%dT%H.%M.%S)
  sample=$(awk '{print $1}' $readset_file)
  echo $sample
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
      genpipes_file=rnaseq_light_${sample}_${timestamp}.sh
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
      genpipes_file=rnaseq_cancer_${sample}_${timestamp}.sh
  elif test $pipeline == tumor_pair; then
      # tumor_pair
      if test $protocol == ensemble; then
        steps="5-13,15-38"
        custom_ini="TP_ensemble.custom.ini"
        genpipes_file="tumor_pair_ensemble_${sample}_${timestamp}.sh"
      elif test $protocol == sv; then
        steps="12-16"
        custom_ini="TP_sv.custom.ini"
        genpipes_file="tumor_pair_sv_${sample}_${timestamp}.sh"
      fi
      echo "$MUGQIC_PIPELINES_HOME/pipelines/tumor_pair/tumor_pair.py \
-t $protocol \
-s $steps \
-c $MUGQIC_PIPELINES_HOME/pipelines/tumor_pair/tumor_pair.base.ini ${beluga_ini}\
$MUGQIC_PIPELINES_HOME/pipelines/common_ini/Homo_sapiens.GRCh38.ini \
$custom_ini \
-j $scheduler \
-r $readset_file \
-p $pair_file \
-g $genpipes_file \
--json-pt"
      chunk_submit=true
  fi
  # Chunking & Submission
  if test $chunk_submit == true && test -f "$genpipes_file"; then
    echo '-> Chunking...'
    $MUGQIC_PIPELINES_HOME/utils/chunk_genpipes.sh $genpipes_file job_chunks 20
    echo '-> Submitting...'
    $MUGQIC_PIPELINES_HOME/utils/submit_genpipes -n 800 job_chunks
  fi
done < ${path}/${input_file}