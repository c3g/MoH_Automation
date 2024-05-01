#!/usr/bin/env python3

import argparse
import json
import re
import os
import logging

logging.basicConfig(format='%(levelname)s: %(asctime)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(prog=os.path.basename(__file__), description="Adds appropriate deliverable tags on files and metrics + metrics thresholds.")
    parser.add_argument('-i', '--input', required=True, help="Input json produced by GenPipes processing.")
    parser.add_argument('-o', '--output', required=False, help="Output json filename (Default: <input_filename>_tagged.json).")
    args = parser.parse_args()

    if not args.output:
        output = os.path.basename(args.input).replace(".json", "_tagged.json")
    else:
        output = args.output

    mark_deliverables(args.input, output)


def mark_deliverables(input_json, output_json):
    # dict for files and metrics to be labelled as deliverables based on job_name (regex)
    # TumorPair.ensemble GenPipes
    deliverables_tumorpair_ensemble = {
        "merge_gatk_variant_annotator.germline.*": {
            "files": [
                ".+.ensemble.germline.vt.annot.vcf.gz",
                ".+.ensemble.germline.vt.annot.vcf.gz.tbi"
                ]
            },
        "merge_gatk_variant_annotator.somatic.*": {
            "files": [
                ".+.ensemble.somatic.vt.annot.vcf.gz",
                ".+.ensemble.somatic.vt.annot.vcf.gz.tbi"
                ]
            },
        "merge_filter_mutect2.*": {
            "files": [
                ".+.mutect2.somatic.vt.vcf.gz",
                ".+.mutect2.somatic.vt.vcf.gz.tbi"
                ]
            },
        "strelka2_paired_germline.filter.*": {
            "files": [
                ".+.strelka2.germline.vt.vcf.gz"
                ]
            },
        "strelka2_paired_somatic.filter.*": {
            "files": [
                ".+.strelka2.somatic.vt.vcf.gz"
                ]
            },
        "merge_filter_paired_vardict.*": {
            "files": [
                ".+.vardict.germline.vt.vcf.gz",
                ".+.vardict.germline.vt.vcf.gz.tbi",
                ".+.vardict.somatic.vt.vcf.gz",
                ".+.vardict.somatic.vt.vcf.gz.tbi"
                ]
            },
        "merge_varscan2.*": {
            "files": [
                ".+.varscan2.germline.vt.vcf.gz",
                ".+.varscan2.germline.vt.vcf.gz.tbi",
                ".+.varscan2.somatic.vt.vcf.gz",
                ".+.varscan2.somatic.vt.vcf.gz.tbi"
                ]
            },
        "cnvkit_batch.call.*": {
            "files": [
                ".+.cnvkit.vcf.gz",
                ".+.cnvkit.vcf.gz.tbi"
                ]
            },
        "gatk_print_reads.*": {
            "files": [
                ".+.sorted.dup.recal.bam",
                ".+.sorted.dup.recal.bam.bai",
                ".+.sorted.dup.recal.bam.md5"
                ]
            },
        "multiqc.*": {
            "files": [
                ".+.multiqc.html"
                ]
            },
        "report_pcgr.*": {
            "files": [
                ".+.pcgr_acmg.grch38.flexdb.html",
                ".+.pcgr_acmg.grch38.maf",
                ".+.pcgr_acmg.grch38.snvs_indels.tiers.tsv",
                ".+.pcgr_acmg.grch38.cna_segments.tsv.gz"
                ]
            },
        "conpair_concordance_contamination.*": {
            "metrics": [
                "concordance",
                "contamination"
                ]
        },
        "picard_collect_multiple_metrics.*": {
            "metrics": [
                "bases_over_q30_percent"
                ]
            },
        "dna_sample_qualimap.*": {
            "metrics": [
                "median_insert_size",
                "mean_insert_size",
                "dedup_coverage",
                "aligned_reads_count"
                ]
            },
        "purple.purity.*": {
            "files": [
                ".+.driver.catalog.somatic.tsv",
                ".+.driver.catalog.germline.tsv",
                ],
            "metrics": [
                "purity"
                ]
            },
    }
    # TumorPair.sv GenPipes
    deliverables_tumorpair_sv = {
        "gridss_paired_somatic.*": {
            "files": [
                ".+.gridss.vcf.gz"
                ]
            },
        "gripss_filter.somatic.*": {
            "files": [
                ".+.gripss.filtered.somatic.vcf.gz"
                ]
            },
        "gripss_filter.germline.*": {
            "files": [
                ".+.gripss.filtered.germline.vcf.gz"
                ]
            },
        "purple.purity.*": {
            "files": [
                ".+.circos.png"
                ]
            },
        "linx_annotations_germline.*": {
            "files": [
                ".+.linx.germline.clusters.tsv",
                ".+.linx.germline.disruption.tsv",
                ".+.linx.germline.driver.catalog.tsv",
                ".+.linx.germline.links.tsv",
                ".+.linx.germline.svs.tsv"
                ]
            },
        "linx_annotations_somatic.*": {
            "files": [
                ".+.linx.breakend.tsv",
                ".+.linx.clusters.tsv",
                ".+.linx.driver.catalog.tsv",
                ".+.linx.drivers.tsv",
                ".+.linx.fusion.tsv",
                ".+.linx.links.tsv",
                ".+.linx.svs.tsv",
                ".+.linx.vis_copy_number.tsv",
                ".+.linx.vis_fusion.tsv",
                ".+.linx.vis_gene_exon.tsv",
                ".+.linx.vis_protein_domain.tsv",
                ".+.linx.vis_segments.tsv",
                ".+.linx.vis_sv_data.tsv"
                ]
            },
        "linx_plot.*": {
            "files": [
                r".+.cluster-\d\d\d.sv\d\d.\d\d\d.png"
            ]
        }
    }
    # RnaSeq.cancer GenPipes
    deliverables_rnaseq_cancer = {
        "filter_gatk.*": {
            "files": [
                ".+.hc.vt.annot.vcf.gz",
                ".+.hc.vt.annot.flt.vcf.gz",
                ".+.hc.vt.annot.flt.vcf.gz.tbi",
                ".+.hc.vt.annot.flt.vcf.gz.md5"
                ]
            },
        "gatk_print_reads.*": {
            "files": [
                ".+.sorted.mdup.split.recal.bam",
                ".+.sorted.mdup.split.recal.bam.bai",
                ".+.sorted.mdup.split.recal.bam.md5"
                ]
            },
        "report_pcgr.*": {
            "files": [
                ".+.pcgr_acmg.grch38.flexdb.html",
                ".+.pcgr_acmg.grch38.maf",
                ".+.pcgr_acmg.grch38.snvs_indels.tiers.tsv"
                ]
            },
        "run_annoFuse.*": {
            "files": [
                ".+.putative_driver_fusions.tsv"
                ]
            },
        "picard_rna_metrics.*": {
            "metrics": [
                "aligned_reads_ratio"
                ]
            },
        "rnaseqc2.*": {
            "metrics": [
                "expression_profiling_efficiency",
                "rrna_rate"
                ]
            },
        "multiqc.*": {
            "files": [
                "multiqc_*.html",
                ".+.multiqc.html"
                ]
            },
    }
    # RnaSeqLight GenPipes
    deliverables_rnaseqlight = {
        "kallisto.*": {
            "files": [
                "abundance_transcripts.tsv",
                "abundance_genes.tsv"
                ],
            "metrics": [
                "mean_insert_size",
                "median_insert_size"
                ]
            },
    }

    with open(input_json, 'r', encoding='utf-8') as j_file:
        current_json_hash = json.load(j_file)

    if current_json_hash["operation_name"] == "GenPipes_TumorPair.ensemble":
        json_hash_out = iterate_json(current_json_hash, deliverables_tumorpair_ensemble)
    elif current_json_hash["operation_name"] == "GenPipes_TumorPair.sv":
        json_hash_out = iterate_json(current_json_hash, deliverables_tumorpair_sv)
    elif current_json_hash["operation_name"] == "GenPipes_RnaSeq.cancer":
        json_hash_out = iterate_json(current_json_hash, deliverables_rnaseq_cancer)
    elif current_json_hash["operation_name"] == "GenPipes_RnaSeqLight":
        json_hash_out = iterate_json(current_json_hash, deliverables_rnaseqlight)

    if json_hash_out:
        with open(output_json, 'w', encoding='utf-8') as j_file:
            json.dump(json_hash_out, j_file, ensure_ascii=False, indent=4)


def iterate_json(current_json_hash, deliverables_metrics_json):
    for sample in current_json_hash["sample"]:
        for readset in sample["readset"]:
            for job in readset["job"]:
                key_to_use = None
                for key in deliverables_metrics_json:
                    result = re.search(key, job["job_name"])
                    if result:
                        key_to_use = key
                if key_to_use:
                    try:
                        for file in job["file"]:
                            result = re.search('|'.join(deliverables_metrics_json[key_to_use]["files"]), file["file_name"])
                            if result:
                                file["file_deliverable"] = True
                    except KeyError:
                        pass
                    try:
                        for metric in job["metric"]:
                            result = re.search('|'.join(deliverables_metrics_json[key_to_use]["metrics"]), metric["metric_name"])
                            if result:
                                metric["metric_deliverable"] = True
                            # Checking thresholds for flagging metrics
                            metric_name = metric["metric_name"]
                            metric_value = metric["metric_value"]
                            metric["metric_flag"] = "PASS"
                            # Concordance
                            if metric_name == "concordance" and float(metric_value) < 99:
                                metric["metric_flag"] = "FAILED"
                            # Contamination
                            if metric_name == "contamination":
                                if float(metric_value) > 0.5:
                                    metric["metric_flag"] = "WARNING"
                                elif float(metric_value) > 0.05:
                                    metric["metric_flag"] = "FAILED"
                            # Median Insert Size
                            if metric_name == "median_insert_size":
                                if int(float(metric_value)) < 300:
                                    metric["metric_flag"] = "WARNING"
                                elif int(float(metric_value)) < 150:
                                    metric["metric_flag"] = "FAILED"
                            # Dedup Coverage DN
                            if metric_name == "dedup_coverage" and float(metric_value) < 30 and sample["sample_name"].endswith("DN"):
                                metric["metric_flag"] = "FAILED"
                            # Dedup Coverage DT
                            if metric_name == "dedup_coverage" and float(metric_value) < 80 and sample["sample_name"].endswith("DT"):
                                metric["metric_flag"] = "FAILED"
                            # Purity
                            if metric_name == "purity" and int(metric_value) < 30:
                                metric["metric_flag"] = "FAILED"
                            # rRNA rate
                            if metric_name == "rrna_rate":
                                if float(metric_value) > 0.35:
                                    metric["metric_flag"] = "FAILED"
                                elif float(metric_value) > 0.1:
                                    metric["metric_flag"] = "WARNING"
                    except KeyError:
                        pass
    return current_json_hash


if __name__ == '__main__':
    main()
