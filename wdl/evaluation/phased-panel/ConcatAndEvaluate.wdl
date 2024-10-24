version 1.0

import "./VcfdistAndOverlapMetricsEvaluation.wdl"
import "./ChromosomePhasedPanelCreationFromHiPhase.wdl"
import "../kage/LeaveOutEvaluation.wdl"
import "../../methods/phasing/HierarchicallyMergeVcfs.wdl"

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

workflow PhasedPanelEvaluation {    # TODO change name later, easier to share configs for now

    input {
        # common inputs
        File reference_fasta
        File reference_fasta_fai
        File reference_fasta_dict
        String output_prefix
        File? monitoring_script

        File hiphase_short_vcf_gz
        File hiphase_short_vcf_gz_tbi
        File hiphase_sv_vcf_gz
        File hiphase_sv_vcf_gz_tbi

        # per chromosome inputs
        Array[File] filter_and_concat_vcf_gzs
        Array[File] filter_and_concat_vcf_gz_tbis
        Array[File] before_shapeit4_collisionless_bcfs
        Array[File] before_shapeit4_collisionless_bcf_csis
        Array[File] phased_vcf_gzs
        Array[File] phased_vcf_gz_tbis
        Array[File] collisionless_bcfs
        Array[File] collisionless_bcf_csis
        Array[File] panel_vcf_gzs
        Array[File] panel_vcf_gz_tbis

        # inputs for LeaveOutEvaluation
        File case_reference_fasta
        File case_reference_fasta_fai
        File case_reference_dict
        Array[File] leave_out_genetic_maps
        File repeat_mask_bed
        File segmental_duplications_bed
        File simple_repeats_bed
        File challenging_medically_relevant_genes_bed
        Array[String] leave_out_chromosomes
        Array[String] leave_out_sample_names
        Int case_average_coverage
        Boolean do_pangenie
        Map[String, File] leave_out_crams
        String leave_out_docker
        String kage_docker
        String pangenie_docker
        Int? cpu_make_count_model
        RuntimeAttributes? leave_out_runtime_attributes
        RuntimeAttributes? leave_out_medium_runtime_attributes
        RuntimeAttributes? leave_out_large_runtime_attributes
        RuntimeAttributes? pangenie_runtime_attributes
        RuntimeAttributes? kage_count_kmers_runtime_attributes
        RuntimeAttributes? kage_genotype_runtime_attributes
        RuntimeAttributes? glimpse_case_chromosome_runtime_attributes
        RuntimeAttributes? glimpse_case_runtime_attributes
        RuntimeAttributes? calculate_metrics_runtime_attributes

        # inputs for HierarchicallyMergeVcfs
        Array[String] hierarchically_merge_regions   # bcftools regions, e.g. ["chr1,chr2,chr3", "chr4,chr5,chr6", ...]
        Int hierarchically_merge_batch_size
        String hierarchically_merge_docker

        # inputs for FixVariantCollisions
        File fix_variant_collisions_java
        Int operation
        String weight_tag
        Int is_weight_format_field

        # inputs for VcfdistAndOverlapMetricsEvaluation
        Array[String] vcfdist_samples
        File vcfdist_truth_vcf
        File vcfdist_truth_vcf_idx
        String evaluation_chromosomes_regions_arg
        File vcfdist_bed_file
        String? vcfdist_extra_args
        String overlap_metrics_docker

        String summarize_evaluations_docker
    }

    call ConcatVcfs as ConcatVcfsFilterAndConcatVcfs { input:
        vcfs = filter_and_concat_vcf_gzs,
        vcf_idxs = filter_and_concat_vcf_gz_tbis,
        prefix = output_prefix + ".filter_and_concat"
    }

    call ConcatVcfs as ConcatVcfsBeforeShapeit4FixVariantCollisions { input:
        vcfs = before_shapeit4_collisionless_bcfs,
        vcf_idxs = before_shapeit4_collisionless_bcf_csis,
        prefix = output_prefix + ".before_shapeit4_collisionless"
    }

    call ConcatVcfs as ConcatVcfsShapeit4 { input:
        vcfs = phased_vcf_gzs,
        vcf_idxs = phased_vcf_gz_tbis,
        prefix = output_prefix + ".phased"
    }

    call ConcatVcfs as ConcatVcfsFixVariantCollisions { input:
        vcfs = collisionless_bcfs,
        vcf_idxs = collisionless_bcf_csis,
        prefix = output_prefix + ".phased.collisionless"
    }

    call ConcatVcfs as ConcatVcfsPanGeniePanelCreation { input:
        vcfs = panel_vcf_gzs,
        vcf_idxs = panel_vcf_gz_tbis,
        prefix = output_prefix + ".panel"
    }

    call LeaveOutEvaluation.LeaveOutEvaluation { input:
        input_vcf_gz = ConcatVcfsPanGeniePanelCreation.vcf_gz,
        input_vcf_gz_tbi = ConcatVcfsPanGeniePanelCreation.vcf_gz_tbi,
        case_reference_fasta = case_reference_fasta,
        case_reference_fasta_fai = case_reference_fasta_fai,
        case_reference_dict = case_reference_dict,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        genetic_maps = leave_out_genetic_maps,
        repeat_mask_bed = repeat_mask_bed,
        segmental_duplications_bed = segmental_duplications_bed,
        simple_repeats_bed = simple_repeats_bed,
        challenging_medically_relevant_genes_bed = challenging_medically_relevant_genes_bed,
        output_prefix = output_prefix,
        chromosomes = leave_out_chromosomes,
        leave_out_sample_names = leave_out_sample_names,
        case_average_coverage = case_average_coverage,
        do_pangenie = do_pangenie,
        leave_out_crams = leave_out_crams,
        docker = leave_out_docker,
        kage_docker = kage_docker,
        pangenie_docker = pangenie_docker,
        monitoring_script = monitoring_script,
        cpu_make_count_model = cpu_make_count_model,
        runtime_attributes = leave_out_runtime_attributes,
        medium_runtime_attributes = leave_out_medium_runtime_attributes,
        large_runtime_attributes = leave_out_large_runtime_attributes,
        pangenie_runtime_attributes = pangenie_runtime_attributes,
        kage_count_kmers_runtime_attributes = kage_count_kmers_runtime_attributes,
        kage_genotype_runtime_attributes = kage_genotype_runtime_attributes,
        glimpse_case_chromosome_runtime_attributes = glimpse_case_chromosome_runtime_attributes,
        glimpse_case_runtime_attributes = glimpse_case_runtime_attributes,
        calculate_metrics_runtime_attributes = calculate_metrics_runtime_attributes
    }

    # merge GLIMPSE VCFs
    call HierarchicallyMergeVcfs.HierarchicallyMergeVcfs as GLIMPSEMergeAcrossSamples { input:
        vcf_gzs = LeaveOutEvaluation.glimpse_vcf_gzs,
        vcf_gz_tbis = LeaveOutEvaluation.glimpse_vcf_gz_tbis,
        regions = hierarchically_merge_regions,
        batch_size = hierarchically_merge_batch_size,
        output_prefix = output_prefix + ".glimpse.merged",
        docker = hierarchically_merge_docker,
        monitoring_script = monitoring_script
    }

    call ChromosomePhasedPanelCreationFromHiPhase.FixVariantCollisions as GenotypingFixVariantCollisions { input:
        phased_bcf = GLIMPSEMergeAcrossSamples.merged_vcf_gz,
        fix_variant_collisions_java = fix_variant_collisions_java,
        operation = operation,
        weight_tag = weight_tag,
        is_weight_format_field = is_weight_format_field,
        output_prefix = output_prefix + ".glimpse.merged.phased.collisionless"
    }

    # evaluate HiPhase short
    call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateHiPhaseShort { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        truth_vcf_idx = vcfdist_truth_vcf_idx,
        eval_vcf = hiphase_short_vcf_gz,
        eval_vcf_idx = hiphase_short_vcf_gz_tbi,
        region = evaluation_chromosomes_regions_arg,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        vcfdist_bed_file = vcfdist_bed_file,
        vcfdist_extra_args = vcfdist_extra_args,
        overlap_phase_tag = "PS",
        overlap_metrics_docker = overlap_metrics_docker
    }

    # evaluate HiPhase SV
    call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateHiPhaseSV { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        truth_vcf_idx = vcfdist_truth_vcf_idx,
        eval_vcf = hiphase_sv_vcf_gz,
        eval_vcf_idx = hiphase_sv_vcf_gz_tbi,
        region = evaluation_chromosomes_regions_arg,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        vcfdist_bed_file = vcfdist_bed_file,
        vcfdist_extra_args = vcfdist_extra_args,
        overlap_phase_tag = "PS",
        overlap_metrics_docker = overlap_metrics_docker
    }

    # evaluate filtered short + SV
    call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateFiltered { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        truth_vcf_idx = vcfdist_truth_vcf_idx,
        eval_vcf = ConcatVcfsFilterAndConcatVcfs.vcf_gz,
        eval_vcf_idx = ConcatVcfsFilterAndConcatVcfs.vcf_gz_tbi,
        region = evaluation_chromosomes_regions_arg,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        vcfdist_bed_file = vcfdist_bed_file,
        vcfdist_extra_args = vcfdist_extra_args,
        overlap_phase_tag = "PS",
        overlap_metrics_docker = overlap_metrics_docker
    }

    # evaluate before Shapeit4 collisionless short + SV
    call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateBeforeShapeit4FixVariantCollisions { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        truth_vcf_idx = vcfdist_truth_vcf_idx,
        eval_vcf = ConcatVcfsBeforeShapeit4FixVariantCollisions.vcf_gz,
        eval_vcf_idx = ConcatVcfsBeforeShapeit4FixVariantCollisions.vcf_gz_tbi,
        region = evaluation_chromosomes_regions_arg,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        vcfdist_bed_file = vcfdist_bed_file,
        vcfdist_extra_args = vcfdist_extra_args,
        overlap_phase_tag = "PS",
        overlap_metrics_docker = overlap_metrics_docker
    }

    # evaluate Shapeit4 short + SV
    call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateShapeit4 { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        truth_vcf_idx = vcfdist_truth_vcf_idx,
        eval_vcf = ConcatVcfsShapeit4.vcf_gz,
        eval_vcf_idx = ConcatVcfsShapeit4.vcf_gz_tbi,
        region = evaluation_chromosomes_regions_arg,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        vcfdist_bed_file = vcfdist_bed_file,
        vcfdist_extra_args = vcfdist_extra_args,
        overlap_phase_tag = "NONE",
        overlap_metrics_docker = overlap_metrics_docker
    }

    # evaluate collisionless short + SV
    call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateFixVariantCollisions { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        truth_vcf_idx = vcfdist_truth_vcf_idx,
        eval_vcf = ConcatVcfsFixVariantCollisions.vcf_gz,
        eval_vcf_idx = ConcatVcfsFixVariantCollisions.vcf_gz_tbi,
        region = evaluation_chromosomes_regions_arg,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        vcfdist_bed_file = vcfdist_bed_file,
        vcfdist_extra_args = vcfdist_extra_args,
        overlap_phase_tag = "NONE",
        overlap_metrics_docker = overlap_metrics_docker
    }

    # evaluate panel short + SV
    call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluatePanel { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        truth_vcf_idx = vcfdist_truth_vcf_idx,
        eval_vcf = ConcatVcfsPanGeniePanelCreation.vcf_gz,
        eval_vcf_idx = ConcatVcfsPanGeniePanelCreation.vcf_gz_tbi,
        region = evaluation_chromosomes_regions_arg,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        vcfdist_bed_file = vcfdist_bed_file,
        vcfdist_extra_args = vcfdist_extra_args,
        overlap_phase_tag = "NONE",
        overlap_metrics_docker = overlap_metrics_docker
    }

    # evaluate GLIMPSE
    call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateGenotyping { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        truth_vcf_idx = vcfdist_truth_vcf_idx,
        eval_vcf = GLIMPSEMergeAcrossSamples.merged_vcf_gz,
        eval_vcf_idx = GLIMPSEMergeAcrossSamples.merged_vcf_gz_tbi,
        region = evaluation_chromosomes_regions_arg,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        vcfdist_bed_file = vcfdist_bed_file,
        vcfdist_extra_args = vcfdist_extra_args,
        overlap_phase_tag = "NONE",
        overlap_metrics_docker = overlap_metrics_docker
    }

    # evaluate collisionless GLIMPSE
    call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateGenotypingFixVariantCollisions { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        truth_vcf_idx = vcfdist_truth_vcf_idx,
        eval_vcf = GenotypingFixVariantCollisions.collisionless_bcf,
        eval_vcf_idx = GenotypingFixVariantCollisions.collisionless_bcf_csi,
        region = evaluation_chromosomes_regions_arg,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        vcfdist_bed_file = vcfdist_bed_file,
        vcfdist_extra_args = vcfdist_extra_args,
        overlap_phase_tag = "NONE",
        overlap_metrics_docker = overlap_metrics_docker
    }

    # summarize GLIMPSE metrics vs. panel

    if (do_pangenie) {
        # merge PanGenie VCFs
        call HierarchicallyMergeVcfs.HierarchicallyMergeVcfs as PanGenieMergeAcrossSamples { input:
            vcf_gzs = select_all(LeaveOutEvaluation.pangenie_vcf_gzs),
            vcf_gz_tbis = select_all(LeaveOutEvaluation.pangenie_vcf_gz_tbis),
            regions = hierarchically_merge_regions,
            batch_size = hierarchically_merge_batch_size,
            output_prefix = output_prefix + ".pangenie.merged",
            docker = hierarchically_merge_docker,
            monitoring_script = monitoring_script
        }

        # summarize PanGenie metrics vs. panel

        # evaluate PanGenie
        call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluatePanGenie { input:
            samples = vcfdist_samples,
            truth_vcf = vcfdist_truth_vcf,
            truth_vcf_idx = vcfdist_truth_vcf_idx,
            eval_vcf = PanGenieMergeAcrossSamples.merged_vcf_gz,
            eval_vcf_idx = PanGenieMergeAcrossSamples.merged_vcf_gz_tbi,
            region = evaluation_chromosomes_regions_arg,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            vcfdist_bed_file = vcfdist_bed_file,
            vcfdist_extra_args = vcfdist_extra_args,
            overlap_phase_tag = "NONE",
            overlap_metrics_docker = overlap_metrics_docker
        }
    }

    Array[String] labels_per_vcf = if do_pangenie then ["HiPhaseShort", "HiPhaseSV", "ConcatAndFiltered", "BeforeShapeit4FixVariantCollisions", "Shapeit4", "FixVariantCollisions", "Panel", "Genotyping", "GenotypingFixVariantCollisions", "PanGenie"] else ["HiPhaseShort", "HiPhaseSV", "ConcatAndFiltered", "BeforeShapeit4FixVariantCollisions", "Shapeit4", "FixVariantCollisions", "Panel", "Genotyping", "GenotypingFixVariantCollisions"]
    call SummarizeEvaluations { input:
        labels_per_vcf = labels_per_vcf,
        vcfdist_outputs_per_vcf_and_sample = select_all([
            EvaluateHiPhaseShort.vcfdist_outputs_per_sample,
            EvaluateHiPhaseSV.vcfdist_outputs_per_sample,
            EvaluateFiltered.vcfdist_outputs_per_sample,
            EvaluateBeforeShapeit4FixVariantCollisions.vcfdist_outputs_per_sample,
            EvaluateShapeit4.vcfdist_outputs_per_sample,
            EvaluateFixVariantCollisions.vcfdist_outputs_per_sample,
            EvaluatePanel.vcfdist_outputs_per_sample,
            EvaluateGenotyping.vcfdist_outputs_per_sample,
            EvaluateGenotypingFixVariantCollisions.vcfdist_outputs_per_sample,
            EvaluatePanGenie.vcfdist_outputs_per_sample
        ]),
        overlap_metrics_outputs_per_vcf = select_all([
            EvaluateHiPhaseShort.overlap_metrics_outputs,
            EvaluateHiPhaseSV.overlap_metrics_outputs,
            EvaluateFiltered.overlap_metrics_outputs,
            EvaluateBeforeShapeit4FixVariantCollisions.overlap_metrics_outputs,
            EvaluateShapeit4.overlap_metrics_outputs,
            EvaluateFixVariantCollisions.overlap_metrics_outputs,
            EvaluatePanel.overlap_metrics_outputs,
            EvaluateGenotyping.overlap_metrics_outputs,
            EvaluateGenotypingFixVariantCollisions.overlap_metrics_outputs,
            EvaluatePanGenie.overlap_metrics_outputs
        ]),
        docker = summarize_evaluations_docker
    }

    output {
        File panel_vcf_gz = ConcatVcfsPanGeniePanelCreation.vcf_gz
        File panel_vcf_gz_tbi = ConcatVcfsPanGeniePanelCreation.vcf_gz_tbi

        File evaluation_summary_tsv = SummarizeEvaluations.evaluation_summary_tsv
    }
}

task ConcatVcfs {

    input {
        Array[File] vcfs
        Array[File]? vcf_idxs
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(vcfs, "GB")) + 1

    command <<<
        set -euxo pipefail
        if ! ~{defined(vcf_idxs)}; then
            for ff in ~{sep=' ' vcfs}; do bcftools index $ff; done
        fi
        bcftools concat ~{sep=" " vcfs} -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File vcf_gz = "~{prefix}.vcf.gz"
        File vcf_gz_tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FixVariantCollisions {

    input {
        File phased_bcf                     # biallelic
        File fix_variant_collisions_java
        Int operation = 1                   # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "UNIT_WEIGHT"   # ID of the weight field; if this field is not found, all weights are set to one; weights are assumed to be non-negative
        Int is_weight_format_field = 0      # given a VCF record in a sample, assign it a weight encoded in the sample column (1) or in the INFO field (0)
        String output_prefix
    }

    command <<<
        set -euxo pipefail

        # convert bcf to vcf.gz
        bcftools view ~{phased_bcf} -Oz -o phased.vcf.gz
        bcftools index -t phased.vcf.gz

        java ~{fix_variant_collisions_java} \
            phased.vcf.gz \
            ~{operation} \
            ~{weight_tag} \
            ~{is_weight_format_field} \
            collisionless.vcf \
            windows.txt \
            histogram.txt \
            null                            # do not output figures

        # replace all missing alleles (correctly) emitted with reference alleles, since this is expected by PanGenie panel-creation script
        bcftools view collisionless.vcf | \
            sed -e 's/\.|0/0|0/g' | sed -e 's/0|\./0|0/g' | sed -e 's/\.|1/0|1/g' | sed -e 's/1|\./1|0/g' | sed -e 's/\.|\./0|0/g' | \
            bcftools view -Oz -o ~{output_prefix}.vcf.gz
        # index and convert via vcf.gz to avoid errors from missing header lines
        bcftools index -t ~{output_prefix}.vcf.gz
        bcftools view ~{output_prefix}.vcf.gz -Ob -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File collisionless_bcf = "~{output_prefix}.bcf"
        File collisionless_bcf_csi = "~{output_prefix}.bcf.csi"
        File windows = "windows.txt"
        File histogram = "histogram.txt"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "16 GiB"
        disks: "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-gatk/gatk:4.6.0.0"
    }
}

task SummarizeEvaluations {
    input {
        Array[String] labels_per_vcf
        Array[Array[VcfdistOutputs]] vcfdist_outputs_per_vcf_and_sample
        Array[OverlapMetricsOutputs] overlap_metrics_outputs_per_vcf

        String docker
        Int disk_size_gb = 50
        Int boot_disk_size_gb = 15
        Int mem_gb = 8
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        python - --labels_per_vcf_txt ~{write_lines(labels_per_vcf)} \
                 --vcfdist_outputs_per_vcf_and_sample_json ~{write_json(vcfdist_outputs_per_vcf_and_sample)} \
                 --overlap_metrics_outputs_per_vcf_json ~{write_json(overlap_metrics_outputs_per_vcf)} \
                 <<-'EOF'
        import argparse
        import json
        import pandas as pd

        def summarize(labels_per_vcf_txt,
                      vcfdist_outputs_per_vcf_and_sample_json,
                      overlap_metrics_outputs_per_vcf_json):
            with open(labels_per_vcf_txt) as f:
                labels = f.read().splitlines()

            with open(vcfdist_outputs_per_vcf_and_sample_json) as f:
                vcfdist_outputs_per_vcf_and_sample = json.load(f)

            with open(overlap_metrics_outputs_per_vcf_json) as f:
                overlap_metrics_outputs_per_vcf = json.load(f)

            summary_dict = {}
            for i, label in enumerate(labels):
                summary_dict[label] = {}
                summary_dict[label]['NUM_VCFDIST_SAMPLES'] = len(vcfdist_outputs_per_vcf_and_sample[i])
                summary_dict[label].update(summarize_vcfdist_outputs_over_samples(vcfdist_outputs_per_vcf_and_sample[i]))
                if overlap_metrics_outputs_per_vcf:
                    summary_dict[label].update(summarize_overlap_metrics_outputs(overlap_metrics_outputs_per_vcf[i]))

            pd.DataFrame.from_dict(summary_dict, orient='index').to_csv('evaluation_summary.tsv', sep='\t', float_format="%.4f")

        def summarize_vcfdist_outputs_over_samples(vcfdist_outputs_per_sample):
            precision_recall_metrics_dict = {}
            for s, vcfdist_outputs in enumerate(vcfdist_outputs_per_sample):
                precision_recall_metrics_dict[s] = {}
                pr_metrics_df = pd.read_csv(vcfdist_outputs['precision_recall_summary_tsv'], sep='\t', index_col=[0, 1])
                for var_type in ['SNP', 'INDEL', 'SV']:
                    var_type_metrics_dict = pr_metrics_df.loc[var_type, 'NONE'][['TRUTH_TP', 'QUERY_TP', 'TRUTH_FN', 'QUERY_FP', 'PREC', 'RECALL', 'F1_SCORE']].add_prefix(f'{var_type}_').to_dict()
                    precision_recall_metrics_dict[s].update(var_type_metrics_dict)
            return pd.DataFrame.from_dict(precision_recall_metrics_dict, orient='index').mean(axis=0)


        def summarize_overlap_metrics_outputs(overlap_metrics_outputs):
            return pd.read_csv(overlap_metrics_outputs['metrics_tsv'], sep='\t').iloc[0].to_dict()

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--labels_per_vcf_txt',
                                type=str)

            parser.add_argument('--vcfdist_outputs_per_vcf_and_sample_json',
                                type=str)

            parser.add_argument('--overlap_metrics_outputs_per_vcf_json',
                                type=str)

            args = parser.parse_args()

            summarize(args.labels_per_vcf_txt,
                      args.vcfdist_outputs_per_vcf_and_sample_json,
                      args.overlap_metrics_outputs_per_vcf_json)

        if __name__ == '__main__':
            main()
        EOF
    >>>

    output {
        File evaluation_summary_tsv = "evaluation_summary.tsv"
    }

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        bootDiskSizeGb: boot_disk_size_gb
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}
