version 1.0

import "../../methods/phasing/PhysicalAndStatisticalPhasing.wdl"
import "../../methods/pangenie/PanGeniePanelCreation.wdl"
import "./VcfdistAndOverlapMetricsEvaluation.wdl"

workflow PhasedPanelEvaluation {

    input {
        # common inputs
        File reference_fasta
        File reference_fasta_fai
        String region
        String gcs_out_root_dir
        String output_prefix

        # inputs for PhysicalAndStatisticalPhasing
        Array[File] sample_bams
        Array[File] sample_bais
        File joint_short_vcf
        File joint_short_vcf_tbi
        File joint_sv_vcf
        File joint_sv_vcf_tbi
        File genetic_mapping_tsv_for_shapeit4
        String chromosome
        Int num_t

        # inputs for PanGeniePanelCreation
        File prepare_vcf_script
        File add_ids_script
        File merge_vcfs_script
        String panel_creation_docker
        File? panel_creation_monitoring_script

        # inputs for VcfdistAndOverlapMetricsEvaluation
        Array[String] vcfdist_samples
        File vcfdist_truth_vcf
        File vcfdist_bed_file
        String? vcfdist_extra_args
        String overlap_metrics_docker
    }

    call PhysicalAndStatisticalPhasing.PhysicalAndStatisticalPhasing { input:
        sample_bams = sample_bams,
        sample_bais = sample_bais,
        joint_short_vcf = joint_short_vcf,
        joint_short_vcf_tbi = joint_short_vcf_tbi,
        joint_sv_vcf = joint_sv_vcf,
        joint_sv_vcf_tbi = joint_sv_vcf_tbi,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        genetic_mapping_tsv_for_shapeit4 = genetic_mapping_tsv_for_shapeit4,
        chromosome = chromosome,
        region = region,
        prefix = output_prefix,
        gcs_out_root_dir = gcs_out_root_dir + "/Phasing",
        num_t = num_t
    }

    call PanGeniePanelCreation.PanGeniePanelCreation { input:
        phased_bcf = PhysicalAndStatisticalPhasing.phased_bcf,
        reference_fasta = reference_fasta,
        prepare_vcf_script = prepare_vcf_script,
        add_ids_script = add_ids_script,
        merge_vcfs_script = merge_vcfs_script,
        output_prefix = output_prefix,
        docker = panel_creation_docker,
        monitoring_script = panel_creation_monitoring_script
    }

    # evaluate HiPhase short
    call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateHiPhaseShort { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        eval_vcf = PhysicalAndStatisticalPhasing.hiphase_short_vcf,
        region = region,
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
        eval_vcf = PhysicalAndStatisticalPhasing.hiphase_sv_vcf,
        region = region,
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
        eval_vcf = PhysicalAndStatisticalPhasing.filtered_vcf,
        region = region,
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
        eval_vcf = PhysicalAndStatisticalPhasing.phased_bcf,
        region = region,
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
        eval_vcf = PanGeniePanelCreation.panel_vcf_gz,
        region = region,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        vcfdist_bed_file = vcfdist_bed_file,
        vcfdist_extra_args = vcfdist_extra_args,
        overlap_phase_tag = "NONE",
        overlap_metrics_docker = overlap_metrics_docker
    }

    output {
    }
}


