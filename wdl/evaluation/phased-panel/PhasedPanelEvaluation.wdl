version 1.0

import "../../methods/phasing/PhysicalAndStatisticalPhasing.wdl"
import "../../methods/pangenie/PanGeniePanelCreation.wdl"
import "../vcfdist/VcfdistEvaluation.wdl"

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

        # inputs for VcfdistEvaluation
        Array[String] vcfdist_samples
        File vcfdist_truth_vcf
        File vcfdist_bed_file
        String vcfdist_docker = "timd1/vcfdist:v2.5.3"
        Int vcfdist_verbosity = 1
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
    call VcfdistEvaluation.VcfdistEvaluation as VcfdistEvaluationHiPhaseShort { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        eval_vcf = PhysicalAndStatisticalPhasing.hiphase_short_vcf,
        bed_file = vcfdist_bed_file,
        region = region,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        docker = vcfdist_docker,
        verbosity = vcfdist_verbosity
    }

    # evaluate HiPhase SV
    call VcfdistEvaluation.VcfdistEvaluation as VcfdistEvaluationHiPhaseSV { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        eval_vcf = PhysicalAndStatisticalPhasing.hiphase_sv_vcf,
        bed_file = vcfdist_bed_file,
        region = region,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        docker = vcfdist_docker,
        verbosity = vcfdist_verbosity
    }

    # evaluate filtered short + SV
    call VcfdistEvaluation.VcfdistEvaluation as VcfdistEvaluationFiltered { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        eval_vcf = PhysicalAndStatisticalPhasing.filtered_vcf,
        bed_file = vcfdist_bed_file,
        region = region,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        docker = vcfdist_docker,
        verbosity = vcfdist_verbosity
    }

    # evaluate Shapeit4 short + SV
    call VcfdistEvaluation.VcfdistEvaluation as VcfdistEvaluationShapeit4 { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        eval_vcf = PhysicalAndStatisticalPhasing.phased_bcf,
        bed_file = vcfdist_bed_file,
        region = region,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        docker = vcfdist_docker,
        verbosity = vcfdist_verbosity
    }

    # evaluate panel short + SV
    call VcfdistEvaluation.VcfdistEvaluation as VcfdistEvaluationPanel { input:
        samples = vcfdist_samples,
        truth_vcf = vcfdist_truth_vcf,
        eval_vcf = PanGeniePanelCreation.panel_vcf_gz,
        bed_file = vcfdist_bed_file,
        region = region,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        docker = vcfdist_docker,
        verbosity = vcfdist_verbosity
    }

    output{
    }
}


