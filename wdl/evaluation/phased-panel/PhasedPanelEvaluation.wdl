version 1.0

import "../../methods/phasing/PhysicalAndStatisticalPhasing.wdl"
import "../vcfdist/VcfdistEvaluation.wdl"

workflow PhasedPanelEvaluation {

    input {
        # common inputs
        File reference_fasta
        File reference_fai
        String region
        String gcs_out_root_dir

        # inputs for phasing workflow
        Array[File] bams_from_all_samples
        Array[File] bais_from_all_samples
        File joint_vcf
        File joint_vcf_tbi
        File joint_sv
        File joint_sv_tbi
        File genetic_mapping_tsv_for_shapeit4
        String chromosome
        String prefix
        Int num_t

        # inputs for vcfdist workflow
        Array[String] evaluation_samples
        File truth_vcf
        File bed_file
        String docker = "timd1/vcfdist:v2.5.3"
        Int verbosity = 1
    }

    call PhysicalAndStatisticalPhasing.PhysicalAndStatisticalPhasing as PhasingWorkflow { input:
        bams_from_all_samples = bams_from_all_samples,
        bais_from_all_samples = bais_from_all_samples,
        joint_vcf = joint_vcf,
        joint_vcf_tbi = joint_vcf_tbi,
        joint_sv = joint_sv,
        joint_sv_tbi = joint_sv_tbi,
        reference = reference_fasta,
        reference_index = reference_fai,
        genetic_mapping_tsv_for_shapeit4 = genetic_mapping_tsv_for_shapeit4,
        chromosome = chromosome,
        region = region,
        prefix = prefix,
        gcs_out_root_dir = gcs_out_root_dir + "/Phasing",
        num_t = num_t
    }

    # evaluate HiPhase short
    call VcfdistEvaluation.VcfdistEvaluation as VcfdistEvaluationHiPhaseShort { input:
        samples = evaluation_samples,
        truth_vcf = truth_vcf,
        eval_vcf = PhasingWorkflow.hiphase_short_vcf,
        bed_file = bed_file,
        region = region,
        reference_fasta = reference_fasta,
        reference_fai = reference_fai,
        docker = docker,
        verbosity = verbosity
    }

    # evaluate HiPhase SV
    call VcfdistEvaluation.VcfdistEvaluation as VcfdistEvaluationHiPhaseSV { input:
        samples = evaluation_samples,
        truth_vcf = truth_vcf,
        eval_vcf = PhasingWorkflow.hiphase_sv_vcf,
        bed_file = bed_file,
        region = region,
        reference_fasta = reference_fasta,
        reference_fai = reference_fai,
        docker = docker,
        verbosity = verbosity
    }

    # evaluate filtered short + SV
    call VcfdistEvaluation.VcfdistEvaluation as VcfdistEvaluationFiltered { input:
        samples = evaluation_samples,
        truth_vcf = truth_vcf,
        eval_vcf = PhasingWorkflow.filtered_vcf,
        bed_file = bed_file,
        region = region,
        reference_fasta = reference_fasta,
        reference_fai = reference_fai,
        docker = docker,
        verbosity = verbosity
    }

    # evaluate Shapeit4 short + SV
    call VcfdistEvaluation.VcfdistEvaluation as VcfdistEvaluationShapeit4 { input:
        samples = evaluation_samples,
        truth_vcf = truth_vcf,
        eval_vcf = PhasingWorkflow.phased_bcf,
        bed_file = bed_file,
        region = region,
        reference_fasta = reference_fasta,
        reference_fai = reference_fai,
        docker = docker,
        verbosity = verbosity
    }

    output{
    }
}


