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
        Int shapeit4_num_threads
        Int merge_num_threads
        Int hiphase_memory
        Int shapeit4_memory
        String shapeit4_extra_args
        String hiphase_extra_args

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
        shapeit4_num_threads = shapeit4_num_threads,
        merge_num_threads = merge_num_threads,
        hiphase_memory = hiphase_memory,
        shapeit4_memory = shapeit4_memory,
        shapeit4_extra_args = shapeit4_extra_args,
        hiphase_extra_args = hiphase_extra_args
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

    call SummarizeEvaluations { input:
        labels_per_vcf = ["HiPhaseShort", "HiPhaseSV", "ConcatAndFiltered", "Shapeit4", "Panel"],
        vcfdist_outputs_per_vcf_and_sample = [
            EvaluateHiPhaseShort.vcfdist_outputs_per_sample,
            EvaluateHiPhaseSV.vcfdist_outputs_per_sample,
            EvaluateFiltered.vcfdist_outputs_per_sample,
            EvaluateShapeit4.vcfdist_outputs_per_sample,
            EvaluatePanel.vcfdist_outputs_per_sample
        ],
        overlap_metrics_outputs_per_vcf = [
            EvaluateHiPhaseShort.overlap_metrics_outputs,
            EvaluateHiPhaseSV.overlap_metrics_outputs,
            EvaluateFiltered.overlap_metrics_outputs,
            EvaluateShapeit4.overlap_metrics_outputs,
            EvaluatePanel.overlap_metrics_outputs
        ]
    }

    output {
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
