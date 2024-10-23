version 1.0

import "./ChromosomePhasedPanelCreationFromHiPhase.wdl"
import "./ConcatAndEvaluate.wdl"

workflow PhasedPanelEvaluation {

    input {
        # common inputs
        File reference_fasta
        File reference_fasta_fai
        File reference_fasta_dict
        Array[String] chromosomes
        String gcs_out_root_dir
        String output_prefix
        Int merge_num_threads
        File? monitoring_script

        # inputs for PhysicalAndStatisticalPhasing
        File hiphase_short_vcf_gz
        File hiphase_short_vcf_gz_tbi
        File hiphase_sv_vcf_gz
        File hiphase_sv_vcf_gz_tbi
        Boolean subset_short_to_sv_windows
        Int window_padding
        String? subset_filter_args
        String? filter_and_concat_short_filter_args
        String? filter_and_concat_sv_filter_args
        String? extra_chunk_args
        File genetic_mapping_tsv_for_shapeit4
        Int shapeit4_num_threads
        Int shapeit4_memory
        String shapeit4_extra_args

        # inputs for FixVariantCollisions
        File fix_variant_collisions_java
        Int operation
        String weight_tag
        Int is_weight_format_field

        # inputs for PanGeniePanelCreation
        File prepare_vcf_script
        File add_ids_script
        File merge_vcfs_script
        Float frac_missing
        String panel_creation_docker

        # inputs for VcfdistAndOverlapMetricsEvaluation
        Array[String] vcfdist_samples
        File vcfdist_truth_vcf
        File vcfdist_truth_vcf_idx
        Array[String] evaluation_chromosomes
        File vcfdist_bed_file
        String? vcfdist_extra_args
        String overlap_metrics_docker

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

        String summarize_evaluations_docker
    }

    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)

    call Sep as SepChromosomes { input:
        strs = chromosomes
    }
    String chromosomes_regions_arg = SepChromosomes.str

    call Sep as SepEvaluationChromosomes { input:
        strs = evaluation_chromosomes
    }
    String evaluation_chromosomes_regions_arg = SepEvaluationChromosomes.str

    scatter (i in range(length(chromosomes))) {
        String chromosome = chromosomes[i]

        call ChromosomePhasedPanelCreationFromHiPhase.PhasedPanelEvaluation as ChromosomePhasedPanelCreationFromHiPhase { input:
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            reference_fasta_dict = reference_fasta_dict,
            chromosome = chromosome,
            output_prefix = output_prefix,
            monitoring_script = monitoring_script,
            hiphase_short_vcf_gz = hiphase_short_vcf_gz,
            hiphase_short_vcf_gz_tbi = hiphase_short_vcf_gz_tbi,
            hiphase_sv_vcf_gz = hiphase_sv_vcf_gz,
            hiphase_sv_vcf_gz_tbi = hiphase_sv_vcf_gz_tbi,
            subset_short_to_sv_windows = subset_short_to_sv_windows,
            window_padding = window_padding,
            subset_filter_args = subset_filter_args,
            filter_and_concat_short_filter_args = filter_and_concat_short_filter_args,
            filter_and_concat_sv_filter_args = filter_and_concat_sv_filter_args,
            extra_chunk_args = extra_chunk_args,
            genetic_mapping_tsv_for_shapeit4 = genetic_mapping_tsv_for_shapeit4,
            shapeit4_num_threads = shapeit4_num_threads,
            shapeit4_memory = shapeit4_memory,
            shapeit4_extra_args = shapeit4_extra_args,
            fix_variant_collisions_java = fix_variant_collisions_java,
            operation = operation,
            weight_tag = weight_tag,
            is_weight_format_field = is_weight_format_field,
            prepare_vcf_script = prepare_vcf_script,
            add_ids_script = add_ids_script,
            merge_vcfs_script = merge_vcfs_script,
            frac_missing = frac_missing,
            panel_creation_docker = panel_creation_docker
        }
    }

    call ConcatAndEvaluate.PhasedPanelEvaluation as ConcatAndEvaluate { input:
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        reference_fasta_dict = reference_fasta_dict,
        output_prefix = output_prefix,
        monitoring_script = monitoring_script,
        hiphase_short_vcf_gz = hiphase_short_vcf_gz,
        hiphase_short_vcf_gz_tbi = hiphase_short_vcf_gz_tbi,
        hiphase_sv_vcf_gz = hiphase_sv_vcf_gz,
        hiphase_sv_vcf_gz_tbi = hiphase_sv_vcf_gz_tbi,
        filter_and_concat_vcf_gzs = ChromosomePhasedPanelCreationFromHiPhase.filter_and_concat_vcf_gz,
        filter_and_concat_vcf_gz_tbis = ChromosomePhasedPanelCreationFromHiPhase.filter_and_concat_vcf_gz_tbi,
        before_shapeit4_collisionless_bcfs = ChromosomePhasedPanelCreationFromHiPhase.before_shapeit4_collisionless_bcf,
        before_shapeit4_collisionless_bcf_csis = ChromosomePhasedPanelCreationFromHiPhase.before_shapeit4_collisionless_bcf_csi,
        phased_vcf_gzs = ChromosomePhasedPanelCreationFromHiPhase.phased_vcf_gz,
        phased_vcf_gz_tbis = ChromosomePhasedPanelCreationFromHiPhase.phased_vcf_gz_tbi,
        collisionless_bcfs = ChromosomePhasedPanelCreationFromHiPhase.collisionless_bcf,
        collisionless_bcf_csis = ChromosomePhasedPanelCreationFromHiPhase.collisionless_bcf_csi,
        panel_vcf_gzs = ChromosomePhasedPanelCreationFromHiPhase.panel_vcf_gz,
        panel_vcf_gz_tbis = ChromosomePhasedPanelCreationFromHiPhase.panel_vcf_gz_tbi,
        case_reference_fasta = case_reference_fasta,
        case_reference_fasta_fai = case_reference_fasta_fai,
        case_reference_dict = case_reference_dict,
        leave_out_genetic_maps = leave_out_genetic_maps,
        repeat_mask_bed = repeat_mask_bed,
        segmental_duplications_bed = segmental_duplications_bed,
        simple_repeats_bed = simple_repeats_bed,
        challenging_medically_relevant_genes_bed = challenging_medically_relevant_genes_bed,
        leave_out_chromosomes = leave_out_chromosomes,
        leave_out_sample_names = leave_out_sample_names,
        case_average_coverage = case_average_coverage,
        do_pangenie = do_pangenie,
        leave_out_crams = leave_out_crams,
        leave_out_docker = leave_out_docker,
        kage_docker = kage_docker,
        pangenie_docker = pangenie_docker,
        cpu_make_count_model = cpu_make_count_model,
        leave_out_runtime_attributes = leave_out_runtime_attributes,
        leave_out_medium_runtime_attributes = leave_out_medium_runtime_attributes,
        leave_out_large_runtime_attributes = leave_out_large_runtime_attributes,
        pangenie_runtime_attributes = pangenie_runtime_attributes,
        kage_count_kmers_runtime_attributes = kage_count_kmers_runtime_attributes,
        kage_genotype_runtime_attributes = kage_genotype_runtime_attributes,
        glimpse_case_chromosome_runtime_attributes = glimpse_case_chromosome_runtime_attributes,
        glimpse_case_runtime_attributes = glimpse_case_runtime_attributes,
        calculate_metrics_runtime_attributes = calculate_metrics_runtime_attributes,
        hierarchically_merge_regions = hierarchically_merge_regions,
        hierarchically_merge_batch_size = hierarchically_merge_batch_size,
        hierarchically_merge_docker = hierarchically_merge_docker,
        fix_variant_collisions_java = fix_variant_collisions_java,
        operation = operation,
        weight_tag = weight_tag,
        is_weight_format_field = is_weight_format_field,
        vcfdist_samples = vcfdist_samples,
        vcfdist_truth_vcf = vcfdist_truth_vcf,
        vcfdist_truth_vcf_idx = vcfdist_truth_vcf_idx,
        evaluation_chromosomes_regions_arg = evaluation_chromosomes_regions_arg,
        vcfdist_bed_file = vcfdist_bed_file,
        vcfdist_extra_args = vcfdist_extra_args,
        overlap_metrics_docker = overlap_metrics_docker,
        summarize_evaluations_docker = summarize_evaluations_docker
    }

    output {
    }
}

task Sep {
    input {
        Array[String] strs
    }

    command {
        echo "~{sep=',' strs}" > str.txt
    }

    output {
        String str = read_string("str.txt")
    }

    runtime {
        cpu: 1
        memory: "3 GiB"
        disks: "local-disk 10 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }
}
