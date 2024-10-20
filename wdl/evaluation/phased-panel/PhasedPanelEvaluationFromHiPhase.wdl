version 1.0

import "../../methods/phasing/PhysicalAndStatisticalPhasing.wdl"
import "../../methods/pangenie/PanGeniePanelCreation.wdl"
import "./VcfdistAndOverlapMetricsEvaluation.wdl"
import "../kage/LeaveOutEvaluation.wdl"
import "../../methods/phasing/Helper.wdl"

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
        String? extra_filter_args
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
        Array[File] genetic_maps
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

    if (subset_short_to_sv_windows) {
        call SubsetVcfShortInSVWindows { input:
            short_vcf_gz = hiphase_short_vcf_gz,
            short_vcf_tbi = hiphase_short_vcf_gz_tbi,
            sv_vcf_gz = hiphase_sv_vcf_gz,
            sv_vcf_tbi = hiphase_sv_vcf_gz_tbi,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            reference_fasta_dict = reference_fasta_dict,
            prefix = output_prefix + ".short.subset.windowed",
            region = chromosomes_regions_arg,
            window_padding = window_padding,
            filter_args = subset_filter_args
        }
    }

    call FilterAndConcatVcfs { input:
        short_vcf = select_first([SubsetVcfShortInSVWindows.subset_short_vcf_gz, hiphase_short_vcf_gz]),
        short_vcf_tbi = select_first([SubsetVcfShortInSVWindows.subset_short_vcf_gz_tbi, hiphase_short_vcf_gz_tbi]),
        sv_vcf = hiphase_sv_vcf_gz,
        sv_vcf_tbi = hiphase_sv_vcf_gz_tbi,
        region = chromosomes_regions_arg,
        extra_filter_args = extra_filter_args,
        prefix = output_prefix + ".filter_and_concat"
    }

    call FixVariantCollisions as BeforeShapeit4FixVariantCollisions { input:
        phased_bcf = FilterAndConcatVcfs.filter_and_concat_vcf,
        fix_variant_collisions_java = fix_variant_collisions_java,
        operation = operation,
        weight_tag = weight_tag,
        is_weight_format_field = is_weight_format_field,
        output_prefix = output_prefix
    }

    scatter (i in range(length(chromosomes))) {
        String chromosome = chromosomes[i]

        call CreateShapeit4Chunks { input:
            vcf = BeforeShapeit4FixVariantCollisions.phased_collisionless_bcf,
            tbi = BeforeShapeit4FixVariantCollisions.phased_collisionless_bcf_csi,
            region = chromosomes[i],
            prefix = output_prefix + "." + chromosome,
            extra_chunk_args = extra_chunk_args
        }

        Array[String] region_list = read_lines(CreateShapeit4Chunks.chunks)

        scatter (j in range(length(region_list))) {
            call Helper.Shapeit4 as Shapeit4 { input:
                vcf_input = BeforeShapeit4FixVariantCollisions.phased_collisionless_bcf,
                vcf_index = BeforeShapeit4FixVariantCollisions.phased_collisionless_bcf_csi,
                mappingfile = genetic_mapping_dict[chromosome],
                region = region_list[j],
                prefix = output_prefix + "." + chromosome + ".shard-" + j + ".filter_and_concat.phased",
                num_threads = shapeit4_num_threads,
                memory = shapeit4_memory,
                extra_args = shapeit4_extra_args
            }
        }

        call LigateVcfs as LigateVcfsChromosome { input:
            vcfs = Shapeit4.phased_bcf,
            prefix = output_prefix + "." + chromosome + ".phased.ligated"
        }
    }

    call LigateVcfs { input:
        vcfs = LigateVcfsChromosome.ligated_vcf,
        prefix = output_prefix + ".phased.ligated"
    }

    call FixVariantCollisions { input:
        phased_bcf = LigateVcfs.ligated_vcf,
        fix_variant_collisions_java = fix_variant_collisions_java,
        operation = operation,
        weight_tag = weight_tag,
        is_weight_format_field = is_weight_format_field,
        output_prefix = output_prefix
    }

    call PanGeniePanelCreation.PanGeniePanelCreation { input:
        phased_bcf = FixVariantCollisions.phased_collisionless_bcf,
        reference_fasta = reference_fasta,
        prepare_vcf_script = prepare_vcf_script,
        add_ids_script = add_ids_script,
        merge_vcfs_script = merge_vcfs_script,
        frac_missing = frac_missing,
        output_prefix = output_prefix,
        docker = panel_creation_docker,
        monitoring_script = monitoring_script
    }

    call LeaveOutEvaluation.LeaveOutEvaluation { input:
        input_vcf_gz = PanGeniePanelCreation.panel_vcf_gz,
        input_vcf_gz_tbi = PanGeniePanelCreation.panel_vcf_gz_tbi,
        case_reference_fasta = case_reference_fasta,
        case_reference_fasta_fai = case_reference_fasta_fai,
        case_reference_dict = case_reference_dict,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        genetic_maps = genetic_maps,
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
    call Helper.MergePerChrVcfWithBcftools as GLIMPSEMergeAcrossSamples { input:
        vcf_input = LeaveOutEvaluation.glimpse_vcf_gzs,
        tbi_input = LeaveOutEvaluation.glimpse_vcf_gz_tbis,
        pref = output_prefix + ".glimpse.merged",
        threads_num = merge_num_threads
    }

    call FixVariantCollisions as GenotypingFixVariantCollisions { input:
        phased_bcf = GLIMPSEMergeAcrossSamples.merged_vcf,
        fix_variant_collisions_java = fix_variant_collisions_java,
        operation = operation,
        weight_tag = weight_tag,
        is_weight_format_field = is_weight_format_field,
        output_prefix = output_prefix
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
        eval_vcf = FilterAndConcatVcfs.filter_and_concat_vcf,
        eval_vcf_idx = FilterAndConcatVcfs.filter_and_concat_vcf_tbi,
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
        eval_vcf = BeforeShapeit4FixVariantCollisions.phased_collisionless_bcf,
        eval_vcf_idx = BeforeShapeit4FixVariantCollisions.phased_collisionless_bcf_csi,
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
        eval_vcf = LigateVcfs.ligated_vcf,
        eval_vcf_idx = LigateVcfs.ligated_vcf_tbi,
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
        eval_vcf = FixVariantCollisions.phased_collisionless_bcf,
        eval_vcf_idx = FixVariantCollisions.phased_collisionless_bcf_csi,
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
        eval_vcf = PanGeniePanelCreation.panel_vcf_gz,
        eval_vcf_idx = PanGeniePanelCreation.panel_vcf_gz_tbi,
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
        eval_vcf = GLIMPSEMergeAcrossSamples.merged_vcf,
        eval_vcf_idx = GLIMPSEMergeAcrossSamples.merged_tbi,
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
        eval_vcf = GenotypingFixVariantCollisions.phased_collisionless_bcf,
        eval_vcf_idx = GenotypingFixVariantCollisions.phased_collisionless_bcf_csi,
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
        call Helper.MergePerChrVcfWithBcftools as PanGenieMergeAcrossSamples { input:
            vcf_input = select_all(LeaveOutEvaluation.pangenie_vcf_gzs),
            tbi_input = select_all(LeaveOutEvaluation.pangenie_vcf_gz_tbis),
            pref = output_prefix + ".pangenie.merged",
            threads_num = merge_num_threads
        }

        # summarize PanGenie metrics vs. panel

        # evaluate PanGenie
        call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluatePanGenie { input:
            samples = vcfdist_samples,
            truth_vcf = vcfdist_truth_vcf,
            truth_vcf_idx = vcfdist_truth_vcf_idx,
            eval_vcf = PanGenieMergeAcrossSamples.merged_vcf,
            eval_vcf_idx = PanGenieMergeAcrossSamples.merged_tbi,
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
        ])
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

task SubsetVcfShortInSVWindows {
    input {
        File short_vcf_gz
        File short_vcf_tbi
        File sv_vcf_gz
        File sv_vcf_tbi
        File reference_fasta
        File reference_fasta_fai
        File reference_fasta_dict
        String prefix
        String region
        Int window_padding
        String? filter_args

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([short_vcf_gz, short_vcf_tbi], "GB")) + 2*ceil(size([sv_vcf_gz, sv_vcf_tbi], "GB")) + 1

    command {
        set -euxo pipefail

        bcftools view ~{sv_vcf_gz} \
            -r ~{region} \
            -Oz -o sv.region.vcf.gz
        bcftools index -t sv.region.vcf.gz
        gatk PreprocessIntervals \
            -L sv.region.vcf.gz \
            --reference ~{reference_fasta} \
            --padding ~{window_padding} \
            --bin-length 0 \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{prefix}.windows.interval_list
        gatk IntervalListToBed \
            -I ~{prefix}.windows.interval_list \
            -O ~{prefix}.windows.bed
        bcftools view ~{short_vcf_gz} \
            -R ~{prefix}.windows.bed \
            ~{filter_args} \
            -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    }

    output {
        File subset_short_vcf_gz = "~{prefix}.vcf.gz"
        File subset_short_vcf_gz_tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:"us.gcr.io/broad-gatk/gatk:4.6.0.0"
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

# filter out singletons (i.e., keep MAC >= 2) and concatenate with deduplication
task FilterAndConcatVcfs {

    input {
        File short_vcf         # multiallelic
        File short_vcf_tbi
        File sv_vcf            # biallelic
        File sv_vcf_tbi
        String prefix
        String region
        String? extra_filter_args = "-i 'MAC>=2'"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([short_vcf, short_vcf_tbi], "GB")) + 2*ceil(size([sv_vcf, sv_vcf_tbi], "GB")) + 1

    command {
        set -euxo pipefail

        # filter SV singletons
        bcftools view ~{extra_filter_args} ~{sv_vcf} \
            -r ~{region} \
            --write-index -Oz -o ~{prefix}.SV.vcf.gz

        # filter short singletons and split to biallelic
        bcftools view ~{extra_filter_args} ~{short_vcf} \
            -r ~{region} | \
            bcftools norm -m-any --do-not-normalize \
            --write-index -Oz -o ~{prefix}.short.vcf.gz

        # concatenate with deduplication; providing SV VCF as first argument preferentially keeps those records
        bcftools concat \
            ~{prefix}.SV.vcf.gz \
            ~{prefix}.short.vcf.gz \
            --allow-overlaps --remove-duplicates \
            -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    }

    output {
        File filter_and_concat_vcf = "~{prefix}.vcf.gz"
        File filter_and_concat_vcf_tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
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

task CreateShapeit4Chunks {

    input {
        File vcf
        File tbi
        String region
        String prefix
        String? extra_chunk_args = "--thread $(nproc) --window-size 5000000 --buffer-size 500000"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf, tbi], "GB")) + 1

    command <<<
        set -euxo pipefail

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_chunk_static
        chmod +x GLIMPSE_chunk_static

        ./GLIMPSE_chunk_static \
            -I ~{vcf} \
            --region ~{region} \
            ~{extra_chunk_args} \
            -O chunks.txt

        # cut chunks + buffers
        cut -f 3 chunks.txt > chunks.regions.txt
    >>>

    output {
        File chunks = "chunks.regions.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
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

task LigateVcfs {

    input {
        Array[File] vcfs
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(vcfs, "GB")) + 1

    command <<<
        set -euxo pipefail
        for ff in ~{sep=' ' vcfs}; do bcftools index $ff; done
        bcftools concat --ligate  ~{sep=" " vcfs} -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File ligated_vcf = "~{prefix}.vcf.gz"
        File ligated_vcf_tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
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
            bcftools view -Oz -o ~{output_prefix}.phased.collisionless.vcf.gz
        # index and convert via vcf.gz to avoid errors from missing header lines
        bcftools index -t ~{output_prefix}.phased.collisionless.vcf.gz
        bcftools view ~{output_prefix}.phased.collisionless.vcf.gz -Ob -o ~{output_prefix}.phased.collisionless.bcf
        bcftools index ~{output_prefix}.phased.collisionless.bcf
    >>>

    output {
        File phased_collisionless_bcf = "~{output_prefix}.phased.collisionless.bcf"
        File phased_collisionless_bcf_csi = "~{output_prefix}.phased.collisionless.bcf.csi"
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
