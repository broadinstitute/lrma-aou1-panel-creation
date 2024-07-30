version 1.0

import "./Helper.wdl" as H


workflow PhysicalAndStatisticalPhasing {
    parameter_meta {
        genetic_mapping_tsv_for_shapeit4: "path to the tsv file for the genetic mapping file address per chromosome"
        chromosome: "string for chromosome to be processed"
        prefix: "output file prefix, usually using chromosome"
        num_t: "integer for threads"
    }

    input {
        Array[File] bams_from_all_samples
        Array[File] bais_from_all_samples
        File joint_vcf
        File joint_vcf_tbi
        File joint_sv
        File joint_sv_tbi
        File reference
        File reference_index
        File genetic_mapping_tsv_for_shapeit4
        String chromosome
        String region
        String prefix
        String gcs_out_root_dir
        Int num_t
        Int merge_num_t = 4
        Int hiphase_memory
        Int shapeit4_memory
    }

    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)
    Int data_length = length(bams_from_all_samples)
    Array[Int] indexes= range(data_length)

    call H.SubsetVCF as SubsetSNPsJoint { input:
        vcf_gz = joint_vcf,
        vcf_tbi = joint_vcf_tbi,
        locus = region
    }

    call H.SubsetVCF as SubsetSVsJoint { input:
        vcf_gz = joint_sv,
        vcf_tbi = joint_sv_tbi,
        locus = region
    }

    scatter (idx in indexes)  {
        File all_chr_bam = bams_from_all_samples[idx]
        File all_chr_bai = bais_from_all_samples[idx]

        call H.SubsetBam as SubsetBam { input:
            bam = all_chr_bam,
            bai = all_chr_bai,
            locus = region
        }

        call H.InferSampleName { input: 
            bam = all_chr_bam, 
            bai = all_chr_bai
        }

        String sample_id = InferSampleName.sample_name

        call H.SplitVCFbySample as SP { input:
            joint_vcf = SubsetSNPsJoint.subset_vcf,
            region = region,
            samplename = sample_id
        }

        call H.SplitVCFbySample as SV_split { input:
            joint_vcf = SubsetSVsJoint.subset_vcf,
            region = region,
            samplename = sample_id
        }

        call convert as cleanvcf{
            input:
                vcf = SV_split.single_sample_vcf,
                samplename = sample_id
                
        }

        call H.HiphaseSVs as HP_SV { input:
            bam = SubsetBam.subset_bam,
            bai = SubsetBam.subset_bai,
            unphased_snp_vcf = SP.single_sample_vcf,
            unphased_snp_tbi = SP.single_sample_vcf_tbi,
            unphased_sv_vcf = cleanvcf.subset_vcf,
            unphased_sv_tbi = cleanvcf.subset_tbi,
            ref_fasta = reference,
            ref_fasta_fai = reference_index,
            samplename = sample_id,
            memory = hiphase_memory
        }
    }

    call H.MergePerChrVcfWithBcftools as MergeAcrossSamples { input:
        vcf_input = HP_SV.phased_snp_vcf,
        tbi_input = HP_SV.phased_snp_vcf_tbi,
        pref = prefix + "short",
        threads_num = merge_num_t
    }

    call SplitMultiallelicCalls as SplitMultiallelicCallsShort { input:
        bcftools_merged_vcf = MergeAcrossSamples.merged_vcf,
        bcftools_merged_vcf_tbi = MergeAcrossSamples.merged_tbi,
        prefix = prefix + "short.split"
    }

    call H.MergePerChrVcfWithBcftools as MergeAcrossSamplesSV { input:
        vcf_input = HP_SV.phased_sv_vcf,
        tbi_input = HP_SV.phased_sv_vcf_tbi,
        pref = prefix + "SV",
        threads_num = merge_num_t
    }

    call ConcatVCFs { input:
        bcftools_small_vcf = SplitMultiallelicCallsShort.normed_vcf,
        bcftools_small_vcf_tbi = SplitMultiallelicCallsShort.normed_vcf_tbi,
        bcftools_sv_vcf = MergeAcrossSamplesSV.merged_vcf,
        bcftools_sv_vcf_tbi = MergeAcrossSamplesSV.merged_tbi,
        prefix = prefix + "concat_short_and_SV"
    }

    call FilterVariants as FilterShortAndSV { input:
        bcftools_vcf = ConcatVCFs.concat_vcf,
        bcftools_vcf_tbi = ConcatVCFs.concat_vcf_tbi,
        prefix = prefix + "concat_short_and_SV.filter"
    }

    call H.Shapeit4 as Shapeit4ShortAndSV { input:
        vcf_input = FilterShortAndSV.filtered_vcf,
        vcf_index = FilterShortAndSV.filtered_vcf_tbi,
        mappingfile = genetic_mapping_dict[chromosome],
        region = region,
        prefix = prefix + "concat_short_and_SV.filter.phased",
        num_threads = num_t,
        memory = shapeit4_memory
    }

    output {
        File hiphase_short_vcf = MergeAcrossSamples.merged_vcf
        File hiphase_short_tbi = MergeAcrossSamples.merged_tbi
        File hiphase_sv_vcf = MergeAcrossSamplesSV.merged_vcf
        File hiphase_sv_tbi = MergeAcrossSamplesSV.merged_tbi
        File filtered_vcf = FilterShortAndSV.filtered_vcf
        File filtered_tbi = FilterShortAndSV.filtered_vcf_tbi
        File phased_bcf = Shapeit4ShortAndSV.phased_bcf
    }
}

task convert {
    meta {
        description: ""
    }
    parameter_meta {
        vcf: "VCF file to be cleaned"
    }
    input {
        File vcf
        String samplename
    }

    Int disk_size = 2*ceil(size([vcf], "GB")) + 1
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cp ~{docker_dir}/convert_lower_case.py ~{work_dir}/convert_lower_case.py
        cd ~{work_dir}

        python convert_lower_case.py -i ~{vcf} -o ~{samplename}_sv_cleaned.vcf
        bgzip ~{samplename}_sv_cleaned.vcf ~{samplename}_sv_cleaned.vcf.gz
        tabix -p vcf ~{samplename}_sv_cleaned.vcf.gz
    >>>

    output {
        File subset_vcf = "~{work_dir}/~{samplename}_sv_cleaned.vcf.gz"
        File subset_tbi = "~{work_dir}/~{samplename}_sv_cleaned.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 2
        memory:  "32 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"hangsuunc/cleanvcf:v1"
    }

}

task SplitMultiallelicCalls {

    meta {
        description: "using bcftools to splic multiallelic calls"
    }

    parameter_meta {
        bcftools_merged_vcf: "bcftools merged vcf files"
    }

    input {
        File bcftools_merged_vcf
        File bcftools_merged_vcf_tbi
        String prefix
    }

    command <<<
        set -euxo pipefail
        bcftools norm -m -any ~{bcftools_merged_vcf} | bgzip > ~{prefix}.normed.vcf.gz
        tabix -p vcf ~{prefix}.normed.vcf.gz
        
    >>>

    output {
        File normed_vcf = " ~{prefix}.normed.vcf.gz"
        File normed_vcf_tbi = "~{prefix}.normed.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }

}

task FilterVariants {

    meta {
        description: "using bcftools to filter variants by F-missing and MAC"
    }

    parameter_meta {
        bcftools_vcf: "bcftools merged vcf files"
    }

    input {
        File bcftools_vcf
        File bcftools_vcf_tbi
        String prefix
    }

    command <<<
        set -euxo pipefail
        bcftools view -i 'MAC>=2' ~{bcftools_vcf} | bgzip > ~{prefix}.filtered.vcf.gz # F_MISSING < 0.05 & 
        tabix -p vcf ~{prefix}.filtered.vcf.gz
        
    >>>

    output {
        File filtered_vcf = "~{prefix}.filtered.vcf.gz"
        File filtered_vcf_tbi = "~{prefix}.filtered.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }

}

task ConcatVCFs {

    meta {
        description: "using bcftools to concatenate small variants and SVs"
    }

    parameter_meta {
        bcftools_small_vcf: "bcftools merged vcf files"
    }

    input {
        File bcftools_small_vcf
        File bcftools_small_vcf_tbi
        File bcftools_sv_vcf
        File bcftools_sv_vcf_tbi
        String prefix
    }

    command <<<
        set -euxo pipefail
        bcftools concat ~{bcftools_small_vcf} ~{bcftools_sv_vcf} -Oz -o ~{prefix}_integrated.vcf.gz
        bcftools sort ~{prefix}_integrated.vcf.gz -O z -o ~{prefix}_integrated_sorted.vcf.gz
        tabix -p vcf ~{prefix}_integrated_sorted.vcf.gz
        
    >>>

    output {
        File concat_vcf = "~{prefix}_integrated_sorted.vcf.gz"
        File concat_vcf_tbi = "~{prefix}_integrated_sorted.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }

}
