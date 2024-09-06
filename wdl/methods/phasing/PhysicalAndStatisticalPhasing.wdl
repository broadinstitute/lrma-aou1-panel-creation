version 1.0

import "./Helper.wdl" as H


workflow PhysicalAndStatisticalPhasing {

    input {
        Array[File] sample_bams
        Array[File] sample_bais
        File joint_short_vcf
        File joint_short_vcf_tbi
        File joint_sv_vcf
        File joint_sv_vcf_tbi
        File reference_fasta
        File reference_fasta_fai
        File genetic_mapping_tsv_for_shapeit4
        String chromosome
        String region
        String prefix
        String gcs_out_root_dir
        Int shapeit4_num_threads
        Int merge_num_threads = 4
        Int hiphase_memory
        Int shapeit4_memory
        String shapeit4_extra_args = "--use-PS 0.0001" # expected error rate in phase sets derived from physical phasing
        String hiphase_extra_args
    }

    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)
    Int data_length = length(sample_bams)
    Array[Int] indexes= range(data_length)

    # call H.SubsetVCF as SubsetVcfShort { input:
    #     vcf_gz = joint_short_vcf,
    #     vcf_tbi = joint_short_vcf_tbi,
    #     locus = region
    # }

    # call H.SubsetVCF as SubsetVcfSV { input:
    #     vcf_gz = joint_sv_vcf,
    #     vcf_tbi = joint_sv_vcf_tbi,
    #     locus = region
    # }

    call UnphaseGenotypes as UnphaseSVGenotypes { input:
        vcf = joint_sv_vcf,
        vcf_tbi = joint_sv_vcf_tbi,
        prefix = prefix + ".unphased"
    }

    scatter (idx in indexes)  {
        File all_chr_bam = sample_bams[idx]
        File all_chr_bai = sample_bais[idx]

        # call H.SubsetBam { input:
        #     bam = all_chr_bam,
        #     bai = all_chr_bai,
        #     locus = region
        # }

        call H.InferSampleName { input: 
            bam = all_chr_bam, 
            bai = all_chr_bai
        }

        String sample_id = InferSampleName.sample_name

        # call H.SplitVCFbySample as SplitVcfbySampleShort { input:
        #     joint_vcf = joint_short_vcf,
        #     region = region,
        #     samplename = sample_id
        # }

        # call H.SplitVCFbySample as SplitVcfbySampleSV { input:
        #     joint_vcf = UnphaseSVGenotypes.unphased_vcf,
        #     region = region,
        #     samplename = sample_id
        # }

        call ConvertLowerCase {
            input:
                vcf = UnphaseSVGenotypes.unphased_vcf,
                prefix = sample_id + ".uppercased_sv_cleaned"
                
        }

        call H.HiPhase { input:
            bam = all_chr_bam,
            bai = all_chr_bai,
            unphased_snp_vcf = joint_short_vcf,
            unphased_snp_tbi = joint_short_vcf_tbi,
            unphased_sv_vcf = ConvertLowerCase.subset_vcf,
            unphased_sv_tbi = ConvertLowerCase.subset_tbi,
            ref_fasta = reference_fasta,
            ref_fasta_fai = reference_fasta_fai,
            samplename = sample_id,
            memory = hiphase_memory,
            extra_args = hiphase_extra_args
        }
    }

    call H.MergePerChrVcfWithBcftools as MergeAcrossSamplesShort { input:
        vcf_input = HiPhase.phased_snp_vcf,
        tbi_input = HiPhase.phased_snp_vcf_tbi,
        pref = prefix + ".short",
        threads_num = merge_num_threads
    }

    call H.MergePerChrVcfWithBcftools as MergeAcrossSamplesSV { input:
        vcf_input = HiPhase.phased_sv_vcf,
        tbi_input = HiPhase.phased_sv_vcf_tbi,
        pref = prefix + ".SV",
        threads_num = merge_num_threads
    }

    call FilterAndConcatVcfs { input:
        short_vcf = MergeAcrossSamplesShort.merged_vcf,
        short_vcf_tbi = MergeAcrossSamplesShort.merged_tbi,
        sv_vcf = MergeAcrossSamplesSV.merged_vcf,
        sv_vcf_tbi = MergeAcrossSamplesSV.merged_tbi,
        prefix = prefix + ".filter_and_concat"
    }

    call H.Shapeit4 { input:
        vcf_input = FilterAndConcatVcfs.filter_and_concat_vcf,
        vcf_index = FilterAndConcatVcfs.filter_and_concat_vcf_tbi,
        mappingfile = genetic_mapping_dict[chromosome],
        region = region,
        prefix = prefix + ".filter_and_concat.phased",
        num_threads = shapeit4_num_threads,
        memory = shapeit4_memory,
        extra_args = shapeit4_extra_args
    }

    output {
        File hiphase_short_vcf = MergeAcrossSamplesShort.merged_vcf
        File hiphase_short_tbi = MergeAcrossSamplesShort.merged_tbi
        File hiphase_sv_vcf = MergeAcrossSamplesSV.merged_vcf
        File hiphase_sv_tbi = MergeAcrossSamplesSV.merged_tbi
        File filtered_vcf = FilterAndConcatVcfs.filter_and_concat_vcf
        File filtered_tbi = FilterAndConcatVcfs.filter_and_concat_vcf_tbi
        File phased_bcf = Shapeit4.phased_bcf
    }
}

task ConvertLowerCase {
    input {
        File vcf
        String prefix
    }

    Int disk_size = 2*ceil(size([vcf], "GB")) + 1
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cp ~{docker_dir}/convert_lower_case.py ~{work_dir}/convert_lower_case.py
        cd ~{work_dir}

        python convert_lower_case.py -i ~{vcf} -o ~{prefix}.vcf
        bgzip ~{prefix}.vcf ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{work_dir}/~{prefix}.vcf.gz"
        File subset_tbi = "~{work_dir}/~{prefix}.vcf.gz.tbi"
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

# filter out singletons (i.e., keep MAC >= 2) and concatenate with deduplication
task FilterAndConcatVcfs {

    input {
        File short_vcf         # multiallelic
        File short_vcf_tbi
        File sv_vcf            # biallelic
        File sv_vcf_tbi
        String prefix
    }

    command <<<
        set -euxo pipefail

        # filter SV singletons
        bcftools view -i 'MAC>=2' ~{sv_vcf} \
            --write-index -Oz -o ~{prefix}.SV.vcf.gz

        # filter short singletons and split to biallelic
        bcftools view -i 'MAC>=2' ~{short_vcf} | \
            bcftools norm -m-any --do-not-normalize \
            --write-index -Oz -o ~{prefix}.short.vcf.gz

        # concatenate with deduplication; providing SV VCF as first argument preferentially keeps those records
        bcftools concat \
            ~{prefix}.SV.vcf.gz \
            ~{prefix}.short.vcf.gz \
            --allow-overlaps --remove-duplicates \
            -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File filter_and_concat_vcf = "~{prefix}.vcf.gz"
        File filter_and_concat_vcf_tbi = "~{prefix}.vcf.gz.tbi"
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

task UnphaseGenotypes {

    input {
        File vcf
        File vcf_tbi
        String prefix
    }

    command <<<
        set -euxo pipefail

        # set (a)ll genotypes to (u)nphased and sort by allele (e.g., 1|0 becomes 0/1)
        bcftools +setGT ~{vcf} -Oz -o ~{prefix}.vcf.gz -- --target-gt a --new-gt u
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File unphased_vcf = "~{prefix}.vcf.gz"
        File unphased_vcf_tbi = "~{prefix}.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}
