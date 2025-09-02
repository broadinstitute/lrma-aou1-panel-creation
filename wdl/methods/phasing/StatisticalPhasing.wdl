version 1.0

import "./Helper.wdl" as H


workflow StatisticalPhasing {

    input {

        File joint_short_vcf
        File joint_short_vcf_tbi
        File joint_sv_vcf
        File joint_sv_vcf_tbi
        File reference_fasta
        File reference_fasta_fai
        File genetic_mapping_tsv_for_shapeit
        String chromosome
        String region
        String prefix
        String gcs_out_root_dir
        Int shapeit_num_threads
        Int merge_num_threads = 4

        Int shapeit_memory
        String shapeit_extra_args 
    }

    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit)


    call H.SubsetVCF as SubsetVcfShort { input:
        vcf_gz = joint_short_vcf,
        vcf_tbi = joint_short_vcf_tbi,
        locus = region
    }

    call H.SubsetVCF as SubsetVcfSV { input:
        vcf_gz = joint_sv_vcf,
        vcf_tbi = joint_sv_vcf_tbi,
        locus = region
    }

    call FilterAndConcatVcfs { input:
        short_vcf = SubsetVcfShort.subset_vcf,
        short_vcf_tbi = SubsetVcfShort.subset_tbi,
        sv_vcf = SubsetVcfSV.subset_vcf,
        sv_vcf_tbi = SubsetVcfSV.subset_tbi,
        prefix = prefix + ".concat"
    }

    call H.shapeit5_phase_common as Shapeit5_phase_common { input:
        vcf_input = FilterAndConcatVcfs.filter_and_concat_vcf,
        vcf_index = FilterAndConcatVcfs.filter_and_concat_vcf_tbi,
        mappingfile = genetic_mapping_dict[chromosome],
        region = region,
        prefix = prefix + ".filter_and_concat.phased",
        num_threads = shapeit_num_threads,
        memory = shapeit_memory,
        extra_args = shapeit_extra_args
    }

    output {
        File filtered_vcf = FilterAndConcatVcfs.filter_and_concat_vcf
        File filtered_tbi = FilterAndConcatVcfs.filter_and_concat_vcf_tbi
        File phased_common_bcf = Shapeit5_phase_common.scaffold_vcf
        File phased_common_tbi = Shapeit5_phase_common.scaffold_vcf_index
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
