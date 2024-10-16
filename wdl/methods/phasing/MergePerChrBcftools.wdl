version 1.0

import "./Helper.wdl" as H


workflow MergePerchrBcftools {

    input {
        Array[File] hiphased_snp_vcf
        Array[File] hiphased_snp_vcf_tbi
        Array[File] hiphased_sv_vcf
        Array[File] hiphased_sv_vcf_tbi

        String prefix
        String gcs_out_root_dir
        Int merge_num_threads = 4
        Int batchsize

    }


    call H.MergePerChrVcfWithBcftools as MergeAcrossSamplesShort { input:
        vcf_input = hiphased_snp_vcf,
        tbi_input = hiphased_snp_vcf_tbi,
        pref = prefix + ".short",
        threads_num = merge_num_threads,
        batch_size = batchsize
    }

    call H.MergePerChrVcfWithBcftools as MergeAcrossSamplesSV { input:
        vcf_input = hiphased_sv_vcf,
        tbi_input = hiphased_sv_vcf_tbi,
        pref = prefix + ".SV",
        threads_num = merge_num_threads,
        batch_size = batchsize
    }

    call FilterAndConcatVcfs { input:
        short_vcf = MergeAcrossSamplesShort.merged_vcf,
        short_vcf_tbi = MergeAcrossSamplesShort.merged_tbi,
        sv_vcf = MergeAcrossSamplesSV.merged_vcf,
        sv_vcf_tbi = MergeAcrossSamplesSV.merged_tbi,
        prefix = prefix + ".filter_and_concat"
    }
    output {
        File filtered_vcf = FilterAndConcatVcfs.filter_and_concat_vcf
        File filtered_tbi = FilterAndConcatVcfs.filter_and_concat_vcf_tbi
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
        cpu: 4
        memory:  "16 GiB"
        disks: "local-disk 375 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }
}

