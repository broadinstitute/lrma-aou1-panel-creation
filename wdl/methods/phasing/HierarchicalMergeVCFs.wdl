version 1.0

import "./Helper.wdl" as H


workflow HierarchicalMergeVCFs {

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

    output {
        File hiphase_short_vcf = MergeAcrossSamplesShort.merged_vcf
        File hiphase_short_tbi = MergeAcrossSamplesShort.merged_tbi
        File hiphase_sv_vcf = MergeAcrossSamplesSV.merged_vcf
        File hiphase_sv_tbi = MergeAcrossSamplesSV.merged_tbi
        
    }
}


task MergePerChrVcfWithBcftools {
    parameter_meta {
        vcf_input: {localization_optional: true}
        tbi_input: {localization_optional: true}
    }
    input{
        Array[File] vcf_input
        Array[File] tbi_input
        String pref
        Int threads_num
    }

    command <<<
        set -eux

        # we do single-sample phased VCFs localization ourselves
        mkdir -p ssp_vcfs
        time \
        gcloud storage cp ~{sep=" " vcf_input} /cromwell_root/ssp_vcfs/

        time \
        gcloud storage cp ~{sep=" " tbi_input} /cromwell_root/ssp_vcfs/

        # then merge, and safely assume all ssp-VCFs are sorted in the same order, on one chr
        cd ssp_vcfs
        ls *.vcf.gz > my_vcfs.txt

        bcftools merge \
            --threads ~{threads_num} \
            --merge none \
            -l my_vcfs.txt \
            -O z \
            -o ~{pref}.AllSamples.vcf.gz

        tabix -@ ~{threads_num} -p vcf ~{pref}.AllSamples.vcf.gz

        # move result files to the correct location for cromwell to de-localize
        mv ~{pref}.AllSamples.vcf.gz ~{pref}.AllSamples.vcf.gz.tbi /cromwell_root/
    >>>

    output{
        File merged_vcf = "~{pref}.AllSamples.vcf.gz"
        File merged_tbi = "~{pref}.AllSamples.vcf.gz.tbi"
    }

    runtime {
        cpu: 16
        memory: "32 GiB"
        disks: "local-disk 375 LOCAL"
        preemptible: 1
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}
