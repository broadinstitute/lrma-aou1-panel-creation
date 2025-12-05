version 1.0

import "./Helper.wdl" as H

workflow SplitCohortVcf {

    input {
        File joint_short_vcf
        File joint_short_vcf_tbi
        String outputdirectory
    }

    call SplitVcf { input:
        joint_vcf = joint_short_vcf,
        joint_vcf_tbi = joint_short_vcf_tbi
    }

    call H.FinalizeToDir { input:
        files = SplitVcf.vcf_by_sample,
        outdir = outputdirectory
    }

    output {
        Array[File] splitted_vcf = SplitVcf.vcf_by_sample
    }
}


task SplitVcf {

    input {
        File joint_vcf
        File joint_vcf_tbi
        Int memory
    }

    command <<<
        set -euxo pipefail
        mkdir output
        bcftools +split -Oz -o output ~{joint_vcf}
        # cd output
        # for vcf in $(find . -name "*.vcf.gz"); do
        #     tabix -p vcf "$vcf"
        # done
        # cd -

    >>>

    output {
        Array[File] vcf_by_sample = glob("output/*vcf.gz")
        # Array[File] vcf_by_sample_tbi = glob("output/*vcf.gz.tbi")

    }
    ###################
    runtime {
        cpu: 4
        memory: memory + " GiB"
        disks: "local-disk 1000 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}