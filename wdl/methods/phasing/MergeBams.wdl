version 1.0


workflow MergeBams {

    input {

        Array[File] bams
        Array[File] bais
        String prefix
        Int threads_num
    }

    call MergeBams { input:
        bam_input = bams,
        bai_input = bais,
        prefix = prefix,
        threads_num = threads_num
    }

    output {
        File merged_bam = MergeBams.merged_bam
        File merged_bai = MergeBams.merged_bai
    }
}


task MergeBams {
    parameter_meta {
        bam_input: {localization_optional: true}
        bai_input: {localization_optional: true}
    }
    input{
        Array[File] bam_input
        Array[File] bai_input
        String prefix
        Int threads_num
    }

    command <<<
        set -eux

        # we do single-sample phased VCFs localization ourselves
        mkdir -p bams
        time \
        gcloud storage cp ~{sep=" " bam_input} /mnt/disks/cromwell_root/bams/

        time \
        gcloud storage cp ~{sep=" " bai_input} /mnt/disks/cromwell_root/bams/

        # then merge, and safely assume all ssp-VCFs are sorted in the same order, on one chr
        cd bams
        ls *.bam > my_bams.txt

        samtools merge \
            --threads ~{threads_num} \
            -b my_bams.txt \
            -O BAM \
            -o ~{prefix}.bam
        samtools index ~{prefix}.bam

    >>>

    output{
        File merged_bam = "bams/~{prefix}.bam"
        File merged_bai = "bams/~{prefix}.bam.bai"
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