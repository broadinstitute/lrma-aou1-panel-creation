version 1.0


workflow TRGTMerge {

    input {
        Array[File] trgt_vcfs
        Array[File] trgt_vcfs_tbi
        String prefix
    }


    call TRGTMerge{ input:
        vcfs = trgt_vcfs,
        vcfs_tbi = trgt_vcfs_tbi,
        prefix = prefix + ".trgtmerged"
    }


    output {
        File merged_vcf = TRGTMerge.merged_vcf
        File merged_tbi = TRGTMerge.merged_tbi
    }
}

task TRGTMerge {
    input {
        Array[File] vcfs
        Array[File] vcfs_tbi
        String prefix
    }

    command <<<
        set -euxo pipefail
        wget https://github.com/PacificBiosciences/trgt/releases/download/v1.5.0/trgt-v1.5.0-x86_64-unknown-linux-gnu.tar.gz
        tar -xvf trgt-v1.5.0-x86_64-unknown-linux-gnu.tar.gz
         # use sep to seperate vcfs array
        ./trgt-v1.5.0-x86_64-unknown-linux-gnu/trgt merge --vcf ~{sep=" " vcfs} -o ~{prefix}.vcf.gz 
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_tbi = "~{prefix}.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 2
        memory:  "8 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}
