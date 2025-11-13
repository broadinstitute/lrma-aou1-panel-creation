version 1.0

workflow SplitCohortVcf {

    input {
        File joint_vcf
        File joint_vcf_tbi
        Int batch_size
        String outputdirectory
    }

    call get_sample_list { input:
        joint_vcf = joint_vcf,
        joint_vcf_tbi = joint_vcf_tbi,
        batch_size = batch_size
    }

    Array[File] sample_batch = get_sample_list.sample_lists
    scatter (index in range(length(sample_batch))) {
        File sample_list = sample_batch[index]
        call split_sample_vcfs { input:
            joint_vcf = joint_vcf,
            joint_vcf_tbi = joint_vcf_tbi,
            sample_list = sample_list,
            prefix = "sample_batch_" + index,
        }
        call SplitVcf{ input:
            joint_vcf = split_sample_vcfs.subset_vcf,
            joint_vcf_tbi = split_sample_vcfs.subset_tbi,
            gcs_output_dir = outputdirectory
        }
    }

    output {
    }
}


task get_sample_list {

    input {
        File joint_vcf
        File joint_vcf_tbi
        Int batch_size
        Int memory = 16
    }

    command <<<
        set -euxo pipefail
        bcftools query -l ~{joint_vcf} | split -l ~{batch_size} - sample_batch_
    >>>

    output {
        Array[File] sample_lists = glob("sample_batch_*")
    }
    ###################
    runtime {
        cpu: 4
        memory: memory + " GiB"
        disks: "local-disk 100 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}

task split_sample_vcfs {

    input {
        File joint_vcf
        File joint_vcf_tbi
        File sample_list
        
        String prefix
        Int memory = 16
    }

    Array[String] sl = read_lines(sample_list)

    command <<<
        set -euxo pipefail
        bcftools view -s ~{sep="," sl} ~{joint_vcf} -Oz -o ~{prefix}.split.vcf.gz
        bcftools index -t ~{prefix}.split.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.split.vcf.gz"
        File subset_tbi = "~{prefix}.split.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 4
        memory: memory + " GiB"
        disks: "local-disk 100 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}

task SplitVcf {

    input {
        File joint_vcf
        File joint_vcf_tbi
        String gcs_output_dir
        Int memory
    }

    command <<<
        set -euxo pipefail
        mkdir output
        bcftools +split -Oz -o output ~{joint_vcf}
        cd output
        for vcf in $(find . -name "*.vcf.gz"); do
            tabix -p vcf "$vcf"
            gsutil cp "$vcf" ~{gcs_output_dir}
        done
        cd -

    >>>

    output {

    }
    ###################
    runtime {
        cpu: 4
        memory: memory + " GiB"
        disks: "local-disk 500 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}