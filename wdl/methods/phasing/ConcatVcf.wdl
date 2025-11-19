version 1.0

workflow ConcatVCF {
    input {
        Array[File] vcfs
        Array[File]? vcf_tbis
        String prefix
        String output_dir


    }
    call bcftools_concat {
        input:
            vcfs = vcfs,
            vcf_tbis = vcf_tbis,
            prefix = prefix,
            output_dir = output_dir
    }


    output {
    }
}

task bcftools_concat {
    input {
        Array[File] vcfs
        Array[File]? vcf_tbis
        String prefix
        String output_dir
    }

    command <<<
        set -euxo pipefail

        # Index all input VCF files if no precomputed tbis provided.
        if ! ~{defined(vcf_tbis)}; then
            for vcf in ~{sep=" " vcfs}; do
                bcftools index "$vcf"
            done
        fi

        bcftools concat \
            ~{sep=" " vcfs} \
            -n \
            -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz

        gsutil cp ~{prefix}.vcf.gz ~{output_dir}
        gsutil cp ~{prefix}.vcf.gz.tbi ~{output_dir}
    >>>

    output {
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 1000 SSD"
    }
}
