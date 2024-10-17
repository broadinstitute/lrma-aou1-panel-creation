version 1.0

struct RuntimeAttributes {
    Int? cpu
    Int? command_mem_gb
    Int? additional_mem_gb
    Int? disk_size_gb
    Int? boot_disk_size_gb
    Boolean? use_ssd
    Int? preemptible
    Int? max_retries
}

workflow HierarchicallyMergeVcfs {
    input {
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        Int batch_size
        String output_prefix
        String? extra_args

        String docker
        File? monitoring_script
    }

    call CreateBatches {
        input:
            vcf_gzs = vcf_gzs,
            vcf_gz_tbis = vcf_gz_tbis,
            batch_size = batch_size,
            docker = docker
    }

    scatter (i in range(length(CreateBatches.vcf_gz_batch_fofns))) {
        call MergeVcfs as MergeVcfsSingleBatch {
            input:
                vcf_gzs = read_lines(CreateBatches.vcf_gz_batch_fofns[i]),
                vcf_gz_tbis = read_lines(CreateBatches.vcf_gz_tbi_batch_fofns[i]),
                output_prefix = output_prefix + "." + i,
                extra_args = extra_args,
                docker = docker,
                monitoring_script = monitoring_script
        }
    }

    call MergeVcfs as MergeVcfsAllBatches {
        input:
            vcf_gzs = MergeVcfsSingleBatch.merged_vcf_gz,
            vcf_gz_tbis = MergeVcfsSingleBatch.merged_vcf_gz_tbi,
            output_prefix = output_prefix,
            extra_args = extra_args,
            docker = docker,
            monitoring_script = monitoring_script
    }

    output {
        File merged_vcf_gz = MergeVcfsAllBatches.merged_vcf_gz
        File merged_vcf_gz_tbi = MergeVcfsAllBatches.merged_vcf_gz_tbi
    }
}

task CreateBatches {
    input {
        Array[String] vcf_gzs
        Array[String] vcf_gz_tbis
        Int batch_size

        String docker
        RuntimeAttributes runtime_attributes = {}
    }

    command <<<
        set -euox pipefail

        cat ~{write_lines(vcf_gzs)} | split -l ~{batch_size} - vcf_gz_batch_
        cat ~{write_lines(vcf_gz_tbis)} | split -l ~{batch_size} - vcf_gz_tbi_batch_
    >>>

    output{
        Array[File] vcf_gz_batch_fofns = glob("vcf_gz_batch_*")
        Array[File] vcf_gz_tbi_batch_fofns = glob("vcf_gz_tbi_batch_*")
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }
}

task MergeVcfs {
    input{
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        String output_prefix
        String? extra_args

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools merge \
            -l ~{write_lines(vcf_gzs)} \
            ~{extra_args} \
            -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    }

    output {
        File monitoring_log = "monitoring.log"
        File merged_vcf_gz = "~{output_prefix}.vcf.gz"
        File merged_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }
}
