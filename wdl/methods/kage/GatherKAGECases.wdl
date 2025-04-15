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

workflow GatherKAGECases {
    input {
        File sample_by_chromosome_kage_vcf_gzs_tsv
        File sample_by_chromosome_kage_vcf_gz_tbis_tsv
        File sample_names_file

        # per chromosome
        Array[String]+ chromosomes

        Int batch_size
        String output_prefix

        String kage_docker
        File? monitoring_script
    }

    Array[Array[String]] sample_by_chromosome_kage_vcf_gzs = read_tsv(sample_by_chromosome_kage_vcf_gzs_tsv)
    Array[Array[String]] sample_by_chromosome_kage_vcf_gzs_kage_vcf_gz_tbis = read_tsv(sample_by_chromosome_kage_vcf_gz_tbis_tsv)
    Array[String] sample_names = read_lines(sample_names_file)

    call CreateBatches {
        input:
            sample_by_chromosome_vcf_gzs = sample_by_chromosome_kage_vcf_gzs,
            sample_by_chromosome_vcf_gz_tbis = sample_by_chromosome_kage_vcf_gzs_kage_vcf_gz_tbis,
            sample_names = sample_names,
            batch_size = batch_size,
            docker = kage_docker
    }

    scatter (b in range(length(CreateBatches.chromosome_vcf_gz_batch_files))) {
        scatter (j in range(length(chromosomes))) {
            String chromosome = chromosomes[j]
            Array[String] chromosome_kage_vcf_gzs = transpose(read_tsv(CreateBatches.chromosome_vcf_gz_batch_files[b]))[j]
            Array[String] chromosome_kage_vcf_gz_tbis = transpose(read_tsv(CreateBatches.chromosome_vcf_gz_tbi_batch_files[b]))[j]

            call Ivcfmerge as ChromosomeKAGEMergeAcrossSamples { input:
                vcf_gzs = chromosome_kage_vcf_gzs,
                vcf_gz_tbis = chromosome_kage_vcf_gz_tbis,
                sample_names = read_lines(CreateBatches.sample_name_batch_files[b]),
                output_prefix = output_prefix + ".batch-" + b + "." + chromosome + ".kage",
                docker = kage_docker,
                monitoring_script = monitoring_script
            }
        }

        call ConcatVcfs as KAGEConcatVcfs {
            input:
                vcf_gzs = ChromosomeKAGEMergeAcrossSamples.merged_vcf_gz,
                vcf_gz_tbis = ChromosomeKAGEMergeAcrossSamples.merged_vcf_gz_tbi,
                output_prefix = output_prefix + ".batch-" + b + ".kage",
                docker = kage_docker,
                monitoring_script = monitoring_script
        }
    }

    call Ivcfmerge as KAGEMergeAcrossSamples { input:
        vcf_gzs = KAGEConcatVcfs.vcf_gz,
        vcf_gz_tbis = KAGEConcatVcfs.vcf_gz_tbi,
        sample_names = sample_names,
        output_prefix = output_prefix + ".kage",
        docker = kage_docker,
        monitoring_script = monitoring_script
    }

    output {
        Array[File] batch_kage_vcf_gzs = KAGEConcatVcfs.vcf_gz
        Array[File] batch_kage_vcf_gz_tbis = KAGEConcatVcfs.vcf_gz_tbi

        File kage_vcf_gz = KAGEMergeAcrossSamples.merged_vcf_gz
        File kage_vcf_gz_tbi = KAGEMergeAcrossSamples.merged_vcf_gz_tbi
    }
}

task CreateBatches {
    input {
        Array[Array[String]] sample_by_chromosome_vcf_gzs
        Array[Array[String]] sample_by_chromosome_vcf_gz_tbis
        Array[String] sample_names
        Int batch_size

        String docker
        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -euox pipefail

        cat ~{write_tsv(sample_by_chromosome_vcf_gzs)} | split -l ~{batch_size} - chromosome_vcf_gz_batch_
        cat ~{write_tsv(sample_by_chromosome_vcf_gz_tbis)} | split -l ~{batch_size} - chromosome_vcf_gz_tbi_batch_
        cat ~{write_lines(sample_names)} | split -l ~{batch_size} - sample_name_batch_
    }

    output {
        Array[File] chromosome_vcf_gz_batch_files = glob("chromosome_vcf_gz_batch_*")
        Array[File] chromosome_vcf_gz_tbi_batch_files = glob("chromosome_vcf_gz_tbi_batch_*")
        Array[File] sample_name_batch_files = glob("sample_name_batch_*")
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 3]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 10]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }
}

# assumes all VCFs have identical variants
task Ivcfmerge {
    input{
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        Array[String] sample_names
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {"use_ssd": true}
    }

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        wget https://github.com/iqbal-lab-org/ivcfmerge/archive/refs/tags/v1.0.0.tar.gz
        tar -xvf v1.0.0.tar.gz

        mkdir decompressed
        cat ~{write_lines(vcf_gzs)} | xargs -I % sh -c 'bcftools annotate --no-version -x INFO % -Ov -o decompressed/$(basename % .gz)'
        time python ivcfmerge-1.0.0/ivcfmerge.py <(ls decompressed/*.vcf) ~{output_prefix}.vcf
        bcftools annotate --no-version -S ~{write_lines(sample_names)} -x FORMAT/FT ~{output_prefix}.vcf -Oz -o ~{output_prefix}.vcf.gz
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
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 250]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }
}

task ConcatVcfs {
    input {
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    Int disk_size_gb = 3 * ceil(size(vcf_gzs, "GB"))

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools concat ~{sep=" " vcf_gzs} --naive -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, disk_size_gb]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File vcf_gz = "~{output_prefix}.vcf.gz"
        File vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }
}
