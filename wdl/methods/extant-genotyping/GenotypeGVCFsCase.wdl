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

workflow GenotypeGVCFsCase {
    input {
        File gvcf
        File gvcf_tbi
        File? intervals
        File? intervals_tbi
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String output_prefix

        Boolean use_bcftools = true

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
        RuntimeAttributes runtime_attributes = {}
        File? monitoring_script
    }

    call GenotypeGVCFs {
        input:
            gvcf = gvcf,
            gvcf_tbi = gvcf_tbi,
            intervals = intervals,
            intervals_tbi = intervals_tbi,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            output_prefix = output_prefix,
            use_bcftools = use_bcftools,
            gatk_docker = gatk_docker,
            runtime_attributes = runtime_attributes,
            monitoring_script = monitoring_script
    }

    output {
        File genotyped_vcf_gz = GenotypeGVCFs.genotyped_vcf_gz
        File genotyped_vcf_gz_tbi = GenotypeGVCFs.genotyped_vcf_gz_tbi
    }
}

# modified from https://github.com/broadinstitute/warp/blob/b83303ee552c3d925ca40d5213b3219dd7fb1967/tasks/broad/JointGenotypingTasks.wdl#L141
task GenotypeGVCFs {

    input {
        File gvcf
        File gvcf_tbi
        File? intervals
        File? intervals_tbi
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String output_prefix
        Boolean use_bcftools
        String? extra_args = "--allow-old-rms-mapping-quality-annotation-data"     # This is needed for gVCFs generated with GATK3 HaplotypeCaller

        String gatk_docker
        RuntimeAttributes runtime_attributes = {}

        File? monitoring_script
    }

    command <<<
        set -euo pipefail

        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        if [ ~{use_bcftools} ]; then
            bcftools norm --threads $(nproc) -m-any ~{"-R " + intervals} ~{gvcf} | \
                bcftools annotate -x ^FORMAT/GT,^FORMAT/GQ,^FORMAT/PL,QUAL,INFO -e 'ALT="<NON_REF>"' \
                    -Oz -o ~{output_prefix}.vcf.gz
            bcftools index -t ~{output_prefix}.vcf.gz
        else
            gatk --java-options "-Xmx~{default=6 runtime_attributes.command_mem_gb}G" \
                GenotypeGVCFs \
                -V ~{gvcf} \
                ~{"-L " + intervals} \
                -R ~{ref_fasta} \
                -O ~{output_prefix}.vcf.gz \
                -G StandardAnnotation \
                ~{extra_args}
        fi
    >>>

    runtime {
        docker: gatk_docker
        cpu: select_first([runtime_attributes.cpu, 2])
        memory: select_first([runtime_attributes.command_mem_gb, 3]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 10]) + if select_first([runtime_attributes.use_ssd, true]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File? monitoring_log = "monitoring.log"
        File genotyped_vcf_gz = "~{output_prefix}.vcf.gz"
        File genotyped_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }
}
