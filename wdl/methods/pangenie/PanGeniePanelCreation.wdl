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

workflow PanGeniePanelCreation {
    input {
        File phased_bcf
        File reference_fasta
        File prepare_vcf_script
        File add_ids_script
        File merge_vcfs_script
        Float frac_missing = 0.2
        String output_prefix

        String docker
        File? monitoring_script
    }

    call PanGeniePanelCreation {
        input:
            phased_bcf = phased_bcf,
            reference_fasta = reference_fasta,
            prepare_vcf_script = prepare_vcf_script,
            add_ids_script = add_ids_script,
            merge_vcfs_script = merge_vcfs_script,
            frac_missing = frac_missing,
            output_prefix = output_prefix,
            docker = docker,
            monitoring_script = monitoring_script
    }

    output {
        File panel_vcf_gz = PanGeniePanelCreation.panel_vcf_gz
        File panel_vcf_gz_tbi = PanGeniePanelCreation.panel_vcf_gz_tbi
    }
}


# TODO consider piping more steps
task PanGeniePanelCreation {
    input {
        File phased_bcf
        File reference_fasta
        String output_prefix

        File prepare_vcf_script
        File add_ids_script
        File merge_vcfs_script
        Float frac_missing

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -euxo pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # validate variants against reference
        bcftools norm --check-ref e --fasta-ref ~{reference_fasta} ~{phased_bcf} &> validate-vcf.log

        # atomize variants and run PanGenie prepare-vcf script
        bcftools norm -a ~{phased_bcf} | \
            python3 ~{prepare_vcf_script} \
                --missing ~{frac_missing} \
            2> prepare-vcf.log \
            1> prepare.vcf

        bcftools stats prepare.vcf > ~{output_prefix}.prepare.stats.txt

        # run PanGenie add-ids script
        cat prepare.vcf | \
            python3 ~{add_ids_script} \
            2> add-ids.log \
            1> prepare.id.vcf

        # split to biallelic
        bcftools norm -m- prepare.id.vcf \
            --threads $(nproc) \
            2> split.log \
            1> prepare.id.split.vcf

        # run PanGenie merge script
        pip install pyfaidx
        python3 ~{merge_vcfs_script} merge \
            -vcf prepare.id.split.vcf \
            -r ~{reference_fasta} \
            -ploidy 2  \
            2> merge-haplotypes.log \
            1> prepare.id.split.mergehap.vcf

        bcftools view prepare.id.split.mergehap.vcf \
            -Oz -o ~{output_prefix}.prepare.id.split.mergehap.vcf.gz
        bcftools index -t ~{output_prefix}.prepare.id.split.mergehap.vcf.gz

        bcftools stats ~{output_prefix}.prepare.id.split.mergehap.vcf.gz > ~{output_prefix}.prepare.id.split.mergehap.stats.txt
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        Array[File] logs = glob("*.log")
        File prepare_stats = "~{output_prefix}.prepare.stats.txt"
        File panel_stats = "~{output_prefix}.prepare.id.split.mergehap.stats.txt"
        File panel_vcf_gz = "~{output_prefix}.prepare.id.split.mergehap.vcf.gz"
        File panel_vcf_gz_tbi = "~{output_prefix}.prepare.id.split.mergehap.vcf.gz.tbi"
    }
}
