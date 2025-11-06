version 1.0

import "Helper.wdl"

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

workflow StatisticalPhasingFromPanGenie {

    input {
        File pangenie_vcf_gz        # cohort, whole genome, multiallelic
        File pangenie_vcf_gz_tbi
        Array[String]+ chromosomes
        String output_prefix

        String docker
        File? monitoring_script

        String? extra_chunk_args
        File genetic_mapping_tsv_for_shapeit4
        Int shapeit4_num_threads
        Int shapeit4_memory
        String shapeit4_extra_args
    }

    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)

    scatter (j in range(length(chromosomes))) {
        String chromosome = chromosomes[j]
        
        call CreateShapeit4Chunks { input:
            vcf = pangenie_vcf_gz,
            tbi = pangenie_vcf_gz_tbi,
            region = chromosome,
            prefix = output_prefix + "." + chromosome,
            extra_chunk_args = extra_chunk_args
        }

        Array[String] region_list = read_lines(CreateShapeit4Chunks.chunks)

        scatter (j in range(length(region_list))) {
            call Helper.Shapeit4 as Shapeit4 { input:
                vcf_input = pangenie_vcf_gz,
                vcf_index = pangenie_vcf_gz_tbi,
                mappingfile = genetic_mapping_dict[chromosome],
                region = region_list[j],
                prefix = output_prefix + "." + chromosome + ".shard-" + j + ".phased",
                num_threads = shapeit4_num_threads,
                memory = shapeit4_memory,
                extra_args = shapeit4_extra_args
            }
        }

        call LigateVcfs { input:
            vcfs = Shapeit4.phased_bcf,
            prefix = output_prefix + "." + chromosome + ".phased.ligated"
        }
    }

    call ConcatVcfs {
        input:
            vcf_gzs = LigateVcfs.ligated_vcf_gz,
            vcf_gz_tbis = LigateVcfs.ligated_vcf_gz_tbi,
            output_prefix = output_prefix + ".phased",
            docker = docker,
            monitoring_script = monitoring_script
    }


    output {
        File phased_vcf_gz = ConcatVcfs.vcf_gz
        File phased_vcf_gz_tbi = ConcatVcfs.vcf_gz_tbi
    }
}

task CreateShapeit4Chunks {

    input {
        File vcf
        File tbi
        String region
        String prefix
        String? extra_chunk_args = "--thread $(nproc) --window-size 5000000 --buffer-size 500000"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf, tbi], "GB")) + 1

    command <<<
        set -euxo pipefail

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_chunk_static
        chmod +x GLIMPSE_chunk_static

        ./GLIMPSE_chunk_static \
            -I ~{vcf} \
            --region ~{region} \
            ~{extra_chunk_args} \
            -O chunks.txt

        # cut chunks + buffers
        cut -f 3 chunks.txt > chunks.regions.txt
    >>>

    output {
        File chunks = "chunks.regions.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task LigateVcfs {

    input {
        Array[File] vcfs
        Array[File]? vcf_idxs
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(vcfs, "GB")) + 1

    command <<<
        set -euxo pipefail
        if ! ~{defined(vcf_idxs)}; then
            for ff in ~{sep=' ' vcfs}; do bcftools index $ff; done
        fi

        wget https://github.com/odelaneau/shapeit5/releases/download/v5.1.1/ligate_static
        chmod +x ligate_static

        ./ligate_static --input ~{write_lines(vcfs)} --output ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File ligated_vcf_gz = "~{prefix}.vcf.gz"
        File ligated_vcf_gz_tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ConcatVcfs {
    input {
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {"use_ssd": true}
    }

    Int disk_size_gb = 3 * ceil(size(vcf_gzs, "GB"))

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        mkdir inputs
        mv ~{sep=' ' vcf_gzs} inputs
        mv ~{sep=' ' vcf_gz_tbis} inputs

        # TODO NOTE USE OF SORT, ENSURE THIS GIVES DESIRED CHROMOSOME ORDER (COULD TAKE IN ARRAY OF CHROMOSOMES INSTEAD)
        if [ $(ls inputs/*.vcf.gz | wc -l) == 1 ]
        then
            cp $(ls inputs/*.vcf.gz) ~{output_prefix}.vcf.gz
            cp $(ls inputs/*.vcf.gz.tbi) ~{output_prefix}.vcf.gz.tbi
        else
            bcftools concat $(ls inputs/*.vcf.gz | sort -V -d) --naive -Oz -o ~{output_prefix}.vcf.gz
            bcftools index -t ~{output_prefix}.vcf.gz
        fi
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
