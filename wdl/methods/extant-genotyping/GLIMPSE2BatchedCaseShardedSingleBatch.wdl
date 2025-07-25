version 1.0

import "../../methods/phasing/HierarchicallyMergeVcfs.wdl" as HierarchicallyMergeVcfs

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

workflow GLIMPSE2BatchedCaseShardedSingleBatch {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        Array[String] sample_names

        # per chromosome
        Array[String]+ chromosomes
        Array[File]+ genetic_maps
        Array[File] panel_split_vcf_gz
        Array[File] panel_split_vcf_gz_tbi

        String extra_chunk_args = "--thread $(nproc) --window-mb 5 --buffer-mb 0.5 --sequential"
        String extra_phase_args = "--impute-reference-only-variants --keep-monomorphic-ref-sites"
        String output_prefix

        # inputs for FixVariantCollisions
        File fix_variant_collisions_java
        Int operation
        String weight_tag
        Int is_weight_format_field

        String docker
        File? monitoring_script

        RuntimeAttributes concat_runtime_attributes = {"use_ssd": true}
        RuntimeAttributes glimpse2_phase_runtime_attributes = {}
        RuntimeAttributes glimpse2_sample_runtime_attributes = {}
        Map[String, Int]? chromosome_to_glimpse2_command_mem_gb      # for running per-chromosome by choosing large chunk size; this will override glimpse2_phase_runtime_attributes
    }

    scatter (j in range(length(chromosomes))) {
        String chromosome = chromosomes[j]

        call GLIMPSE2Chunk as ChromosomeGLIMPSE2Chunk {
            input:
                vcf = panel_split_vcf_gz[j],
                tbi = panel_split_vcf_gz_tbi[j],
                region = chromosome,
                prefix = output_prefix + "." + chromosome,
                extra_chunk_args = extra_chunk_args
        }

        Array[String] input_regions = read_lines(ChromosomeGLIMPSE2Chunk.input_regions)
        Array[String] output_regions = read_lines(ChromosomeGLIMPSE2Chunk.output_regions)

        if (defined(chromosome_to_glimpse2_command_mem_gb)) {
            Int command_mem_gb = select_first([chromosome_to_glimpse2_command_mem_gb])[chromosome]
        }
        scatter (k in range(length(input_regions))) {
            call GLIMPSE2Phase as ChunkedGLIMPSE2Phase {
                input:
                    input_vcf_gz = input_vcf_gz,
                    input_vcf_gz_tbi = input_vcf_gz_tbi,
                    panel_split_vcf_gz = panel_split_vcf_gz[j],
                    panel_split_vcf_gz_tbi = panel_split_vcf_gz_tbi[j],
                    input_region = input_regions[k],
                    output_region = output_regions[k],
                    genetic_map = genetic_maps[j],
                    output_prefix = output_prefix + "." + chromosome + ".shard-" + k + ".phased",
                    extra_phase_args = extra_phase_args,
                    docker = docker,
                    monitoring_script = monitoring_script,
                    runtime_attributes = glimpse2_phase_runtime_attributes,
                    command_mem_gb = command_mem_gb
            }
        }

        call GLIMPSE2Ligate as ChromosomeGLIMPSE2Ligate {
            input:
                vcfs = ChunkedGLIMPSE2Phase.phased_vcf_gz,
                vcf_idxs = ChunkedGLIMPSE2Phase.phased_vcf_gz_tbi,
                prefix = output_prefix + "." + chromosome + ".ligated",
                docker = docker,
                monitoring_script = monitoring_script
        }

        call GLIMPSE2Sample as ChromosomeGLIMPSE2Sample {
            input:
                ligated_vcf_gz = ChromosomeGLIMPSE2Ligate.ligated_vcf_gz,
                ligated_vcf_gz_tbi = ChromosomeGLIMPSE2Ligate.ligated_vcf_gz_tbi,
                output_prefix = output_prefix + "." + chromosome + ".sampled",
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = glimpse2_sample_runtime_attributes
        }

        call FixVariantCollisions as ChromosomePhasedCollisionless { input:
            vcf_gz = ChromosomeGLIMPSE2Sample.glimpse2_vcf_gz,
            vcf_gz_tbi = ChromosomeGLIMPSE2Sample.glimpse2_vcf_gz_tbi,
            fix_variant_collisions_java = fix_variant_collisions_java,
            operation = operation,
            weight_tag = weight_tag,
            is_weight_format_field = is_weight_format_field,
            output_prefix = output_prefix + "." + chromosome + ".glimpse2.collisionless",
            monitoring_script = monitoring_script
        }
    }

    call ConcatVcfs as GLIMPSE2PosteriorsConcatVcfs {
        input:
            vcf_gzs = ChromosomeGLIMPSE2Ligate.ligated_vcf_gz,
            vcf_gz_tbis = ChromosomeGLIMPSE2Ligate.ligated_vcf_gz_tbi,
            output_prefix = output_prefix + ".glimpse2.posteriors",
            docker = docker,
            monitoring_script = monitoring_script
    }

    call ConcatVcfs as GLIMPSE2ConcatVcfs {
        input:
            vcf_gzs = ChromosomeGLIMPSE2Sample.glimpse2_vcf_gz,
            vcf_gz_tbis = ChromosomeGLIMPSE2Sample.glimpse2_vcf_gz_tbi,
            output_prefix = output_prefix + ".glimpse2",
            docker = docker,
            monitoring_script = monitoring_script
    }

    call ConcatVcfs as PhasedCollisionlessConcatVcfs {
        input:
            vcf_gzs = ChromosomePhasedCollisionless.collisionless_vcf_gz,
            vcf_gz_tbis = ChromosomePhasedCollisionless.collisionless_vcf_gz_tbi,
            output_prefix = output_prefix + ".glimpse2.collisionless",
            docker = docker,
            monitoring_script = monitoring_script
    }

    output {
        File glimpse_posteriors_vcf_gz = GLIMPSE2PosteriorsConcatVcfs.vcf_gz
        File glimpse_posteriors_vcf_gz_tbi = GLIMPSE2PosteriorsConcatVcfs.vcf_gz_tbi
        File glimpse_vcf_gz = GLIMPSE2ConcatVcfs.vcf_gz
        File glimpse_vcf_gz_tbi = GLIMPSE2ConcatVcfs.vcf_gz_tbi
        File phased_collisionless_vcf_gz = PhasedCollisionlessConcatVcfs.vcf_gz
        File phased_collisionless_vcf_gz_tbi = PhasedCollisionlessConcatVcfs.vcf_gz_tbi
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

        # TODO FIX LEXICOGRAPHICAL BUG!
        mkdir inputs
        mv ~{sep=' ' vcf_gzs} inputs
        mv ~{sep=' ' vcf_gz_tbis} inputs

        if [ $(ls inputs/*.vcf.gz | wc -l) == 1 ]
        then
            cp $(ls inputs/*.vcf.gz) ~{output_prefix}.vcf.gz
            cp $(ls inputs/*.vcf.gz.tbi) ~{output_prefix}.vcf.gz.tbi
        else
            bcftools concat $(ls inputs/*.vcf.gz) --naive -Oz -o ~{output_prefix}.vcf.gz
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

task GLIMPSE2Chunk {
    input {
        File vcf
        File tbi
        String region
        String prefix
        String? extra_chunk_args

        RuntimeAttributes runtime_attributes = {}
    }

    Int disk_size_gb = 2*ceil(size([vcf, tbi], "GB")) + 1

    command <<<
        set -euxo pipefail

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_chunk_static
        chmod +x GLIMPSE2_chunk_static

        ./GLIMPSE2_chunk_static \
            -I ~{vcf} \
            --region ~{region} \
            ~{extra_chunk_args} \
            -O chunks.txt

        # cut chunks + buffers
        cut -f 3 chunks.txt > input-regions.txt
        cut -f 4 chunks.txt > output-regions.txt
    >>>

    output {
        File input_regions = "input-regions.txt"
        File output_regions = "output-regions.txt"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
        cpu: select_first([runtime_attributes.cpu, 4])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, disk_size_gb]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 10])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }
}

task GLIMPSE2Phase {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File panel_split_vcf_gz
        File panel_split_vcf_gz_tbi
        String input_region
        String output_region
        File genetic_map
        String output_prefix
        String? extra_phase_args

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
        Int? command_mem_gb = 7
    }

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # TODO skeep only biallelic SNVs for now
        bcftools view --no-version -r ~{input_region},~{output_region} -m2 -M2 -v snps ~{input_vcf_gz} \
            -Ob -o ~{output_prefix}.biSNV.bcf
        bcftools index ~{output_prefix}.biSNV.bcf

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_phase_static
        chmod +x GLIMPSE2_phase_static

        ./GLIMPSE2_phase_static \
            -I ~{output_prefix}.biSNV.bcf \
            -R ~{panel_split_vcf_gz} \
            --input-region ~{input_region} \
            --output-region ~{output_region} \
            --map ~{genetic_map} \
            --input-GL \
            --thread $(nproc) \
            ~{extra_phase_args} \
            --output ~{output_prefix}.raw.vcf.gz

        # take input VCF header and add GLIMPSE INFO and FORMAT lines (GLIMPSE header only contains a single chromosome and breaks bcftools concat --naive)
        bcftools view --no-version -h ~{input_vcf_gz} | grep '^##' > input.header.txt
        bcftools view --no-version -h ~{output_prefix}.raw.vcf.gz | grep -E '^##INFO|^##FORMAT|^##NMAIN|^##FPLOIDY' > glimpse2.header.txt
        bcftools view --no-version -h ~{input_vcf_gz} | grep '^#CHROM' > input.columns.txt
        cat input.header.txt glimpse2.header.txt input.columns.txt > header.txt
        bcftools reheader -h header.txt ~{output_prefix}.raw.vcf.gz > ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, command_mem_gb]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File phased_vcf_gz = "~{output_prefix}.vcf.gz"          # note that this actually contains posteriors */* GTs; these are converted to phased *|* GTs later on in the sample step
        File phased_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }
}

task GLIMPSE2Ligate {
    input {
        Array[File] vcfs
        Array[File]? vcf_idxs
        String prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    Int disk_size_gb = 2*ceil(size(vcfs, "GB")) + 1

    command <<<
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        if ! ~{defined(vcf_idxs)}; then
            for ff in ~{sep=' ' vcfs}; do bcftools index $ff; done
        fi

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_ligate_static
        chmod +x GLIMPSE2_ligate_static

        ./GLIMPSE2_ligate_static --input ~{write_lines(vcfs)} --output ~{prefix}.vcf.gz --thread $(nproc)
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File monitoring_log = "monitoring.log"
        File ligated_vcf_gz = "~{prefix}.vcf.gz"
        File ligated_vcf_gz_tbi = "~{prefix}.vcf.gz.tbi"
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 2])
        memory: select_first([runtime_attributes.command_mem_gb, 7]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, disk_size_gb]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 10])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }
}

task GLIMPSE2Sample {
    input {
        File ligated_vcf_gz
        File ligated_vcf_gz_tbi
        String output_prefix

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

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_sample_static
        chmod +x GLIMPSE2_sample_static

        ./GLIMPSE2_sample_static --input ~{ligated_vcf_gz} \
            --solve \
            --output ~{output_prefix}.vcf.gz \
            --log ~{output_prefix}.sample.log
        bcftools index -t ~{output_prefix}.vcf.gz
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

    output {
        File monitoring_log = "monitoring.log"
        File glimpse2_vcf_gz = "~{output_prefix}.vcf.gz"
        File glimpse2_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }
}

task FixVariantCollisions {

    input {
        File vcf_gz                         # biallelic
        File vcf_gz_tbi
        File fix_variant_collisions_java
        Int operation = 1                   # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "UNIT_WEIGHT"   # ID of the weight field; if this field is not found, all weights are set to one; weights are assumed to be non-negative
        Int is_weight_format_field = 0      # given a VCF record in a sample, assign it a weight encoded in the sample column (1) or in the INFO field (0)
        String output_prefix

        File? monitoring_script
    }

    command <<<
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        java ~{fix_variant_collisions_java} \
            ~{vcf_gz} \
            ~{operation} \
            ~{weight_tag} \
            ~{is_weight_format_field} \
            collisionless.vcf \
            windows.txt \
            histogram.txt \
            null                            # do not output figures

        # replace all missing alleles (correctly) emitted with reference alleles, since this is expected by PanGenie panel-creation script
         bcftools view collisionless.vcf | \
            sed -e 's/\.|0/0|0/g' | sed -e 's/0|\./0|0/g' | sed -e 's/\.|1/0|1/g' | sed -e 's/1|\./1|0/g' | sed -e 's/\.|\./0|0/g' | \
            bcftools view -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    >>>

    output {
        File monitoring_log = "monitoring.log"
        File collisionless_vcf_gz = "~{output_prefix}.vcf.gz"
        File collisionless_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
        File windows = "windows.txt"
        File histogram = "histogram.txt"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "16 GiB"
        disks: "local-disk 100 SSD"
        bootDiskSizeGb: 10
        preemptible:     3
        max_retries:           2
        docker:"us.gcr.io/broad-gatk/gatk:4.6.0.0"
    }
}
