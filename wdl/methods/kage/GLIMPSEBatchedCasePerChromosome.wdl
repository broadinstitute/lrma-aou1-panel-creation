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

workflow GLIMPSEBatchedCasePerChromosome {
    input {
        File sample_by_chromosome_kage_vcf_gzs_tsv
        File sample_by_chromosome_kage_vcf_gz_tbis_tsv
        File sample_names_file
        File reference_fasta
        File reference_fasta_fai
        File reference_dict

        # per chromosome
        Array[String]+ chromosomes
        Array[File]+ genetic_maps
        Array[File] panel_split_vcf_gz # for GLIMPSE
        Array[File] panel_split_vcf_gz_tbi
        Map[String, Int] chromosome_to_glimpse_command_mem_gb

        Int batch_size
        String output_prefix

        # inputs for FixVariantCollisions
        File fix_variant_collisions_java
        Int operation
        String weight_tag
        Int is_weight_format_field

        String kage_docker
        File? monitoring_script

        RuntimeAttributes glimpse_case_chromosome_runtime_attributes = {}
        RuntimeAttributes glimpse_case_runtime_attributes = {}
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

            call GLIMPSECaseChromosome as GLIMPSEBatchedCaseChromosome {
                input:
                    kage_vcf_gz = ChromosomeKAGEMergeAcrossSamples.merged_vcf_gz,
                    kage_vcf_gz_tbi = ChromosomeKAGEMergeAcrossSamples.merged_vcf_gz_tbi,
                    panel_split_vcf_gz = panel_split_vcf_gz[j],
                    panel_split_vcf_gz_tbi = panel_split_vcf_gz_tbi[j],
                    reference_fasta_fai = reference_fasta_fai,
                    chromosome = chromosome,
                    genetic_map = genetic_maps[j],
                    output_prefix = output_prefix + ".batch-" + b,
                    docker = kage_docker,
                    monitoring_script = monitoring_script,
                    command_mem_gb = chromosome_to_glimpse_command_mem_gb[chromosome],
                    runtime_attributes = glimpse_case_chromosome_runtime_attributes
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

        call GLIMPSECase as GLIMPSEBatchedCase {
            input:
                chromosome_glimpse_vcf_gzs = GLIMPSEBatchedCaseChromosome.chromosome_glimpse_vcf_gz,
                chromosome_glimpse_vcf_gz_tbis = GLIMPSEBatchedCaseChromosome.chromosome_glimpse_vcf_gz_tbi,
                output_prefix = output_prefix + ".batch-" + b,
                docker = kage_docker,
                monitoring_script = monitoring_script,
                runtime_attributes = glimpse_case_runtime_attributes
        }

        call FixVariantCollisions as GLIMPSEFixVariantCollisions { input:
            vcf_gz = GLIMPSEBatchedCase.glimpse_vcf_gz,
            vcf_gz_tbi = GLIMPSEBatchedCase.glimpse_vcf_gz_tbi,
            fix_variant_collisions_java = fix_variant_collisions_java,
            operation = operation,
            weight_tag = weight_tag,
            is_weight_format_field = is_weight_format_field,
            output_prefix = output_prefix + ".batch-" + b + ".kage.glimpse.collisionless",
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

    call Ivcfmerge as GLIMPSEUnphasedMergeAcrossSamples { input:
        vcf_gzs = GLIMPSEBatchedCase.glimpse_unphased_vcf_gz,
        vcf_gz_tbis = GLIMPSEBatchedCase.glimpse_unphased_vcf_gz_tbi,
        sample_names = sample_names,
        output_prefix = output_prefix + ".kage.glimpse.unphased",
        docker = kage_docker,
        monitoring_script = monitoring_script
    }

    call Ivcfmerge as GLIMPSEMergeAcrossSamples { input:
        vcf_gzs = GLIMPSEBatchedCase.glimpse_vcf_gz,
        vcf_gz_tbis = GLIMPSEBatchedCase.glimpse_vcf_gz_tbi,
        sample_names = sample_names,
        output_prefix = output_prefix + ".kage.glimpse",
        docker = kage_docker,
        monitoring_script = monitoring_script
    }

    call Ivcfmerge as PhasedCollisionlessMergeAcrossSamples { input:
        vcf_gzs = GLIMPSEFixVariantCollisions.collisionless_vcf_gz,
        vcf_gz_tbis = GLIMPSEFixVariantCollisions.collisionless_vcf_gz_tbi,
        sample_names = sample_names,
        output_prefix = output_prefix + ".kage.glimpse.collisionless",
        docker = kage_docker,
        monitoring_script = monitoring_script
    }

    output {
        Array[File] batch_kage_vcf_gzs = KAGEConcatVcfs.vcf_gz
        Array[File] batch_kage_vcf_gz_tbis = KAGEConcatVcfs.vcf_gz_tbi
        Array[File] batch_glimpse_unphased_vcf_gzs = GLIMPSEBatchedCase.glimpse_unphased_vcf_gz
        Array[File] batch_glimpse_unphased_vcf_gz_tbis = GLIMPSEBatchedCase.glimpse_unphased_vcf_gz_tbi
        Array[File] batch_glimpse_vcf_gzs = GLIMPSEBatchedCase.glimpse_vcf_gz
        Array[File] batch_glimpse_vcf_gz_tbis = GLIMPSEBatchedCase.glimpse_vcf_gz_tbi
        Array[File] batch_phased_collisionless_vcf_gzs = GLIMPSEFixVariantCollisions.collisionless_vcf_gz
        Array[File] batch_phased_collisionless_vcf_gz_tbis = GLIMPSEFixVariantCollisions.collisionless_vcf_gz_tbi

        File kage_vcf_gz = KAGEMergeAcrossSamples.merged_vcf_gz
        File kage_vcf_gz_tbi = KAGEMergeAcrossSamples.merged_vcf_gz_tbi
        File glimpse_unphased_vcf_gz = GLIMPSEUnphasedMergeAcrossSamples.merged_vcf_gz
        File glimpse_unphased_vcf_gz_tbi = GLIMPSEUnphasedMergeAcrossSamples.merged_vcf_gz_tbi
        File glimpse_vcf_gz = GLIMPSEMergeAcrossSamples.merged_vcf_gz
        File glimpse_vcf_gz_tbi = GLIMPSEMergeAcrossSamples.merged_vcf_gz_tbi
        File phased_collisionless_vcf_gz = PhasedCollisionlessMergeAcrossSamples.merged_vcf_gz
        File phased_collisionless_vcf_gz_tbi = PhasedCollisionlessMergeAcrossSamples.merged_vcf_gz_tbi
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

task GLIMPSECaseChromosome {
    input {
        File kage_vcf_gz
        File kage_vcf_gz_tbi
        File panel_split_vcf_gz       # for GLIMPSE
        File panel_split_vcf_gz_tbi
        File reference_fasta_fai
        String chromosome
        File genetic_map
        String output_prefix

        String docker
        File? monitoring_script

        Int command_mem_gb
        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools view --no-version -r ~{chromosome} ~{kage_vcf_gz} | \
            sed -e 's/nan/-1000000.0/g' | sed -e 's/-inf/-1000000.0/g' | sed -e 's/inf/-1000000.0/g' | \
            bcftools view -Ob -o ~{output_prefix}.kage.nonan.~{chromosome}.bcf
        bcftools index ~{output_prefix}.kage.nonan.~{chromosome}.bcf

        # TODO update to GLIMPSE2; first figure out why it complains about AC/AN and GT being inconsistent?
        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_phase_static
        chmod +x GLIMPSE_phase_static

        CHROMOSOME_LENGTH=$(grep -P "~{chromosome}\t" ~{reference_fasta_fai} | cut -f 2)
        ./GLIMPSE_phase_static \
            -I ~{output_prefix}.kage.nonan.~{chromosome}.bcf \
            -R ~{panel_split_vcf_gz} \
            --input-region ~{chromosome}:1-$CHROMOSOME_LENGTH \
            --output-region ~{chromosome}:1-$CHROMOSOME_LENGTH \
            --map ~{genetic_map} \
            --input-GL \
            --thread $(nproc) \
            --output ~{output_prefix}.~{chromosome}.kage.glimpse.raw.vcf.gz

        # take KAGE VCF header and add GLIMPSE INFO and FORMAT lines (GLIMPSE header only contains a single chromosome and breaks bcftools concat --naive)
        bcftools view --no-version -h ~{kage_vcf_gz} | grep '^##' > kage.header.txt
        bcftools view --no-version -h ~{output_prefix}.~{chromosome}.kage.glimpse.raw.vcf.gz | grep -E '^##INFO|^##FORMAT|^##NMAIN|^##FPLOIDY' > glimpse.header.txt
        bcftools view --no-version -h ~{kage_vcf_gz} | grep '^#CHROM' > kage.columns.txt
        cat kage.header.txt glimpse.header.txt kage.columns.txt > header.txt
        bcftools reheader -h header.txt ~{output_prefix}.~{chromosome}.kage.glimpse.raw.vcf.gz > ~{output_prefix}.~{chromosome}.kage.glimpse.vcf.gz
        bcftools index -t ~{output_prefix}.~{chromosome}.kage.glimpse.vcf.gz
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
        File chromosome_glimpse_vcf_gz = "~{output_prefix}.~{chromosome}.kage.glimpse.vcf.gz"
        File chromosome_glimpse_vcf_gz_tbi = "~{output_prefix}.~{chromosome}.kage.glimpse.vcf.gz.tbi"
    }
}

task GLIMPSECase {
    input {
        Array[File] chromosome_glimpse_vcf_gzs
        Array[File] chromosome_glimpse_vcf_gz_tbis
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

        bcftools concat \
            -f ~{write_lines(chromosome_glimpse_vcf_gzs)} \
            --naive \
            -Oz -o ~{output_prefix}.kage.glimpse.unphased.vcf.gz
        bcftools index -t ~{output_prefix}.kage.glimpse.unphased.vcf.gz

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_sample_static
        chmod +x GLIMPSE_sample_static

        ./GLIMPSE_sample_static --input ~{output_prefix}.kage.glimpse.unphased.vcf.gz \
            --solve \
            --output ~{output_prefix}.kage.glimpse.vcf.gz \
            --log ~{output_prefix}.sample.log
        bcftools index -t ~{output_prefix}.kage.glimpse.vcf.gz
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
        File glimpse_unphased_vcf_gz = "~{output_prefix}.kage.glimpse.unphased.vcf.gz"
        File glimpse_unphased_vcf_gz_tbi = "~{output_prefix}.kage.glimpse.unphased.vcf.gz.tbi"
        File glimpse_vcf_gz = "~{output_prefix}.kage.glimpse.vcf.gz"
        File glimpse_vcf_gz_tbi = "~{output_prefix}.kage.glimpse.vcf.gz.tbi"
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
        disks: "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-gatk/gatk:4.6.0.0"
    }
}
