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

workflow KAGECase {
    input {
        File input_cram
        File panel_index
        File panel_kmer_index_only_variants_with_revcomp
        File panel_split_vcf_gz # for GLIMPSE
        File panel_split_vcf_gz_tbi
        File panel_multi_split_vcf_gz # for filling in biallelic-only VCFs produced by KAGE
        File panel_multi_split_vcf_gz_tbi
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        Array[String]+ chromosomes
        Boolean subset_reads = true
        String sample_name
        Float average_coverage

        String docker
        String kage_docker
        File? monitoring_script
    }

    call IndexCaseReads {
        # TODO we require the alignments to subset by chromosome; change to start from raw reads
        input:
            input_cram = input_cram,
            docker = docker,
            monitoring_script = monitoring_script
    }

    call KAGECountKmers {
        input:
            input_cram = input_cram,
            input_cram_idx = IndexCaseReads.cram_idx,
            panel_kmer_index_only_variants_with_revcomp = panel_kmer_index_only_variants_with_revcomp,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            reference_dict = reference_dict,
            chromosomes = chromosomes,
            subset_reads = subset_reads,
            output_prefix = sample_name,
            docker = kage_docker,
            monitoring_script = monitoring_script
    }

    call KAGEGenotype {
        input:
            kmer_counts = KAGECountKmers.kmer_counts,
            panel_index = panel_index,
            panel_multi_split_vcf_gz = panel_multi_split_vcf_gz,
            panel_multi_split_vcf_gz_tbi = panel_multi_split_vcf_gz_tbi,
            output_prefix = sample_name,
            sample_name = sample_name,
            average_coverage = average_coverage,
            docker = kage_docker,
            monitoring_script = monitoring_script
    }

    scatter (j in range(length(chromosomes))) {
        call GLIMPSECaseChromosome {
            input:
                kage_vcf_gz = KAGEGenotype.kage_vcf_gz,
                kage_vcf_gz_tbi = KAGEGenotype.kage_vcf_gz_tbi,
                panel_split_vcf_gz = panel_split_vcf_gz,
                panel_split_vcf_gz_tbi = panel_split_vcf_gz_tbi,
                reference_fasta_fai = reference_fasta_fai,
                chromosome = chromosomes[j],
                output_prefix = sample_name,
                docker = kage_docker,
                monitoring_script = monitoring_script
        }
    }

    call GLIMPSECaseGather {
        input:
            chromosome_glimpse_vcf_gzs = GLIMPSECaseChromosome.chromosome_glimpse_vcf_gz,
            chromosome_glimpse_vcf_gz_tbis = GLIMPSECaseChromosome.chromosome_glimpse_vcf_gz_tbi,
            output_prefix = sample_name,
            docker = kage_docker,
            monitoring_script = monitoring_script
    }

    output {
        File kmer_counts = KAGECountKmers.kmer_counts
        File kage_vcf_gz = KAGEGenotype.kage_vcf_gz
        File kage_vcf_gz_tbi = KAGEGenotype.kage_vcf_gz_tbi
        File glimpse_vcf_gz = GLIMPSECaseGather.glimpse_vcf_gz
        File glimpse_vcf_gz_tbi = GLIMPSECaseGather.glimpse_vcf_gz_tbi
    }
}

# some 1000G CRAM indices have issues due to htsjdk version, so we reindex; see e.g. https://github.com/broadinstitute/gatk/issues/7076
task IndexCaseReads {
    input {
        File input_cram

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    # if we instead simply use File cram_idx = "~{input_cram}.crai" in the output block,
    # Terra tries to localize an index adjacent to the CRAM in PreprocessCaseReads???
    String output_prefix = basename(input_cram, ".cram")

    command {
        set -eou pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        samtools index -@ $(nproc) ~{input_cram} ~{output_prefix}.cram.crai
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 2])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, true]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File cram_idx = "~{output_prefix}.cram.crai"
    }
}

task KAGECountKmers {
    input {
        File input_cram
        File input_cram_idx
        File panel_kmer_index_only_variants_with_revcomp
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        Array[String]+ chromosomes
        Boolean subset_reads
        String output_prefix

        String docker
        File? monitoring_script

        String kmer_mapper_args = "-c 10000000"

        RuntimeAttributes runtime_attributes = {}
        Int? cpu = 8
    }

    # explicitly specify CPU and memory separately for KAGE commands;
    # using nproc w/ automatic CPU/memory scaling of custom instances on Terra can be suboptimal or cause OOM
    Int cpu_resolved = select_first([runtime_attributes.cpu, cpu])

    command <<<
        set -eou pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        date
        mkfifo ~{output_prefix}.preprocessed.fa
        if ~{subset_reads}; then
            # hacky way to get chromosomes into bed file
            grep -P '~{sep="\\t|" chromosomes}\t' ~{reference_fasta_fai} | cut -f 1,2 | sed -e 's/\t/\t1\t/g' > chromosomes.bed

            echo "Subsetting reads..."
            samtools view --reference ~{reference_fasta} -@ $(nproc) ~{if subset_reads then "--regions-file chromosomes.bed" else ""} -u -X ~{input_cram} ~{input_cram_idx} | \
                samtools fasta --reference ~{reference_fasta} -@ $(nproc) > ~{output_prefix}.preprocessed.fa &
        else
            echo "Not subsetting reads..."
            samtools fasta --reference ~{reference_fasta} -@ $(nproc) -X ~{input_cram} ~{input_cram_idx} > ~{output_prefix}.preprocessed.fa &
        fi

        kmer_mapper map \
            ~{kmer_mapper_args} \
            -t ~{cpu_resolved} \
            -i ~{panel_kmer_index_only_variants_with_revcomp} \
            -f ~{output_prefix}.preprocessed.fa \
            -o ~{output_prefix}.kmer_counts.npy

        rm ~{output_prefix}.preprocessed.fa
    >>>

    runtime {
        docker: docker
        cpu: cpu_resolved
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File kmer_counts = "~{output_prefix}.kmer_counts.npy"
    }
}

task KAGEGenotype {
    input {
        File kmer_counts
        File panel_index
        File panel_multi_split_vcf_gz # for filling in biallelic-only VCFs produced by KAGE
        File panel_multi_split_vcf_gz_tbi
        String output_prefix
        String sample_name
        Float average_coverage

        String docker
        File? monitoring_script

        Boolean? ignore_helper_model = false
        String? kage_genotype_extra_args

        RuntimeAttributes runtime_attributes = {}
        Int? cpu = 8
    }

    # explicitly specify CPU and memory separately for KAGE commands;
    # using nproc w/ automatic CPU/memory scaling of custom instances on Terra can be suboptimal or cause OOM
    Int cpu_resolved = select_first([runtime_attributes.cpu, cpu])

    command {
        set -eou pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        kage genotype \
            -i ~{panel_index} \
            -c ~{kmer_counts} \
            --average-coverage ~{average_coverage} \
            -s ~{sample_name} \
            ~{true='-I' false='' ignore_helper_model} \
            ~{kage_genotype_extra_args} \
            -o ~{output_prefix}.kage.bi.vcf

        # we need to add split multiallelics to biallelic-only KAGE VCF
        # create single-sample header from LOO panel w/ split multiallelics
        bcftools view --no-version -h -G ~{panel_multi_split_vcf_gz} | \
            sed 's/##fileformat=VCFv4.2/##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods.">/g' | \
            sed 's/INFO$/INFO\tFORMAT\t~{sample_name}/g' > ~{output_prefix}.multi.split.header.txt
        # create single-sample missing genotypes from LOO panel w/ split multiallelics
        bcftools view --no-version -H -G ~{panel_multi_split_vcf_gz} | \
            sed 's/$/\tGT:GL\t.\/.:nan,nan,nan/g' > ~{output_prefix}.multi.split.GT.txt
        # create single-sample VCF w/ split multiallelics
        bgzip -c <(cat ~{output_prefix}.multi.split.header.txt ~{output_prefix}.multi.split.GT.txt) > ~{output_prefix}.multi.split.vcf.gz
        bcftools index -t ~{output_prefix}.multi.split.vcf.gz

        bgzip -c ~{output_prefix}.kage.bi.vcf > ~{output_prefix}.kage.bi.vcf.gz
        bcftools index -t ~{output_prefix}.kage.bi.vcf.gz

        bcftools concat --no-version -a ~{output_prefix}.kage.bi.vcf.gz ~{output_prefix}.multi.split.vcf.gz -Oz -o ~{output_prefix}.kage.vcf.gz
        bcftools index -t ~{output_prefix}.kage.vcf.gz
    }

    runtime {
        docker: docker
        cpu: cpu_resolved
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 10]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File kage_vcf_gz = "~{output_prefix}.kage.vcf.gz"
        File kage_vcf_gz_tbi = "~{output_prefix}.kage.vcf.gz.tbi"
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
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -eou pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools view --no-version -r ~{chromosome} ~{kage_vcf_gz} | \
            sed -e 's/nan/-1000000.0/g' | sed -e 's/-inf/-1000000.0/g' | sed -e 's/inf/-1000000.0/g' | bgzip > ~{output_prefix}.kage.nonan.~{chromosome}.vcf.gz
        bcftools index -t ~{output_prefix}.kage.nonan.~{chromosome}.vcf.gz

        # TODO update to GLIMPSE2; first figure out why it complains about AC/AN and GT being inconsistent?
        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_phase_static
        chmod +x GLIMPSE_phase_static

        CHROMOSOME_LENGTH=$(grep -P "~{chromosome}\t" ~{reference_fasta_fai} | cut -f 2)
        ./GLIMPSE_phase_static \
            -I ~{output_prefix}.kage.nonan.~{chromosome}.vcf.gz \
            -R ~{panel_split_vcf_gz} \
            --input-region ~{chromosome}:1-$CHROMOSOME_LENGTH \
            --output-region ~{chromosome}:1-$CHROMOSOME_LENGTH \
            --input-GL \
            --thread $(nproc) \
            --output ~{output_prefix}.kage.glimpse.~{chromosome}.vcf.gz
        bcftools index -t ~{output_prefix}.kage.glimpse.~{chromosome}.vcf.gz
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
        File chromosome_glimpse_vcf_gz = "~{output_prefix}.kage.glimpse.~{chromosome}.vcf.gz"
        File chromosome_glimpse_vcf_gz_tbi = "~{output_prefix}.kage.glimpse.~{chromosome}.vcf.gz.tbi"
    }
}

task GLIMPSECaseGather {
    input {
        Array[File] chromosome_glimpse_vcf_gzs
        Array[File] chromosome_glimpse_vcf_gz_tbis
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -eou pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools concat --no-version ~{sep=" " chromosome_glimpse_vcf_gzs} -Oz -o ~{output_prefix}.kage.glimpse.vcf.gz
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
        File glimpse_vcf_gz = "~{output_prefix}.kage.glimpse.vcf.gz"
        File glimpse_vcf_gz_tbi = "~{output_prefix}.kage.glimpse.vcf.gz.tbi"
    }
}