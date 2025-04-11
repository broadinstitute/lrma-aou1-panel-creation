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

workflow KAGECasePerChromosome {
    input {
        File input_cram
        String sample_name
        File reference_fasta
        File reference_fasta_fai
        File reference_dict

        # per chromosome
        Array[String]+ chromosomes
        Array[File] panel_index
        Array[File] panel_kmer_index_only_variants_with_revcomp
        Array[File] panel_multi_split_vcf_gz # for filling in biallelic-only VCFs produced by KAGE
        Array[File] panel_multi_split_vcf_gz_tbi

        Float average_coverage

        String samtools_docker
        String kage_docker
        File? monitoring_script

        RuntimeAttributes kage_count_kmers_runtime_attributes = {}
        RuntimeAttributes kage_genotype_runtime_attributes = {}
    }

    call IndexCaseReads {
        # TODO we require the alignments to subset by chromosome; change to start from raw reads
        input:
            input_cram = input_cram,
            docker = samtools_docker,
            monitoring_script = monitoring_script
    }

    scatter (j in range(length(chromosomes))) {
        String chromosome = chromosomes[j]

        call KAGECountKmers {
            input:
                input_cram = input_cram,
                input_cram_idx = IndexCaseReads.cram_idx,
                panel_kmer_index_only_variants_with_revcomp = panel_kmer_index_only_variants_with_revcomp[j],
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                reference_dict = reference_dict,
                chromosomes = [chromosome],
                subset_reads = true,
                output_prefix = sample_name + "." + chromosome,
                docker = kage_docker,
                monitoring_script = monitoring_script,
                runtime_attributes = kage_count_kmers_runtime_attributes
        }

        call KAGEGenotype {
            input:
                kmer_counts = KAGECountKmers.kmer_counts,
                panel_index = panel_index[j],
                panel_multi_split_vcf_gz = panel_multi_split_vcf_gz[j],
                panel_multi_split_vcf_gz_tbi = panel_multi_split_vcf_gz_tbi[j],
                output_prefix = sample_name + "." + chromosome,
                sample_name = sample_name,
                average_coverage = average_coverage,
                docker = kage_docker,
                monitoring_script = monitoring_script,
                runtime_attributes = kage_genotype_runtime_attributes
        }
    }

    output {
        File cram_idx = IndexCaseReads.cram_idx
        Array[File] chromosome_kmer_counts = KAGECountKmers.kmer_counts
        Array[File] chromosome_kage_vcf_gzs = KAGEGenotype.kage_vcf_gz
        Array[File] chromosome_kage_vcf_gz_tbis = KAGEGenotype.kage_vcf_gz_tbi
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
        set -euox pipefail

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
        preemptible: select_first([runtime_attributes.preemptible, 3])
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
        set -euox pipefail

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
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, true]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 3])
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

        Boolean? ignore_helper_model = true
        String? kage_genotype_extra_args

        RuntimeAttributes runtime_attributes = {}
        Int? cpu = 8
    }

    # explicitly specify CPU and memory separately for KAGE commands;
    # using nproc w/ automatic CPU/memory scaling of custom instances on Terra can be suboptimal or cause OOM
    Int cpu_resolved = select_first([runtime_attributes.cpu, cpu])

    command {
        set -euox pipefail

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
            ~{true='-I true' false='-I false' ignore_helper_model} \
            ~{kage_genotype_extra_args} \
            -o ~{output_prefix}.kage.bi.vcf
        bcftools view ~{output_prefix}.kage.bi.vcf -Ob -o ~{output_prefix}.kage.bi.bcf
        bcftools index ~{output_prefix}.kage.bi.bcf

        # we need to add split multiallelics to biallelic-only KAGE BCF
        # create single-sample header from LOO panel w/ split multiallelics
        bcftools view --no-version -h -G ~{panel_multi_split_vcf_gz} | \
            sed 's/##fileformat=VCFv4.2/##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods.">/g' | \
            sed 's/INFO$/INFO\tFORMAT\t~{sample_name}/g' > ~{output_prefix}.multi.split.header.txt
        # create single-sample missing genotypes from LOO panel w/ split multiallelics
        bcftools view --no-version -H -G ~{panel_multi_split_vcf_gz} | \
            sed 's/$/\tGT:GL\t.\/.:nan,nan,nan/g' > ~{output_prefix}.multi.split.GT.txt
        # create single-sample BCF w/ split multiallelics
        bcftools view <(cat ~{output_prefix}.multi.split.header.txt ~{output_prefix}.multi.split.GT.txt) -Ob -o ~{output_prefix}.multi.split.bcf
        bcftools index ~{output_prefix}.multi.split.bcf

        bcftools concat --no-version -a ~{output_prefix}.kage.bi.bcf ~{output_prefix}.multi.split.bcf -Oz -o ~{output_prefix}.kage.vcf.gz
        bcftools index -t ~{output_prefix}.kage.vcf.gz
    }

    runtime {
        docker: docker
        cpu: cpu_resolved
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 10]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 3])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File kage_vcf_gz = "~{output_prefix}.kage.vcf.gz"
        File kage_vcf_gz_tbi = "~{output_prefix}.kage.vcf.gz.tbi"
    }
}
