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
        File? input_crai
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

        String kage_docker
        File? monitoring_script

        RuntimeAttributes kage_runtime_attributes = {}
    }

    if (!defined(input_crai)) {
        call IndexCaseReads {
            # TODO we require the alignments to subset by chromosome; change to start from raw reads
            input:
                input_cram = input_cram,
                monitoring_script = monitoring_script
        }
    }

    call KAGE {
        input:
            input_cram = input_cram,
            input_cram_idx = select_first([input_crai, IndexCaseReads.cram_idx]),
            panel_kmer_index_only_variants_with_revcomp = panel_kmer_index_only_variants_with_revcomp,
            panel_index = panel_index,
            panel_multi_split_vcf_gz = panel_multi_split_vcf_gz,
            panel_multi_split_vcf_gz_tbi = panel_multi_split_vcf_gz_tbi,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            reference_dict = reference_dict,
            chromosomes = chromosomes,
            subset_reads = true,
            average_coverage = average_coverage,
            output_prefix = sample_name,
            docker = kage_docker,
            monitoring_script = monitoring_script,
            runtime_attributes = kage_runtime_attributes
    }

    output {
        Array[File] chromosome_kmer_counts = KAGE.chromosome_kmer_counts
        Array[File] chromosome_kage_vcf_gzs = KAGE.chromosome_kage_vcf_gzs
        Array[File] chromosome_kage_vcf_gz_tbis = KAGE.chromosome_kage_vcf_gz_tbis
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
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File cram_idx = "~{output_prefix}.cram.crai"
    }
}

task KAGE {
    input {
        File input_cram
        File input_cram_idx
        Array[File] panel_kmer_index_only_variants_with_revcomp
        Array[File] panel_index
        Array[File] panel_multi_split_vcf_gz # for filling in biallelic-only VCFs produced by KAGE
        Array[File] panel_multi_split_vcf_gz_tbi
        String output_prefix
        String sample_name
        Float average_coverage
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        Array[String]+ chromosomes
        Boolean subset_reads
        String output_prefix

        String docker
        File? monitoring_script

        String kmer_mapper_args = "-c 10000000"

        Boolean? ignore_helper_model = true
        String? kage_genotype_extra_args

        RuntimeAttributes runtime_attributes = {}
        Int? cpu = 8
    }

    # explicitly specify CPU and memory separately for KAGE commands;
    # using nproc w/ automatic CPU/memory scaling of custom instances on Terra can be suboptimal or cause OOM
    Int cpu_resolved = select_first([runtime_attributes.cpu, cpu])

    Int num_chromosomes = length(chromosomes)

    command <<<
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        CHROMOSOMES=(~{sep=" " chromosomes})
        KMER_INDICES=(~{sep=" " panel_kmer_index_only_variants_with_revcomp})
        INDICES=(~{sep=" " panel_index})
        MULTI_SPLITS=(~{sep=" " panel_multi_split_vcf_gz})

        NUM_CHROMOSOMES=~{num_chromosomes}
        NUM_DIGITS=${#NUM_CHROMOSOMES}

        for (( c = 0; c < ~{num_chromosomes}; c++ ))
        do
            CHROMOSOME=${CHROMOSOMES[$c]}
            KMER_INDEX=${KMER_INDICES[$c]}
            INDEX=${INDICES[$c]}
            MULTI_SPLIT=${MULTI_SPLITS[$c]}

            C_WITH_LEADING_ZEROS=$(printf "%0${NUM_DIGITS}d" $c)

            date
            mkfifo ~{output_prefix}.preprocessed.fa
            if ~{subset_reads}; then
                # hacky way to get chromosomes into bed file
                grep -P $CHROMOSOME ~{reference_fasta_fai} | cut -f 1,2 | sed -e 's/\t/\t1\t/g' > chromosome.bed

                echo "Subsetting reads..."
                samtools view --reference ~{reference_fasta} -@ $(nproc) ~{if subset_reads then "--regions-file chromosome.bed" else ""} -u -X ~{input_cram} ~{input_cram_idx} | \
                    samtools fasta --reference ~{reference_fasta} -@ $(nproc) > ~{output_prefix}.preprocessed.fa &
            else
                echo "Not subsetting reads..."
                samtools fasta --reference ~{reference_fasta} -@ $(nproc) -X ~{input_cram} ~{input_cram_idx} > ~{output_prefix}.preprocessed.fa &
            fi

            kmer_mapper map \
                ~{kmer_mapper_args} \
                -t ~{cpu_resolved} \
                -i $KMER_INDEX \
                -f ~{output_prefix}.preprocessed.fa \
                -o ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.kmer_counts.npy

            rm ~{output_prefix}.preprocessed.fa

            kage genotype \
                -i $INDEX \
                -c ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.kmer_counts.npy \
                --average-coverage ~{average_coverage} \
                -s ~{sample_name} \
                ~{true='-I true' false='-I false' ignore_helper_model} \
                ~{kage_genotype_extra_args} \
                -o ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.kage.bi.vcf
            bcftools view ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.kage.bi.vcf --write-index -Ob -o ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.kage.bi.bcf

            # we need to add split multiallelics to biallelic-only KAGE BCF
            # create single-sample header from LOO panel w/ split multiallelics
            bcftools view --no-version -h -G $MULTI_SPLIT | \
                sed 's/##fileformat=VCFv4.2/##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods.">/g' | \
                sed 's/INFO$/INFO\tFORMAT\t~{sample_name}/g' > ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.multi.split.header.txt
            # create single-sample missing genotypes from LOO panel w/ split multiallelics
            bcftools view --no-version -H -G $MULTI_SPLIT | \
                sed 's/$/\tGT:GL\t.\/.:nan,nan,nan/g' > ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.multi.split.GT.txt
            # create single-sample BCF w/ split multiallelics
            bcftools view <(cat ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.multi.split.header.txt ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.multi.split.GT.txt) --write-index -Ob -o ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.multi.split.bcf

            bcftools concat --no-version -a ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.kage.bi.bcf ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.multi.split.bcf -Oz -o ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.kage.vcf.gz
            bcftools index -t ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.kage.vcf.gz

            rm ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.multi.split.header.txt ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.multi.split.GT.txt ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.kage.bi.bcf ~{output_prefix}.shard-$C_WITH_LEADING_ZEROS.$CHROMOSOME.multi.split.bcf
        done
    >>>

    runtime {
        docker: docker
        cpu: cpu_resolved
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 150]) + if select_first([runtime_attributes.use_ssd, true]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        Array[File] chromosome_kmer_counts = glob("~{output_prefix}.*.kmer_counts.npy")
        Array[File] chromosome_kage_vcf_gzs = glob("~{output_prefix}.*.kage.vcf.gz")
        Array[File] chromosome_kage_vcf_gz_tbis = glob("~{output_prefix}.*.kage.vcf.gz.tbi")
    }
}
