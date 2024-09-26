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

workflow KAGEPlusGLIMPSEBatchedCase {
    input {
        Array[File] input_crams
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
        Array[File]+ genetic_maps
        Boolean subset_reads = true
        Array[String] sample_names
        Float average_coverage
        String output_prefix

        String docker
        String kage_docker
        File? monitoring_script
    }

    scatter (i in range(length(input_crams))) {
        call IndexCaseReads {
            # TODO we require the alignments to subset by chromosome; change to start from raw reads
            input:
                input_cram = input_crams[i],
                docker = docker,
                monitoring_script = monitoring_script
        }

        call KAGECountKmers {
            input:
                input_cram = input_crams[i],
                input_cram_idx = IndexCaseReads.cram_idx,
                panel_kmer_index_only_variants_with_revcomp = panel_kmer_index_only_variants_with_revcomp,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                reference_dict = reference_dict,
                chromosomes = chromosomes,
                subset_reads = subset_reads,
                output_prefix = sample_names[i],
                docker = kage_docker,
                monitoring_script = monitoring_script
        }

        call KAGEGenotype {
            input:
                kmer_counts = KAGECountKmers.kmer_counts,
                panel_index = panel_index,
                panel_multi_split_vcf_gz = panel_multi_split_vcf_gz,
                panel_multi_split_vcf_gz_tbi = panel_multi_split_vcf_gz_tbi,
                output_prefix = sample_names[i],
                sample_name = sample_names[i],
                average_coverage = average_coverage,
                docker = kage_docker,
                monitoring_script = monitoring_script
        }
    }

    call MergeSamples {
        input:
            input_bcfs = KAGEGenotype.kage_bcf,
            output_prefix = output_prefix,
            monitoring_script = monitoring_script
    }

    scatter (j in range(length(chromosomes))) {
        call GLIMPSECaseChromosome as GLIMPSEBatchedCaseChromosome {
            input:
                kage_bcf = MergeSamples.merged_bcf,
                panel_split_vcf_gz = panel_split_vcf_gz,
                panel_split_vcf_gz_tbi = panel_split_vcf_gz_tbi,
                reference_fasta_fai = reference_fasta_fai,
                chromosome = chromosomes[j],
                genetic_map = genetic_maps[j],
                output_prefix = output_prefix,
                docker = kage_docker,
                monitoring_script = monitoring_script
        }
    }

    call GLIMPSECase as GLIMPSEBatchedCase {
        input:
            chromosome_glimpse_vcf_gzs = GLIMPSEBatchedCaseChromosome.chromosome_glimpse_vcf_gz,
            chromosome_glimpse_vcf_gz_tbis = GLIMPSEBatchedCaseChromosome.chromosome_glimpse_vcf_gz_tbi,
            output_prefix = output_prefix,
            docker = kage_docker,
            monitoring_script = monitoring_script
    }

    output {
        Array[File] cram_idxs = IndexCaseReads.cram_idx
        Array[File] kmer_counts = KAGECountKmers.kmer_counts
        Array[File] kage_bcfs = KAGEGenotype.kage_bcf
        File kage_bcf = MergeSamples.merged_bcf
        File glimpse_unphased_vcf_gz = GLIMPSEBatchedCase.glimpse_unphased_vcf_gz
        File glimpse_unphased_vcf_gz_tbi = GLIMPSEBatchedCase.glimpse_unphased_vcf_gz_tbi
        File glimpse_bcf = GLIMPSEBatchedCase.glimpse_vcf_gz
        File glimpse_vcf_gz_tbi = GLIMPSEBatchedCase.glimpse_vcf_gz_tbi
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

        Boolean? ignore_helper_model = true
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
            ~{true='-I true' false='-I false' ignore_helper_model} \
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

        bcftools concat --no-version -a ~{output_prefix}.kage.bi.vcf.gz ~{output_prefix}.multi.split.vcf.gz -Ob -o ~{output_prefix}.kage.bcf
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
        File kage_bcf = "~{output_prefix}.kage.bcf"
    }
}

task MergeSamples {
    parameter_meta {
        input_bcfs: {localization_optional: true}
    }

    input {
        Array[File] input_bcfs
        String output_prefix

        File? monitoring_script
    }

    command <<<
        set -eou pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # we do single-sample BCF localization ourselves
        mkdir -p input_bcfs
        time \
        gcloud storage cp ~{sep=" " input_bcfs} /cromwell_root/input_bcfs/

        # then merge, and safely assume all single-sample VCFs are sorted in the same order, on one chr
        ls /cromwell_root/input_bcfs/*.bcf > input_bcfs.txt

        bcftools merge \
            --threads $(nproc) \
            --merge none \
            -l input_bcfs.txt \
            -Ob -o ~{output_prefix}.merged.bcf
    >>>

    output {
        File monitoring_log = "monitoring.log"
        File merged_bcf = "~{output_prefix}.merged.bcf"
    }

    runtime {
        cpu: 8
        memory: "32 GiB"
        disks: "local-disk 250 LOCAL"
        preemptible: 2
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}

task GLIMPSECaseChromosome {
    input {
        File kage_bcf
        File panel_split_vcf_gz       # for GLIMPSE
        File panel_split_vcf_gz_tbi
        File reference_fasta_fai
        String chromosome
        File genetic_map
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

        bcftools view --no-version -r ~{chromosome} ~{kage_bcf} | \
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
            --map ~{genetic_map} \
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
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # TODO update to GLIMPSE2
        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_ligate_static
        chmod +x GLIMPSE_ligate_static

        ./GLIMPSE_ligate_static \
            --input ~{write_lines(chromosome_glimpse_vcf_gzs)} \
            --output ~{output_prefix}.ligate.bcf \
            --log ~{output_prefix}.ligate.log

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_sample_static
        chmod +x GLIMPSE_sample_static

        ./GLIMPSE_sample_static --input ~{output_prefix}.ligate.bcf \
            --solve \
            --output ~{output_prefix}.sample.bcf \
            --log ~{output_prefix}.sample.log

        bcftools view --no-version ~{output_prefix}.ligate.bcf -Oz -o ~{output_prefix}.kage.glimpse.unphased.vcf.gz
        bcftools index -t ~{output_prefix}.kage.glimpse.unphased.vcf.gz

        bcftools view --no-version ~{output_prefix}.sample.bcf -Oz -o ~{output_prefix}.kage.glimpse.vcf.gz
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
        File ligate_log = "~{output_prefix}.ligate.log"
        File glimpse_unphased_vcf_gz = "~{output_prefix}.kage.glimpse.unphased.vcf.gz"
        File glimpse_unphased_vcf_gz_tbi = "~{output_prefix}.kage.glimpse.unphased.vcf.gz.tbi"
        File glimpse_vcf_gz = "~{output_prefix}.kage.glimpse.vcf.gz"
        File glimpse_vcf_gz_tbi = "~{output_prefix}.kage.glimpse.vcf.gz.tbi"
    }
}
