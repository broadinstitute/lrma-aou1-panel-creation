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
        Array[String] sample_names
        File reference_fasta
        File reference_fasta_fai
        File reference_dict

        # per chromosome
        Array[String]+ chromosomes
        Array[File]+ genetic_maps
        Array[File] panel_index
        Array[File] panel_kmer_index_only_variants_with_revcomp
        Array[File] panel_split_vcf_gz # for GLIMPSE
        Array[File] panel_split_vcf_gz_tbi
        Array[File] panel_multi_split_vcf_gz # for filling in biallelic-only VCFs produced by KAGE
        Array[File] panel_multi_split_vcf_gz_tbi

        Float average_coverage
        String output_prefix

        # inputs for FixVariantCollisions
        File fix_variant_collisions_java
        Int operation
        String weight_tag
        Int is_weight_format_field

        String kage_docker
        File? monitoring_script

        RuntimeAttributes kage_count_kmers_runtime_attributes = {}
        RuntimeAttributes kage_genotype_runtime_attributes = {}
        RuntimeAttributes glimpse_case_chromosome_runtime_attributes = {}
        RuntimeAttributes glimpse_case_runtime_attributes = {}
    }

    scatter (i in range(length(input_crams))) {
        call IndexCaseReads {
            # TODO we require the alignments to subset by chromosome; change to start from raw reads
            input:
                input_cram = input_crams[i],
                monitoring_script = monitoring_script
        }
    }

    scatter (j in range(length(chromosomes))) {
        scatter (i in range(length(input_crams))) {
            call KAGECountKmers {
                input:
                    input_cram = input_crams[i],
                    input_cram_idx = IndexCaseReads.cram_idx[i],
                    panel_kmer_index_only_variants_with_revcomp = panel_kmer_index_only_variants_with_revcomp[j],
                    reference_fasta = reference_fasta,
                    reference_fasta_fai = reference_fasta_fai,
                    reference_dict = reference_dict,
                    chromosomes = [chromosomes[j]],
                    subset_reads = true,
                    output_prefix = sample_names[i] + "." + chromosomes[j],
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
                    output_prefix = sample_names[i] + "." + chromosomes[j],
                    sample_name = sample_names[i],
                    average_coverage = average_coverage,
                    docker = kage_docker,
                    monitoring_script = monitoring_script,
                    runtime_attributes = kage_genotype_runtime_attributes
            }
        }

        call Ivcfmerge as KAGEMergeAcrossSamples { input:
            vcf_gzs = KAGEGenotype.kage_vcf_gz,
            vcf_gz_tbis = KAGEGenotype.kage_vcf_gz_tbi,
            sample_names = sample_names,
            output_prefix = output_prefix + "." + chromosomes[j] + ".kage.merged",
            monitoring_script = monitoring_script
        }

        call GLIMPSECaseChromosome as GLIMPSEBatchedCaseChromosome {
            input:
                kage_vcf_gz = KAGEMergeAcrossSamples.merged_vcf_gz,
                kage_vcf_gz_tbi = KAGEMergeAcrossSamples.merged_vcf_gz_tbi,
                panel_split_vcf_gz = panel_split_vcf_gz[j],
                panel_split_vcf_gz_tbi = panel_split_vcf_gz_tbi[j],
                reference_fasta_fai = reference_fasta_fai,
                chromosome = chromosomes[j],
                genetic_map = genetic_maps[j],
                output_prefix = output_prefix + "." + chromosomes[j],
                docker = kage_docker,
                monitoring_script = monitoring_script,
                runtime_attributes = glimpse_case_chromosome_runtime_attributes
        }
    }

    call ConcatVcfs as KAGEConcatVcfs {
        input:
            vcf_gzs = KAGEMergeAcrossSamples.merged_vcf_gz,
            vcf_gz_tbis = KAGEMergeAcrossSamples.merged_vcf_gz_tbi,
            output_prefix = output_prefix + ".kage",
            monitoring_script = monitoring_script
    }

    call GLIMPSECase as GLIMPSEBatchedCase {
        input:
            chromosome_glimpse_vcf_gzs = GLIMPSEBatchedCaseChromosome.chromosome_glimpse_vcf_gz,
            chromosome_glimpse_vcf_gz_tbis = GLIMPSEBatchedCaseChromosome.chromosome_glimpse_vcf_gz_tbi,
            output_prefix = output_prefix,
            docker = kage_docker,
            monitoring_script = monitoring_script,
            runtime_attributes = glimpse_case_runtime_attributes
    }

    call FixVariantCollisions as GenotypingFixVariantCollisions { input:
        vcf_gz = GLIMPSEBatchedCase.glimpse_vcf_gz,
        vcf_gz_tbi = GLIMPSEBatchedCase.glimpse_vcf_gz_tbi,
        fix_variant_collisions_java = fix_variant_collisions_java,
        operation = operation,
        weight_tag = weight_tag,
        is_weight_format_field = is_weight_format_field,
        output_prefix = output_prefix + ".glimpse.collisionless",
        monitoring_script = monitoring_script
    }

    output {
        File kage_vcf_gz = KAGEConcatVcfs.vcf_gz
        File kage_vcf_gz_tbi = KAGEConcatVcfs.vcf_gz_tbi
        File glimpse_unphased_vcf_gz = GLIMPSEBatchedCase.glimpse_unphased_vcf_gz
        File glimpse_unphased_vcf_gz_tbi = GLIMPSEBatchedCase.glimpse_unphased_vcf_gz_tbi
        File glimpse_vcf_gz = GLIMPSEBatchedCase.glimpse_vcf_gz
        File glimpse_vcf_gz_tbi = GLIMPSEBatchedCase.glimpse_vcf_gz_tbi
        File phased_collisionless_vcf_gz = GenotypingFixVariantCollisions.collisionless_vcf_gz
        File phased_collisionless_vcf_gz_tbi = GenotypingFixVariantCollisions.collisionless_vcf_gz_tbi
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
        bcftools view ~{output_prefix}.kage.bi.vcf --write-index -Ob -o ~{output_prefix}.kage.bi.bcf

        # we need to add split multiallelics to biallelic-only KAGE BCF
        # create single-sample header from LOO panel w/ split multiallelics
        bcftools view --no-version -h -G ~{panel_multi_split_vcf_gz} | \
            sed 's/##fileformat=VCFv4.2/##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods.">/g' | \
            sed 's/INFO$/INFO\tFORMAT\t~{sample_name}/g' > ~{output_prefix}.multi.split.header.txt
        # create single-sample missing genotypes from LOO panel w/ split multiallelics
        bcftools view --no-version -H -G ~{panel_multi_split_vcf_gz} | \
            sed 's/$/\tGT:GL\t.\/.:nan,nan,nan/g' > ~{output_prefix}.multi.split.GT.txt
        # create single-sample BCF w/ split multiallelics
        bcftools view <(cat ~{output_prefix}.multi.split.header.txt ~{output_prefix}.multi.split.GT.txt) --write-index -Ob -o ~{output_prefix}.multi.split.bcf

        bcftools concat --no-version -a ~{output_prefix}.kage.bi.bcf ~{output_prefix}.multi.split.bcf -Oz -o ~{output_prefix}.kage.vcf.gz
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
            bcftools view --write-index -Ob -o ~{output_prefix}.kage.nonan.~{chromosome}.bcf

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
            --output ~{output_prefix}.kage.glimpse.~{chromosome}.raw.vcf.gz

        # take KAGE VCF header and add GLIMPSE INFO and FORMAT lines (GLIMPSE header only contains a single chromosome and breaks bcftools concat --naive)
        bcftools view --no-version -h ~{kage_vcf_gz} | grep '^##' > kage.header.txt
        bcftools view --no-version -h ~{output_prefix}.kage.glimpse.~{chromosome}.raw.vcf.gz | grep -E '^##INFO|^##FORMAT|^##NMAIN|^##FPLOIDY' > glimpse.header.txt
        bcftools view --no-version -h ~{kage_vcf_gz} | grep '^#CHROM' > kage.columns.txt
        cat kage.header.txt glimpse.header.txt kage.columns.txt > header.txt
        bcftools reheader -h header.txt ~{output_prefix}.kage.glimpse.~{chromosome}.raw.vcf.gz > ~{output_prefix}.kage.glimpse.~{chromosome}.vcf.gz
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
