version 1.0

import "../methods/kage/KAGEPanelWithPreprocessing.wdl" as KAGEPanelWithPreprocessing
import "../methods/pangenie/PanGenieCase.wdl" as PanGenieCase

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

workflow LeaveOutEvaluation {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File case_reference_fasta
        File case_reference_fasta_fai
        File reference_fasta
        File reference_fasta_fai
        Array[File] genetic_maps
        File repeat_mask_bed
        File segmental_duplications_bed
        File simple_repeats_bed
        File challenging_medically_relevant_genes_bed
        String output_prefix
        Array[String] chromosomes
        Array[Array[String]] leave_out_sample_names_array
        Boolean do_pangenie

        # TODO we require the alignments to subset by chromosome; change to start from raw reads
        Map[String, File] leave_out_crams

        String docker
        String kage_docker
        String pangenie_docker
        File? monitoring_script

        RuntimeAttributes? runtime_attributes
        RuntimeAttributes? medium_runtime_attributes
        RuntimeAttributes? large_runtime_attributes
    }

    call PreprocessPanelVCF {
        # TODO move into KAGEPanel?
        input:
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            repeat_mask_bed = repeat_mask_bed,
            segmental_duplications_bed = segmental_duplications_bed,
            simple_repeats_bed = simple_repeats_bed,
            challenging_medically_relevant_genes_bed = challenging_medically_relevant_genes_bed,
            chromosomes = chromosomes,
            output_prefix = output_prefix,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = runtime_attributes
    }

    scatter (i in range(length(leave_out_sample_names_array))) {
        Array[String] leave_out_sample_names = leave_out_sample_names_array[i]
        String leave_out_output_prefix = output_prefix + ".LO-" + i

        call CreateLeaveOneOutPanelVCF {
            input:
                input_vcf_gz = PreprocessPanelVCF.preprocessed_panel_vcf_gz,
                input_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_vcf_gz_tbi,
                output_prefix = leave_out_output_prefix,
                leave_out_sample_names = leave_out_sample_names,
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = runtime_attributes
        }

        call KAGEPanelWithPreprocessing.KAGEPanelWithPreprocessing as KAGELeaveOneOutPanel {
            input:
                input_vcf_gz = CreateLeaveOneOutPanelVCF.leave_out_panel_bi_vcf_gz,
                input_vcf_gz_tbi = CreateLeaveOneOutPanelVCF.leave_out_panel_bi_vcf_gz_tbi,
                do_preprocessing = false,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                output_prefix = leave_out_output_prefix,
                chromosomes = chromosomes,
                docker = kage_docker,
                monitoring_script = monitoring_script
        }

        scatter (j in range(length(leave_out_sample_names))) {
            String leave_out_sample_name = leave_out_sample_names[j]
            String leave_out_cram = leave_out_crams[leave_out_sample_name]

            call IndexCaseReads {
                # TODO we require the alignments to subset by chromosome; change to start from raw reads
                input:
                    input_cram = leave_out_cram,
                    docker = docker,
                    monitoring_script = monitoring_script
            }

            call PreprocessCaseReads {
                # TODO we require the alignments to subset by chromosome; change to start from raw reads
                input:
                    input_cram = leave_out_cram,
                    input_cram_idx = IndexCaseReads.cram_idx,
                    reference_fasta = case_reference_fasta,
                    reference_fasta_fai = case_reference_fasta_fai,
                    output_prefix = leave_out_sample_name,
                    chromosomes = chromosomes,
                    docker = docker,
                    monitoring_script = monitoring_script
            }

            # KAGE+GLIMPSE case
            call KAGECase {
                input:
                    input_fasta = PreprocessCaseReads.preprocessed_fasta,
                    panel_index = KAGELeaveOneOutPanel.index,
                    panel_kmer_index_only_variants_with_revcomp = KAGELeaveOneOutPanel.kmer_index_only_variants_with_revcomp,
                    panel_multi_split_vcf_gz = CreateLeaveOneOutPanelVCF.leave_out_panel_multi_split_vcf_gz,
                    panel_multi_split_vcf_gz_tbi = CreateLeaveOneOutPanelVCF.leave_out_panel_multi_split_vcf_gz_tbi,
                    reference_fasta_fai = reference_fasta_fai,
                    output_prefix = leave_out_sample_name,
                    sample_name = leave_out_sample_name,
                    docker = kage_docker,
                    monitoring_script = monitoring_script
            }

            scatter (k in range(length(chromosomes))) {
                call GLIMPSECaseChromosome {
                    input:
                        kage_vcf_gz = KAGECase.kage_vcf_gz,
                        kage_vcf_gz_tbi = KAGECase.kage_vcf_gz_tbi,
                        panel_split_vcf_gz = CreateLeaveOneOutPanelVCF.leave_out_panel_split_vcf_gz,
                        panel_split_vcf_gz_tbi = CreateLeaveOneOutPanelVCF.leave_out_panel_split_vcf_gz_tbi,
                        reference_fasta_fai = reference_fasta_fai,
                        chromosome = chromosomes[k],
                        genetic_map = genetic_maps[k],
                        output_prefix = leave_out_sample_name,
                        docker = kage_docker,
                        monitoring_script = monitoring_script
                }
            }

            call GLIMPSECaseGather {
                input:
                    chromosome_glimpse_vcf_gzs = GLIMPSECaseChromosome.chromosome_glimpse_vcf_gz,
                    chromosome_glimpse_vcf_gz_tbis = GLIMPSECaseChromosome.chromosome_glimpse_vcf_gz_tbi,
                    output_prefix = leave_out_sample_name,
                    docker = kage_docker,
                    monitoring_script = monitoring_script
            }

            # KAGE evaluation
            call CalculateMetrics as CalculateMetricsKAGE {
                input:
                    case_vcf_gz = KAGECase.kage_vcf_gz,
                    case_vcf_gz_tbi = KAGECase.kage_vcf_gz_tbi,
                    truth_vcf_gz = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz,
                    truth_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz_tbi,
                    chromosomes = chromosomes,
                    label = "KAGE",
                    sample_name = leave_out_sample_name,
                    docker = docker,
                    monitoring_script = monitoring_script
            }

            # KAGE+GLIMPSE evaluation
            call CalculateMetrics as CalculateMetricsKAGEPlusGLIMPSE {
                input:
                    case_vcf_gz = GLIMPSECaseGather.glimpse_vcf_gz,
                    case_vcf_gz_tbi = GLIMPSECaseGather.glimpse_vcf_gz_tbi,
                    truth_vcf_gz = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz,
                    truth_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz_tbi,
                    chromosomes = chromosomes,
                    label = "KAGE+GLIMPSE",
                    sample_name = leave_out_sample_name,
                    docker = docker,
                    monitoring_script = monitoring_script
            }

            if (do_pangenie) {
                # PanGenie case
                call PanGenieCase.PanGenie as PanGenieCase {
                    input:
                        panel_vcf_gz = CreateLeaveOneOutPanelVCF.leave_out_panel_vcf_gz,
                        panel_vcf_gz_tbi = CreateLeaveOneOutPanelVCF.leave_out_panel_vcf_gz_tbi,
                        input_fasta = PreprocessCaseReads.preprocessed_fasta,
                        reference_fasta = reference_fasta,
                        chromosomes = chromosomes,
                        sample_name = leave_out_sample_name,
                        output_prefix = leave_out_sample_name,
                        docker = pangenie_docker,
                        monitoring_script = monitoring_script
                }

                # PanGenie evaluation
                call CalculateMetrics as CalculateMetricsPanGenie {
                    input:
                        case_vcf_gz = PanGenieCase.genotyping_vcf_gz,
                        case_vcf_gz_tbi = PanGenieCase.genotyping_vcf_gz_tbi,
                        truth_vcf_gz = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz,
                        truth_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz_tbi,
                        chromosomes = chromosomes,
                        label = "PanGenie",
                        sample_name = leave_out_sample_name,
                        docker = docker,
                        monitoring_script = monitoring_script
                }
            }
        }
    }

    output {
    }
}

task PreprocessPanelVCF {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File repeat_mask_bed
        File segmental_duplications_bed
        File simple_repeats_bed
        File challenging_medically_relevant_genes_bed
        Array[String] chromosomes
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

        bcftools view --no-version ~{input_vcf_gz} -r ~{sep="," chromosomes} -Ou | \
            bcftools norm --no-version -m+ -N -Ou | \
            bcftools plugin fill-tags --no-version -Ou -- -t AF,AC,AN | \
            truvari anno svinfo | \
            bcftools annotate --no-version -a ~{repeat_mask_bed} -c CHROM,FROM,TO -m +RM -Ou | \
            bcftools annotate --no-version -a ~{segmental_duplications_bed} -c CHROM,FROM,TO -m +SD -Ou | \
            bcftools annotate --no-version -a ~{simple_repeats_bed} -c CHROM,FROM,TO -m +SR -Ou | \
            bcftools annotate --no-version -a ~{challenging_medically_relevant_genes_bed} -c CHROM,FROM,TO -m +CMRG -Oz -o ~{output_prefix}.preprocessed.vcf.gz
        bcftools index -t ~{output_prefix}.preprocessed.vcf.gz

        bcftools view --no-version --min-alleles 3  ~{output_prefix}.preprocessed.vcf.gz -Ou | \
            bcftools norm --no-version -m- -N -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.multi.split.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.multi.split.vcf.gz

        bcftools norm --no-version -m- -N ~{output_prefix}.preprocessed.vcf.gz -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.split.temp.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.split.temp.vcf.gz

        bcftools annotate --no-version -a ~{output_prefix}.preprocessed.multi.split.vcf.gz ~{output_prefix}.preprocessed.split.temp.vcf.gz -m +MULTIALLELIC \
            -Oz -o ~{output_prefix}.preprocessed.split.vcf.gz
        bcftools index -t ~{output_prefix}.preprocessed.split.vcf.gz
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
        File preprocessed_panel_vcf_gz = "~{output_prefix}.preprocessed.vcf.gz"
        File preprocessed_panel_vcf_gz_tbi = "~{output_prefix}.preprocessed.vcf.gz.tbi"
        File preprocessed_panel_split_vcf_gz = "~{output_prefix}.preprocessed.split.vcf.gz"
        File preprocessed_panel_split_vcf_gz_tbi = "~{output_prefix}.preprocessed.split.vcf.gz.tbi"
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
        set -e

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
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File cram_idx = "~{output_prefix}.cram.crai"
    }
}

task PreprocessCaseReads {
    input {
        File input_cram
        File input_cram_idx
        File reference_fasta
        File reference_fasta_fai
        Array[String] chromosomes
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    String filter_N_regex = "/^>/{N;/^>.*\\n.*N.*/d}"

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # hacky way to get chromosomes into bed file
        grep -P '~{sep="\\t|" chromosomes}\t' ~{reference_fasta_fai} | cut -f 1,2 | sed -e 's/\t/\t1\t/g' > chromosomes.bed

        # filter out read pairs containing N nucleotides
        # TODO move functionality into KAGE code
        samtools view --reference ~{reference_fasta} -@ $(nproc) -L chromosomes.bed -u ~{input_cram} | \
            samtools fasta --reference ~{reference_fasta} -@ $(nproc) | \
            sed -E '~{filter_N_regex}' > ~{output_prefix}.preprocessed.fasta
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 2])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File preprocessed_fasta = "~{output_prefix}.preprocessed.fasta"
    }
}

task CreateLeaveOneOutPanelVCF {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        String output_prefix
        Array[String] leave_out_sample_names

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

        bcftools view --no-version ~{input_vcf_gz} -s ^~{sep=',' leave_out_sample_names} --trim-alt-alleles -Ou | \
            bcftools view --no-version --min-alleles 2 -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LO.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.LO.vcf.gz

        bcftools norm --no-version -m- -N ~{output_prefix}.preprocessed.LO.vcf.gz -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LO.split.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.LO.split.vcf.gz

        # we need to drop multiallelics before LO and trimming, otherwise there may be representation issues in the graph
        bcftools view --no-version --min-alleles 2 --max-alleles 2  ~{input_vcf_gz} -Ou | \
            bcftools view --no-version -s ^~{sep=',' leave_out_sample_names} --trim-alt-alleles -Ou | \
            bcftools view --no-version --min-alleles 2 -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LO.bi.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.LO.bi.vcf.gz

         bcftools view --no-version --min-alleles 3  ~{input_vcf_gz} -Ou | \
            bcftools norm --no-version -m- -N -Ou | \
            bcftools view --no-version -s ^~{sep=',' leave_out_sample_names} --trim-alt-alleles -Ou | \
            bcftools view --no-version --min-alleles 2 -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LO.multi.split.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.LO.multi.split.vcf.gz
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
        File leave_out_panel_vcf_gz = "~{output_prefix}.preprocessed.LO.vcf.gz"
        File leave_out_panel_vcf_gz_tbi = "~{output_prefix}.preprocessed.LO.vcf.gz.tbi"
        File leave_out_panel_split_vcf_gz = "~{output_prefix}.preprocessed.LO.split.vcf.gz"
        File leave_out_panel_split_vcf_gz_tbi = "~{output_prefix}.preprocessed.LO.split.vcf.gz.tbi"
        File leave_out_panel_bi_vcf_gz = "~{output_prefix}.preprocessed.LO.bi.vcf.gz"
        File leave_out_panel_bi_vcf_gz_tbi = "~{output_prefix}.preprocessed.LO.bi.vcf.gz.tbi"
        File leave_out_panel_multi_split_vcf_gz = "~{output_prefix}.preprocessed.LO.multi.split.vcf.gz"
        File leave_out_panel_multi_split_vcf_gz_tbi = "~{output_prefix}.preprocessed.LO.multi.split.vcf.gz.tbi"
    }
}

task KAGECase {
    input {
        File input_fasta
        File panel_index
        File panel_kmer_index_only_variants_with_revcomp
        File panel_multi_split_vcf_gz # for filling in biallelic-only VCFs produced by KAGE
        File panel_multi_split_vcf_gz_tbi
        File reference_fasta_fai
        String output_prefix
        String sample_name
        Float average_coverage

        String docker
        File? monitoring_script

        String kmer_mapper_args = "-c 100000000"
        Boolean? ignore_helper_model = true
        String? kage_genotype_extra_args


        RuntimeAttributes runtime_attributes = {}
        Int? cpu = 8
    }

    Int cpu_resolved = select_first([runtime_attributes.cpu, cpu])

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        kmer_mapper map \
            ~{kmer_mapper_args} \
            -t ~{cpu_resolved} \
            -i ~{panel_kmer_index_only_variants_with_revcomp} \
            -f ~{input_fasta} \
            -o ~{output_prefix}.kmer_counts.npy

        kage genotype \
            -i ~{panel_index} \
            -c ~{output_prefix}.kmer_counts.npy \
            --average-coverage ~{average_coverage} \
            -s ~{sample_name} \
            ~{true='-I' false='' ignore_helper_model} \
            ~{kage_genotype_extra_args} \
            -o ~{output_prefix}.kage.bi.vcf

        # we need to add split multiallelics to biallelic-only KAGE VCF
        # create single-sample header from LO panel w/ split multiallelics
        bcftools view --no-version -h -G ~{panel_multi_split_vcf_gz} | \
            sed 's/##fileformat=VCFv4.2/##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods.">/g' | \
            sed 's/INFO$/INFO\tFORMAT\t~{sample_name}/g' > ~{output_prefix}.multi.split.header.txt
        # create single-sample missing genotypes from LO panel w/ split multiallelics
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
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File kmer_counts = "~{output_prefix}.kmer_counts.npy"
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
        File genetic_map
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
        set -e

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

task CalculateMetrics {
    input {
        File case_vcf_gz        # bi+multi split
        File case_vcf_gz_tbi
        File truth_vcf_gz       # bi+multi split
        File truth_vcf_gz_tbi
        Array[String] chromosomes
        String label
        String sample_name

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command <<<
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # split multiallelics in case (may be redundant)
        bcftools norm --no-version -r ~{sep="," chromosomes} -m- -N ~{case_vcf_gz} -Oz -o case.split.vcf.gz
        bcftools index -t case.split.vcf.gz

        # mark case variants in panel
        bcftools annotate --no-version -r ~{sep="," chromosomes} -a case.split.vcf.gz -m +CASE ~{truth_vcf_gz} -Oz -o panel.annot.vcf.gz
        bcftools index -t panel.annot.vcf.gz

        conda install -y seaborn numpy=1.26.4

        python - --case_vcf_gz case.split.vcf.gz \
                 --truth_vcf_gz panel.annot.vcf.gz \
                 --label ~{label} \
                 --sample_name ~{sample_name} \
                 <<-'EOF'
        import argparse
        import allel
        import numpy as np
        import pandas as pd
        import sklearn.metrics
        import matplotlib
        import matplotlib.pyplot as plt
        import seaborn as sns

        matplotlib.use('Agg')

        def load_callset(callset_file, **kwargs):
            return allel.read_vcf(callset_file, **kwargs)

        def get_gt_vsp(callset):
            return allel.GenotypeDaskArray(callset['calldata/GT'])

        def calculate_metrics_and_plot(case_vcf_gz, truth_vcf_gz, label, sample_name, average='macro'):
            samples = [sample_name]

            callset = load_callset(case_vcf_gz, samples=samples, fields='*', alt_number=1)
            gt_vp = get_gt_vsp(callset)[:, 0, :].compute()

            truth_callset = load_callset(truth_vcf_gz, samples=samples, fields='*', alt_number=1)
            is_case_V = truth_callset['variants/CASE']
            truth_gt_vp = get_gt_vsp(truth_callset)[is_case_V, 0, :].compute()

            is_any_missing_case_v = np.any(gt_vp == -1, axis=1)
            is_any_missing_truth_v = np.any(truth_gt_vp == -1, axis=1)
            is_missing_v = is_any_missing_case_v | is_any_missing_truth_v

            is_snp_v = truth_callset['variants/is_snp'][is_case_V]
            is_multiallelic_v = truth_callset['variants/MULTIALLELIC'][is_case_V]
            altfreq_v = truth_callset['variants/AF'][is_case_V]
            is_altfreq_v = [['[0%, 1%)', (0. <= altfreq_v) & (altfreq_v < 0.01)],
                            ['[1%, 5%)', (0.01 <= altfreq_v) & (altfreq_v < 0.05)],
                            ['[5%, 10%)', (0.05 <= altfreq_v) & (altfreq_v < 0.1)],
                            ['[10%, 50%)', (0.1 <= altfreq_v) & (altfreq_v < 0.50)],
                            ['[50%, 100%]', (0.5 <= altfreq_v) & (altfreq_v <= 1.)]]
            altlen_v = truth_callset['variants/altlen'][is_case_V]
            is_altlen_v = [['(-inf,-500]', altlen_v <= -500],
                           ['(-500,-50]', (-500 < altlen_v) & (altlen_v <= -50)],
                           ['(-50,-1]', (-50 < altlen_v) & (altlen_v <= -1)],
                           ['0 (SNP)', is_snp_v],
                           ['[1,50)', (1 <= altlen_v) & (altlen_v < 50)],
                           ['[50,500)', (50 <= altlen_v) & (altlen_v < 500)],
                           ['[500,5000)', (500 <= altlen_v) & (altlen_v < 5000)],
                           ['[5000,inf)', 5000 <= altlen_v]]
            is_sv_v = (altlen_v <= -50) | (altlen_v >= 50)
            num_i = len(is_altfreq_v)
            num_j = len(is_altlen_v)

            metrics_dicts = []

            for context in ['ALL', 'CMRG', 'US', 'RM', 'SD', 'SR']:

                if context == 'ALL':
                    is_context_v = True
                elif context == 'US':
                    is_context_v = ~(truth_callset[f'variants/RM'][is_case_V] |
                                     truth_callset[f'variants/SD'][is_case_V] |
                                     truth_callset[f'variants/SR'][is_case_V])
                else:
                    is_context_v = truth_callset[f'variants/{context}'][is_case_V]

                for allelic in ['biallelic+multiallelic', 'multiallelic', 'biallelic']:
                    num_evals = np.zeros((num_i, num_j))
                    precisions = np.zeros((num_i, num_j))
                    recalls = np.zeros((num_i, num_j))
                    f1s = np.zeros((num_i, num_j))
                    confusion_matrices = np.zeros((num_i, num_j, 3, 3))

                    if allelic == 'biallelic+multiallelic':
                        is_allelic_v = True
                    elif allelic == 'multiallelic':
                        is_allelic_v = is_multiallelic_v
                    else:
                        is_allelic_v = ~is_multiallelic_v

                    truth_missing_count = (is_any_missing_truth_v & is_allelic_v & is_context_v).sum()
                    case_missing_count = (is_any_missing_case_v & is_allelic_v & is_context_v).sum()
                    missing_count = (is_missing_v & is_allelic_v & is_context_v).sum()

                    for i, (filter_name_i, is_v_i) in enumerate(is_altfreq_v):
                        for j, (filter_name_j, is_v_j) in enumerate(is_altlen_v):
                            is_eval_v = ~is_missing_v & is_v_i & is_v_j & is_allelic_v & is_context_v

                            num_eval = np.sum(is_eval_v)
                            enc_gt_n = np.sum(gt_vp[is_eval_v], axis=1)
                            truth_enc_gt_n = np.sum(truth_gt_vp[is_eval_v], axis=1)
                            if num_eval == 0:
                                precision = np.nan
                                recall = np.nan
                                f1 = np.nan
                            else:
                                precision = sklearn.metrics.precision_score(truth_enc_gt_n, enc_gt_n, average=average)
                                recall = sklearn.metrics.recall_score(truth_enc_gt_n, enc_gt_n, average=average)
                                f1 = sklearn.metrics.f1_score(truth_enc_gt_n, enc_gt_n, average=average)
                            confusion_matrix = sklearn.metrics.confusion_matrix(truth_enc_gt_n, enc_gt_n, labels=[0, 1, 2])

                            num_evals[i][j] = num_eval
                            precisions[i][j] = precision
                            recalls[i][j] = recall
                            f1s[i][j] = f1
                            confusion_matrices[i][j] = confusion_matrix

                            metrics_dicts.append({
                                'LABEL': label,
                                'SAMPLE_NAME': sample_name,
                                'CONTEXT': context,
                                'ALLELIC': allelic,
                                'ALTFREQ': filter_name_i,
                                'ALTLEN': filter_name_j,
                                'NUM_EVAL': num_eval,
                                'PRECISION': precision,
                                'RECALL': recall,
                                'F1': f1,
                                'CONFUSION_MATRIX': confusion_matrix
                            })

                    if num_evals.sum() == 0:
                        continue

                    non_sv_enc_gt_n = np.sum(gt_vp[~is_missing_v & ~is_sv_v & is_allelic_v & is_context_v], axis=1)
                    non_sv_truth_enc_gt_n = np.sum(truth_gt_vp[~is_missing_v & ~is_sv_v & is_allelic_v & is_context_v], axis=1)
                    non_sv_precision = np.nan if non_sv_truth_enc_gt_n.size == 0 else sklearn.metrics.precision_score(non_sv_truth_enc_gt_n, non_sv_enc_gt_n, average=average)
                    non_sv_recall = np.nan if non_sv_truth_enc_gt_n.size == 0 else sklearn.metrics.recall_score(non_sv_truth_enc_gt_n, non_sv_enc_gt_n, average=average)
                    non_sv_f1 = np.nan if non_sv_truth_enc_gt_n.size == 0 else sklearn.metrics.f1_score(non_sv_truth_enc_gt_n, non_sv_enc_gt_n, average=average)
                    non_sv_count = np.sum(~is_sv_v & is_allelic_v & is_context_v)

                    sv_enc_gt_n = np.sum(gt_vp[~is_missing_v & is_sv_v & is_allelic_v & is_context_v], axis=1)
                    sv_truth_enc_gt_n = np.sum(truth_gt_vp[~is_missing_v & is_sv_v & is_allelic_v & is_context_v], axis=1)
                    sv_precision = sklearn.metrics.precision_score(sv_truth_enc_gt_n, sv_enc_gt_n, average=average)
                    sv_recall = sklearn.metrics.recall_score(sv_truth_enc_gt_n, sv_enc_gt_n, average=average)
                    sv_f1 = sklearn.metrics.f1_score(sv_truth_enc_gt_n, sv_enc_gt_n, average=average)
                    sv_count = np.sum(is_sv_v & is_allelic_v & is_context_v)

        #             fig, ax = plt.subplots(5, 1, figsize=(12, 20))
                    fig, ax = plt.subplots(4, 1, figsize=(12, 16))

                    ax[0] = sns.heatmap(num_evals, ax=ax[0], linewidths=1, linecolor='k', annot=True,
                                        norm=matplotlib.colors.LogNorm(), cmap='Blues')
                    cbar = ax[0].collections[0].colorbar
                    cbar.ax.set_ylabel('number of alt alleles', rotation=270, labelpad=20)
                    ax[0].set_title(f'{label}\n{sample_name}\ncontext = {context}, {allelic}\n\n' +
                                    f'alt allele count\n' +
                                    f'non-SV = {non_sv_count}, SV = {sv_count}\n' +
                                    f'(missing: truth = {truth_missing_count}, case = {case_missing_count}, union = {missing_count})')

        #             diag_mask = np.tile(np.eye(3), (num_i, num_j)).astype(bool)
        #             ax[1] = sns.heatmap(np.reshape(confusion_matrices, (3 * num_i, 3 * num_j)),
        #                                 ax=ax[1], linecolor='k',
        # #                                 norm=matplotlib.colors.LogNorm(),
        #                                 mask=~diag_mask, cmap='Greens')
        #             ax[1] = sns.heatmap(np.reshape(confusion_matrices, (3 * num_i, 3 * num_j)),
        #                                 ax=ax[1], linecolor='k',
        # #                                 norm=matplotlib.colors.LogNorm(),
        #                                 mask=diag_mask, cmap='Reds', cbar=False)
        #             cbar = ax[1].collections[1].colorbar
        # #             cbar.ax.set_ylabel(f'number of alt alleles', rotation=270, labelpad=40)
        #             ax[1].set_title(f'normalized confusion matrices')

                    ax[1] = sns.heatmap(precisions, ax=ax[1], linewidths=1, linecolor='k', annot=True,
                                        vmin=0.5, cmap='Greens')
                    cbar = ax[1].collections[0].colorbar
                    cbar.ax.set_ylabel(f'precision ({average})', rotation=270, labelpad=20)
                    ax[1].set_title(f'precision ({average})\n' +
                                    f'non-SV = {non_sv_precision:.4f}, SV = {sv_precision:.4f}')

                    ax[2] = sns.heatmap(recalls, ax=ax[2], linewidths=1, linecolor='k', annot=True,
                                        vmin=0.5, cmap='Greens')
                    cbar = ax[2].collections[0].colorbar
                    cbar.ax.set_ylabel(f'recall ({average})', rotation=270, labelpad=20)
                    ax[2].set_title(f'recall ({average})\n' +
                                    f'non-SV = {non_sv_recall:.4f}, SV = {sv_recall:.4f}')

                    ax[3] = sns.heatmap(f1s, ax=ax[3], linewidths=1, linecolor='k', annot=True,
                                        vmin=0.5, cmap='Greens')
                    cbar = ax[3].collections[0].colorbar
                    cbar.ax.set_ylabel(f'F1 ({average})', rotation=270, labelpad=20)
                    ax[3].set_title(f'F1 ({average})\n' +
                                    f'non-SV = {non_sv_f1:.4f}, SV = {sv_f1:.4f}')

                    for i, a in enumerate(ax):
                        a.set_ylabel('AF')
                        a.tick_params(bottom=False, left=False)
                        a.set_xticklabels([])
        #                 if i != 1:
                        a.set_yticklabels([filter_name for filter_name, _ in is_altfreq_v], rotation=0)
        #                 else:
        #                     a.set_yticklabels(sum([[None, filter_name, None]
        #                                            for filter_name, _ in is_altfreq_v], []),
        #                                       rotation=0)
                    ax[3].set_xlabel('ALT length - REF length (bp)')
                    ax[3].set_xticklabels([filter_name for filter_name, _ in is_altlen_v], rotation=0)

                    plt.savefig(f'{sample_name}.{context}.{allelic}.{label}.metrics.png')
                    plt.show()

            metrics_df = pd.DataFrame.from_dict(metrics_dicts)
            metrics_df.to_csv(f'{sample_name}.{label}.metrics.tsv', sep='\t', index=False)


        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--case_vcf_gz',
                                type=str)

            parser.add_argument('--truth_vcf_gz',
                                type=str)

            parser.add_argument('--label',
                                type=str)

            parser.add_argument('--sample_name',
                                type=str)

            args = parser.parse_args()

            calculate_metrics_and_plot(args.case_vcf_gz,
                                       args.truth_vcf_gz,
                                       args.label,
                                       args.sample_name)

        if __name__ == '__main__':
            main()
        EOF
    >>>

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
        File metrics_tsv = "~{sample_name}.~{label}.metrics.tsv"
        Array[File] metrics_plots = glob("~{sample_name}.*.png")
    }
}
