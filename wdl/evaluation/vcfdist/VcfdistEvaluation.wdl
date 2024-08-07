version 1.0

workflow VcfdistEvaluation {
    input {
        Array[String] samples
        File truth_vcf
        File eval_vcf
        File bed_file
        String region
        File reference_fasta
        File reference_fasta_fai
        String docker
    }

    scatter (sample in samples) {

        call SubsetSampleFromVcf as SubsetSampleFromVcfEval { input:
            vcf = eval_vcf,
            sample = sample,
            region = region,
            reference_fasta_fai = reference_fasta_fai
        }

        call SubsetSampleFromVcf as SubsetSampleFromVcfTruth { input:
            vcf = truth_vcf,
            sample = sample,
            region = region,
            reference_fasta_fai = reference_fasta_fai
        }

        call Vcfdist { input:
            sample = sample,
            eval_vcf = SubsetSampleFromVcfEval.single_sample_vcf,
            truth_vcf = SubsetSampleFromVcfTruth.single_sample_vcf,
            bed_file = bed_file,
            reference_fasta = reference_fasta
        }
    }

    output {
        Array[File] precision_recall_summary_tsvs = Vcfdist.precision_recall_summary_tsv
        Array[File] summary_vcfs = Vcfdist.summary_vcf
        Array[File] precision_recall_tsvs = Vcfdist.precision_recall_tsv
        Array[File] query_tsvs = Vcfdist.query_tsv
        Array[File] truth_tsvs = Vcfdist.truth_tsv
        Array[File] phasing_summary_tsvs = Vcfdist.phasing_summary_tsv
        Array[File] switch_flips_tsvs = Vcfdist.switch_flips_tsv
    }
}

task SubsetSampleFromVcf {
    input {
        File vcf
        String sample
        String region
        File reference_fasta_fai
    }

    Int disk_size = 1 + ceil(2 * (size(vcf, "GiB")))

    command <<<
        set -euxo pipefail

        bcftools index ~{vcf}
        bcftools view ~{vcf} \
            -s ~{sample} \
            -r ~{region} \
            -Oz -o ~{sample}.subset.g.vcf.gz
        bcftools reheader ~{sample}.subset.g.vcf.gz \
            --fai ~{reference_fasta_fai} \
            -Oz -o ~{sample}.subset.reheadered.g.vcf.gz
        bcftools index -t ~{sample}.subset.reheadered.g.vcf.gz
    >>>
    
    output {
        File single_sample_vcf = "~{sample}.subset.reheadered.g.vcf.gz"
        File single_sample_vcf_tbi = "~{sample}.subset.reheadered.g.vcf.gz.tbi"
    }

    runtime {
        cpu: 1
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task Vcfdist {
    input {
        String sample
        File eval_vcf
        File truth_vcf
        File bed_file
        File reference_fasta
        Int verbosity = 1

        String docker = "timd1/vcfdist:v2.5.3"
        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 16
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        vcfdist \
            ~{eval_vcf} \
            ~{truth_vcf} \
            ~{reference_fasta} \
            -b ~{bed_file} \
            -v ~{verbosity}

        for tsv in $(ls *.tsv); do mv $tsv ~{sample}.$tsv; done
    >>>

    output {
        File precision_recall_summary_tsv = "~{sample}.precision-recall-summary.tsv"
        File summary_vcf = "~{sample}.summary.vcf"
        File precision_recall_tsv = "~{sample}.precision-recall.tsv"
        File query_tsv = "~{sample}.query.tsv"
        File truth_tsv = "~{sample}.truth.tsv"
        File phasing_summary_tsv = "~{sample}.phasing-summary.tsv"
        File switch_flips_tsv = "~{sample}.switchflips.tsv"
    }

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}
