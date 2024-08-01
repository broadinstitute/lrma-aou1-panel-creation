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
        String docker = "timd1/vcfdist:v2.5.3"
        Int verbosity = 1
    }

    scatter (sample_id in samples)  {

        call SplitVCFbySample as SP_eval { input:
            joint_vcf = eval_vcf,
            reference_index = reference_fasta_fai,
            region = region,
            samplename = sample_id
        }

        call SplitVCFbySample as SP_truth { input:
            joint_vcf = truth_vcf,
            reference_index = reference_fasta_fai,
            region = region,
            samplename = sample_id,
        }

        call RunVcfdistTask { input:
            sampleid = sample_id,
            truth_vcf = SP_truth.single_sample_vcf,
            eval_vcf = SP_eval.single_sample_vcf,
            bed_file = bed_file,
            reference_fasta = reference_fasta,
            docker = docker,
            verbosity = verbosity
        }
    }

    output {
        Array[File] prs_tsv = RunVcfdistTask.prs_tsv
        Array[File] vcfdistsummary = RunVcfdistTask.summary
        Array[File] precrec_tsv = RunVcfdistTask.precrec_tsv
        Array[File] query_tsv = RunVcfdistTask.query_tsv
        Array[File] truth_tsv = RunVcfdistTask.truth_tsv
        
    }
}

task RunVcfdistTask {
    input {
        File truth_vcf
        File eval_vcf
        File bed_file
        File reference_fasta
        Int verbosity
        String sampleid
        
        String docker
        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 16
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        vcfdist \
            ~{eval_vcf} \
            ~{truth_vcf} \
            ~{reference_fasta} \
            -b ~{bed_file} \
            -v ~{verbosity}

        mv "precision-recall-summary.tsv" "~{sampleid}.precision-recall-summary.tsv"
        mv "phasing-summary.tsv" "~{sampleid}.phasing-summary.tsv"
        mv "summary.vcf" "~{sampleid}.summary.vcf"
        mv "precision-recall.tsv" "~{sampleid}.precision-recall.tsv"
        mv "query.tsv" "~{sampleid}.query.tsv"
        mv "truth.tsv" "~{sampleid}.truth.tsv"
        mv "switchflips.tsv" "~{sampleid}.switchflips.tsv"

    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File prs_tsv = "~{sampleid}.precision-recall-summary.tsv"
        File summary = "~{sampleid}.summary.vcf"
        File precrec_tsv = "~{sampleid}.precision-recall.tsv"
        File query_tsv = "~{sampleid}.query.tsv"
        File truth_tsv = "~{sampleid}.truth.tsv"
        File phasing_tsv = "~{sampleid}.phasing-summary.tsv"
        File switch_flip_tsv = "~{sampleid}.switchflips.tsv"
    }
}

task SplitVCFbySample {
    input{       
        File joint_vcf
        File reference_index
        String region
        String samplename
    }
    
    command <<<
        set -x pipefail

        bcftools index ~{joint_vcf}

        bcftools view -s ~{samplename} ~{joint_vcf} -r ~{region} -o ~{samplename}.subset.g.vcf.gz

        bcftools reheader --fai ~{reference_index} ~{samplename}.subset.g.vcf.gz > ~{samplename}.subset.reheadered.g.vcf.gz

        tabix -p vcf ~{samplename}.subset.reheadered.g.vcf.gz

    >>>
    
    output {
        File single_sample_vcf = "~{samplename}.subset.reheadered.g.vcf.gz"
        File single_sample_vcf_tbi = "~{samplename}.subset.reheadered.g.vcf.gz.tbi"
    }


    Int disk_size = 1 + ceil(2 * (size(joint_vcf, "GiB")))

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
