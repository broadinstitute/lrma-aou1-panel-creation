version 1.0

import "./Helper.wdl" as H


workflow PhysicalAndStatisticalPhasingTRGT {

    input {
        File all_chr_bam
        File all_chr_bai
        File reference_fasta
        File reference_fasta_fai
        File trgt_vcf
        String small_vcfs_directory
        String sv_vcfs_directory
        String prefix
        Int hiphase_memory
        String hiphase_extra_args
        String sample_id
        String region
    }

    call ConvertLowerCase {input:
        vcf = sv_vcfs_directory + "/" + sample_id + '.vcf.gz',
        prefix = sample_id + ".uppercased_sv_cleaned"
            
    }

    call UnphaseGenotypes as UnphaseSVGenotypes { input:
        vcf = ConvertLowerCase.subset_vcf,
        vcf_tbi = ConvertLowerCase.subset_tbi,
        prefix = prefix + ".unphased"
    }

    call H.SubsetVCF as SubsetVcfShort { input:
        vcf_gz = small_vcfs_directory + "/" + sample_id + '.vcf.gz',
        locus = region
    }

    call H.SubsetVCF as SubsetVcfSV { input:
        vcf_gz = UnphaseSVGenotypes.unphased_vcf,
        locus = region
    }

    call H.SubsetVCF as SubsetVcfTRGT { input:
        vcf_gz = trgt_vcf,
        locus = region
    }


    call H.HiphaseSVTrgt as Hiphase { input:
        bam = all_chr_bam,
        bai = all_chr_bai,
        unphased_snp_vcf = SubsetVcfShort.subset_vcf,
        unphased_sv_vcf = SubsetVcfSV.subset_vcf,
        unphased_trgt_vcf = SubsetVcfTRGT.subset_vcf,
        ref_fasta = reference_fasta,
        ref_fasta_fai = reference_fasta_fai,
        samplename = sample_id,
        memory = hiphase_memory,
        extra_args = hiphase_extra_args
    }

    output {
        File hiphase_short_vcf = Hiphase.phased_snp_vcf
        File hiphase_short_tbi = Hiphase.phased_snp_vcf_tbi
        File hiphase_sv_vcf = Hiphase.phased_sv_vcf
        File hiphase_sv_tbi = Hiphase.phased_sv_vcf_tbi
        File hiphase_trgt_vcf = Hiphase.phased_trgt_vcf
        File hiphase_trgt_tbi = Hiphase.phased_trgt_vcf_tbi
    }
}

task ConvertLowerCase {
    input {
        File vcf
        String prefix
    }
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cp ~{docker_dir}/convert_lower_case.py ~{work_dir}/convert_lower_case.py
        cd ~{work_dir}

        python convert_lower_case.py -i ~{vcf} -o ~{prefix}.vcf
        bgzip ~{prefix}.vcf ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{work_dir}/~{prefix}.vcf.gz"
        File subset_tbi = "~{work_dir}/~{prefix}.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 2
        memory:  "8 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"hangsuunc/cleanvcf:v1"
    }
}

task UnphaseGenotypes {

    input {
        File vcf
        File vcf_tbi
        String prefix
    }

    command <<<
        set -euxo pipefail

        # set (a)ll genotypes to (u)nphased and sort by allele (e.g., 1|0 becomes 0/1)
        bcftools +setGT ~{vcf} -Oz -o ~{prefix}.vcf.gz -- --target-gt a --new-gt u
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File unphased_vcf = "~{prefix}.vcf.gz"
        File unphased_vcf_tbi = "~{prefix}.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}
