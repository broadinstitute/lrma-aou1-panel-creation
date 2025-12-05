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
        String prefix
        Int hiphase_memory
        String hiphase_extra_args
        String sample_id
        String region
    }


    call H.SubsetVCF as SubsetVcfShort { input:
        vcf_gz = small_vcfs_directory + "/" + sample_id + '.vcf.gz',
        locus = region
    }

    call H.SubsetVCF as SubsetVcfTRGT { input:
        vcf_gz = trgt_vcf,
        locus = region
    }


    call H.HiphaseShortTrgt as Hiphase { input:
        bam = all_chr_bam,
        bai = all_chr_bai,
        unphased_snp_vcf = SubsetVcfShort.subset_vcf,
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
        File hiphase_trgt_vcf = Hiphase.phased_trgt_vcf
        File hiphase_trgt_tbi = Hiphase.phased_trgt_vcf_tbi
    }
}
