version 1.0

import "./Helper.wdl" as H


workflow LigateVcf {

    input {

        Array[File] phased_vcfs
        Array[File]? phased_vcf_tbis
        File reference_fasta
        File reference_fasta_fai
        String prefix
        String chromosome
        String gcs_out_root_dir
    }

    call H.LigateVcfs as LigateScaffold { input:
        vcfs = phased_vcfs,
        prefix = prefix + "." + chromosome + ".shapeit4.phased.ligated"
    }

    output {
        File phased_vcf = LigateScaffold.ligated_vcf_gz
        File phased_vcf_tbi = LigateScaffold.ligated_vcf_gz_tbi
    }
}
