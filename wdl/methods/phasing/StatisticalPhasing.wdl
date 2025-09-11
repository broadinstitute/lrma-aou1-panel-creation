version 1.0

import "./Helper.wdl" as H


workflow StatisticalPhasing {

    input {

        File joint_short_vcf
        File? joint_short_vcf_tbi
        File? joint_sv_vcf
        File? joint_sv_vcf_tbi
        File reference_fasta
        File reference_fasta_fai
        File genetic_mapping_tsv_for_shapeit
        String chromosome
        String region
        String prefix
        String gcs_out_root_dir
        Int shapeit_num_threads
        String shapeit5_rare_extra_args
        Int merge_num_threads = 4

        Boolean shapeit5 = true

        Int shapeit_memory
        String shapeit5_common_extra_args 
    }

    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit)


    call H.SubsetVCF as SubsetVcfShort { input:
        vcf_gz = joint_short_vcf,
        locus = region
    }
    
    # if no SV VCF provided, use the short VCF as the input
    # concatenate the SV VCF if provided    
    if (defined(joint_sv_vcf) && defined(joint_sv_vcf_tbi)) {
        call H.SubsetVCF as SubsetVcfSV { input:
            vcf_gz = select_first([joint_sv_vcf ]),
            vcf_tbi = select_first([joint_sv_vcf_tbi]),
            locus = region
        }

        call H.FilterAndConcatVcfs { input:
            short_vcf = SubsetVcfShort.subset_vcf,
            short_vcf_tbi = SubsetVcfShort.subset_tbi,
            sv_vcf = SubsetVcfSV.subset_vcf,
            sv_vcf_tbi = SubsetVcfSV.subset_tbi,
            prefix = prefix + ".concat",
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            region = region,
        }
    } 

    call H.CreateChunks as CreateChunks { input:
        vcf = select_first([FilterAndConcatVcfs.filter_and_concat_vcf,SubsetVcfShort.subset_vcf]),
        tbi = select_first([FilterAndConcatVcfs.filter_and_concat_vcf_tbi,SubsetVcfShort.subset_tbi]),
        region = region,
        prefix = prefix + ".chunks",
        extra_chunk_args = "--thread 4 --window-size 50000 --buffer-size 5000"
    }

    Array[String] region_list = read_lines(CreateChunks.chunks)

    scatter (i in range(length(region_list))) {
        if (shapeit5) {
            call H.shapeit5_phase_common as Shapeit5_phase_common { input:
                vcf_input = select_first([FilterAndConcatVcfs.filter_and_concat_vcf,SubsetVcfShort.subset_vcf]),
                vcf_index = select_first([FilterAndConcatVcfs.filter_and_concat_vcf_tbi,SubsetVcfShort.subset_tbi]),
                mappingfile = genetic_mapping_dict[chromosome],
                region = region_list[i],
                prefix = prefix + ".filter_and_concat.phased",
                num_threads = shapeit_num_threads,
                memory = shapeit_memory,
                extra_args = shapeit5_common_extra_args
            }
        } 
        if (!shapeit5) {
            call H.Shapeit4 as Shapeit4 { input:
                vcf_input = select_first([FilterAndConcatVcfs.filter_and_concat_vcf,SubsetVcfShort.subset_vcf]),
                vcf_index = select_first([FilterAndConcatVcfs.filter_and_concat_vcf_tbi,SubsetVcfShort.subset_tbi]),
                mappingfile = genetic_mapping_dict[chromosome],
                region = region_list[i],
                prefix = prefix + ".filter_and_concat.phased",
                num_threads = shapeit_num_threads,
                memory = shapeit_memory
            }
        }
    }

    call H.LigateVcfs as LigateScaffold { input:
        vcfs = select_all(flatten([Shapeit5_phase_common.scaffold_vcf, Shapeit4.phased_bcf])),
        vcf_idxs = select_all(flatten([Shapeit5_phase_common.scaffold_vcf_index, Shapeit4.phased_bcf_index])),
        prefix = prefix + "." + region + ".scaffold.ligated"
    }

    # phase rare
    if (shapeit5) {
        scatter (i in range(length(region_list))) {
            call H.shapeit5_phase_rare as Shapeit5_phase_rare { input:
                vcf_input = select_first([FilterAndConcatVcfs.filter_and_concat_vcf, SubsetVcfShort.subset_vcf]),
                vcf_index = select_first([FilterAndConcatVcfs.filter_and_concat_vcf_tbi, SubsetVcfShort.subset_tbi]),
                scaffold_bcf = LigateScaffold.ligated_vcf_gz,
                scaffold_bcf_index = LigateScaffold.ligated_vcf_gz_tbi,
                mappingfile = genetic_mapping_dict[chromosome],
                chunk_region = region_list[i],
                scaffold_region = region,
                prefix = prefix + ".chunk.phase.rare.phased",
                chunknum = i,
                num_threads = shapeit_num_threads,
                memory = shapeit_memory,
                extra_args = shapeit5_rare_extra_args
            }
        }

        call H.BcftoolsConcatBCFs as ConcatRare { input:
            vcfs = Shapeit5_phase_rare.chunk_vcf,
            vcf_idxs = Shapeit5_phase_rare.chunk_vcf_index,
            prefix = prefix + ".phase.rare.concat"
        }

    }



    output {
        # File phased_scaffold_vcf = LigateScaffold.ligated_vcf_gz
        # File phased_scaffold_vcf_tbi = LigateScaffold.ligated_vcf_gz_tbi
        File phased_vcf = select_first([ConcatRare.concated_bcf,LigateScaffold.ligated_vcf_gz])
        File phasedvcf_tbi = select_first([ConcatRare.concated_bcf_index,LigateScaffold.ligated_vcf_gz_tbi])
    }
}
