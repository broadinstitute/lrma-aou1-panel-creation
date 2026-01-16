version 1.0


workflow StatisticalPhasing {

    input {

        File joint_short_vcf
        File joint_short_vcf_idx
        File joint_sv_vcf
        File joint_sv_vcf_idx
        File reference_fasta
        File reference_fasta_fai
        File genetic_mapping_tsv_for_shapeit
        String chromosome
        String region
        String output_prefix

        Int filter_and_concat_shard_size = 1000000
        String filter_and_concat_short_filter_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))<50'"
        String filter_and_concat_sv_filter_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))>=50'"

        # inputs for FixVariantCollisions
        File fix_variant_collisions_java
        Int? operation
        String? weight_tag
        Int? is_weight_format_field

        String extra_chunk_args = "--thread $(nproc) --window-size 2000000 --buffer-size 200000"

        Boolean do_shapeit5 = true
        String shapeit4_common_extra_args = "--thread $(nproc)"
        String shapeit5_rare_extra_args =  "--thread $(nproc)"
        String filter_common_args = "-i 'MAF>=0.001'"
        String filter_rare_args = "-i 'MAF<0.001'"
    }

    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit)


    call CreateShards { input:
        region = region,
        bin_size = filter_and_concat_shard_size,
        pad_size = 0,
        output_prefix = output_prefix + ".subset_create_chunks"
    }

    scatter (s in range(length(CreateShards.shard_regions))) {
        String shard_region = CreateShards.shard_regions[s]
        
        call SubsetVCFStreaming as SubsetVcfShort { input:
            vcf = joint_short_vcf,
            vcf_idx = joint_short_vcf_idx,
            region = shard_region,
            output_prefix = output_prefix + ".subsetShort.shard-" + s
        }

        call SubsetVCF as SubsetVcfSV { input:
            vcf = select_first([joint_sv_vcf ]),
            vcf_idx = select_first([joint_sv_vcf_idx]),
            region = shard_region,
            output_prefix = output_prefix + ".subsetSV.shard-" + s
        }

        call FilterAndConcatVcfs { input:
            short_vcf = SubsetVcfShort.subset_vcf,
            short_vcf_idx = SubsetVcfShort.subset_idx,
            sv_vcf = SubsetVcfSV.subset_vcf,
            sv_vcf_idx = SubsetVcfSV.subset_idx,
            output_prefix = output_prefix + ".filterAndConcat.shard-" + s,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            region = shard_region,
            filter_and_concat_short_filter_args = filter_and_concat_short_filter_args,
            filter_and_concat_sv_filter_args = filter_and_concat_sv_filter_args
        }
        
        # added variant collision fix step
        call FixVariantCollisions { input:
            phased_bcf = select_first([FilterAndConcatVcfs.filter_and_concat_vcf, SubsetVcfShort.subset_vcf]),
            fix_variant_collisions_java = fix_variant_collisions_java,
            operation = operation,
            weight_tag = weight_tag,
            is_weight_format_field = is_weight_format_field,
            output_prefix = output_prefix + ".collisionless.shard-" + s,
        }
    }

    call BcftoolsConcatNaive { input:
        vcfs = FixVariantCollisions.phased_collisionless_vcf,
        vcf_idxs = FixVariantCollisions.phased_collisionless_vcf_idx,
        output_prefix = output_prefix + ".collisionless"
    }

    call CreateShapeitChunks { input:
        vcf = BcftoolsConcatNaive.concatenated_vcf,
        vcf_idx = BcftoolsConcatNaive.concatenated_vcf_idx,
        region = region,
        extra_chunk_args = extra_chunk_args
    }

#    Array[String] region_list = read_lines(CreateShapeitChunks.common_chunks) # for shapeit4

#    scatter (i in range(length(region_list))) {
#        # phase common using shapeit4
#        if (!shapeit5) {
#            call Shapeit4 as Shapeit4_all { input:
#                vcf_input = ConcatSubsets.concatenated_vcf,
#                vcf_index = ConcatSubsets.concatenated_vcf_tbi,
#                mappingfile = genetic_mapping_dict[chromosome],
#                region = region_list[i],
#                output_prefix = output_prefix + ".filter_and_concat.phased",
#                cpu = shapeit4_cpu,
#                memory = shapeit4_memory,
#                extra_args = shapeit4_common_extra_args
#            }
#        }
#        if (shapeit5) {
#            call FilterCommonandRareVariants { input:
#                vcf_gz = ConcatSubsets.concatenated_vcf,
#                vcf_gz_tbi = ConcatSubsets.concatenated_vcf_tbi,
#                output_prefix = output_prefix + ".filter_common_and_rare",
#                region = region_list[i],
#                filter_common_args = filter_common_args,
#                filter_rare_args = filter_rare_args
#            }
#            call Shapeit4 as Shapeit4_common { input:
#                vcf_input = FilterCommonandRareVariants.filter_common_vcf,
#                vcf_index = FilterCommonandRareVariants.filter_common_vcf_tbi,
#                mappingfile = genetic_mapping_dict[chromosome],
#                region = region_list[i],
#                output_prefix = output_prefix + ".filter_and_concat.common",
#                cpu = shapeit4_cpu,
#                memory = shapeit4_memory,
#                extra_args = shapeit4_common_extra_args
#            }
#        }
#        
#    }

#    call LigateVcfs as LigateScaffold { input:
#        vcfs = select_all(flatten([Shapeit4_common.phased_bcf,Shapeit4_all.phased_bcf])),
#        vcf_idxs = select_all(flatten([Shapeit4_common.phased_bcf_index, Shapeit4_all.phased_bcf_index])),
#        output_prefix = output_prefix + ".scaffold.ligated"
#    }

#    # phase rare
#    if (shapeit5) {
#        Array[String] rare_region_list = read_lines(CreateShapeitChunks.rare_chunks) # for shapeit5
#        scatter (i in range(length(region_list))) {
#            call Shapeit5PhaseRareNew as Shapeit5_phase_rare { input:
#                vcf_input = ConcatSubsets.concatenated_vcf,
#                vcf_index = ConcatSubsets.concatenated_vcf_tbi,
#                scaffold_bcf = LigateScaffold.ligated_vcf_gz,
#                scaffold_bcf_index = LigateScaffold.ligated_vcf_gz_tbi,
#                mappingfile = genetic_mapping_dict[chromosome],
#                chunk_region = rare_region_list[i],
#                scaffold_region = region_list[i],
#                output_prefix = output_prefix + ".chunk.phase.rare.phased",
#                chunknum = i,
#                cpu = shapeit5_cpu,
#                memory = shapeit5_memory,
#                extra_args = shapeit5_rare_extra_args
#            }
#        }

#        # concat rare, following the tutorial
#        call BcftoolsConcatNaive as ConcatRare { input:
#            vcfs = flatten([Shapeit5_phase_rare.chunk_vcf]),
#            vcf_tbis = flatten([Shapeit5_phase_rare.chunk_vcf_index]),
#            output_prefix = output_prefix + ".phase.common.rare.concat"
#        }

#    }

#    output {
#        File phased_vcf = select_first([ConcatRare.concatenated_vcf,LigateScaffold.ligated_vcf_gz])
#        File phased_vcf_tbi = select_first([ConcatRare.concatenated_vcf_tbi,LigateScaffold.ligated_vcf_gz_tbi])
#    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Boolean? use_ssd
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

task CreateShards {
    input {
        String region
        Int bin_size
        Int pad_size
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        python - --region ~{region} \
                 --bin_size ~{bin_size} \
                 --pad_size ~{pad_size} \
                 --output_file ~{output_prefix} \
                 <<-'EOF'
        import argparse

        def split_region(region):
            chromosome, span = region.split(":")
            start, end = span.split("-")
            return(chromosome, int(start), int(end))

        def split_region_to_intervals(region, bin_size, pad_size):
            chromo, start, end = split_region(region)
            bin_num = (end - start)//bin_size
            intervals = [(chromo, start, start + bin_size + pad_size)]
            for i in range(1, bin_num):
                start_pos = start + i*bin_size - pad_size
                end_pos = start_pos + bin_size + pad_size
                intervals.append((chromo, start_pos, end_pos))
            if end > start + bin_num*bin_size:
                intervals.append((chromo, start + bin_num*bin_size - pad_size, end))
            return(intervals)

        def write_output_file(content, output_file):
            with open(output_file, "w") as f:
                for item in content:
                    l = "%s:%d-%d" % (item[0], item[1], item[2])
                    f.write(l+ "\n")

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--region',
                                type=str)

            parser.add_argument('--output_file',
                                type=str)

            parser.add_argument('--bin_size',
                    type=int)

            parser.add_argument('--pad_size',
                    type=int)

            args = parser.parse_args()

            intervals = split_region_to_intervals(args.region, args.bin_size, args.pad_size)
            write_output_file(intervals, args.output_file + ".txt")

        if __name__ == "__main__":
            main()
        EOF
    >>>

    output {
        Array[String] shard_regions = read_lines("~{output_prefix}.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            10,
        boot_disk_gb:       10,
        use_ssd:            false,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task SubsetVCF {
    input {
        File vcf
        File vcf_idx
        String region
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * ceil(size([vcf, vcf_idx], "GB"))

    command <<<
        set -euxo pipefail

        bcftools view --no-version ~{vcf} --regions ~{region} -Ob -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File subset_vcf = "~{output_prefix}.bcf"
        File subset_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            false,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SubsetVCFStreaming {
    parameter_meta {
        vcf: {
            localization_optional: true
        }
        vcf_idx: {
            localization_optional: true
        }
    }

    input {
        File vcf
        File vcf_idx
        String region
        String output_prefix

        Int disk_gb = 100

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        bcftools view --no-version ~{vcf} --regions ~{region} -Ob -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File subset_vcf = "~{output_prefix}.bcf"
        File subset_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            false,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FilterAndConcatVcfs {
    input {
        File short_vcf         # multiallelic
        File short_vcf_idx
        File sv_vcf            # biallelic
        File sv_vcf_idx
        String output_prefix
        String region
        File reference_fasta
        File reference_fasta_fai
        String filter_and_concat_short_filter_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))<50'"
        String filter_and_concat_sv_filter_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))>=50'"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * (ceil(size(short_vcf, "GB")) + size(sv_vcf, "GB"))

    command <<<
        set -euxo pipefail

        # filter SV
        bcftools +fill-tags --no-version -r ~{region} ~{sv_vcf} -- -t AF,AC,AN | \
            bcftools view --no-version ~{filter_and_concat_sv_filter_args} \
                -Ob -o ~{output_prefix}.SV.bcf
        bcftools index ~{output_prefix}.SV.bcf

        # split to biallelic and filter short
        bcftools norm --no-version -r ~{region} -m-any -N -f ~{reference_fasta} ~{short_vcf} | \
            bcftools +fill-tags --no-version -- -t AF,AC,AN | \
            bcftools view --no-version ~{filter_and_concat_short_filter_args} | \
            bcftools sort --no-version -Ob -o ~{output_prefix}.short.bcf
        bcftools index ~{output_prefix}.short.bcf

        # concatenate with deduplication; providing SV VCF as first argument preferentially keeps those records
        bcftools concat --no-version \
            ~{output_prefix}.SV.bcf \
            ~{output_prefix}.short.bcf \
            --allow-overlaps --remove-duplicates | \
            bcftools sort --no-version -Ob -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File filter_and_concat_vcf = "~{output_prefix}.bcf"
        File filter_and_concat_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FixVariantCollisions {
    input {
        File phased_bcf                     # biallelic
        File fix_variant_collisions_java
        Int operation = 1                   # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "UNIT_WEIGHT"   # ID of the weight field; if this field is not found, all weights are set to one; weights are assumed to be non-negative
        Int is_weight_format_field = 0      # given a VCF record in a sample, assign it a weight encoded in the sample column (1) or in the INFO field (0)
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * (ceil(size(phased_bcf, "GB")))

    command <<<
        set -euxo pipefail

        # convert bcf to vcf.gz
        bcftools view --no-version ~{phased_bcf} -Oz -o phased.vcf.gz
        bcftools index -t phased.vcf.gz

        java ~{fix_variant_collisions_java} \
            phased.vcf.gz \
            ~{operation} \
            ~{weight_tag} \
            ~{is_weight_format_field} \
            collisionless.vcf \
            windows.txt \
            histogram.txt \
            null                            # do not output figures

        # replace all missing alleles (correctly) emitted with reference alleles, since this is expected by PanGenie panel-creation script
        bcftools view --no-version collisionless.vcf | \
            sed -e 's/\.|0/0|0/g' | sed -e 's/0|\./0|0/g' | sed -e 's/\.|1/0|1/g' | sed -e 's/1|\./1|0/g' | sed -e 's/\.|\./0|0/g' | \
            bcftools view --no-version -Oz -o ~{output_prefix}.phased.collisionless.vcf.gz
        # index and convert via vcf.gz to avoid errors from missing header lines
        bcftools index -t ~{output_prefix}.phased.collisionless.vcf.gz
        bcftools view --no-version ~{output_prefix}.phased.collisionless.vcf.gz -Ob -o ~{output_prefix}.phased.collisionless.bcf
        bcftools index ~{output_prefix}.phased.collisionless.bcf
    >>>

    output {
        File phased_collisionless_vcf = "~{output_prefix}.phased.collisionless.bcf"
        File phased_collisionless_vcf_idx = "~{output_prefix}.phased.collisionless.bcf.csi"
        File windows = "windows.txt"
        File histogram = "histogram.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CreateShapeitChunks {
    input {
        File vcf
        File vcf_idx
        String region
        String extra_chunk_args = "--thread $(nproc) --window-size 5000000 --buffer-size 500000"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 2 * ceil(size([vcf, vcf_idx], "GB")) + 1

    command <<<
        set -euxo pipefail

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_chunk_static
        chmod +x GLIMPSE_chunk_static

        ./GLIMPSE_chunk_static \
            -I ~{vcf} \
            --region ~{region} \
            ~{extra_chunk_args} \
            -O chunks.txt

        # cut chunks + buffers
        cut -f 3 chunks.txt > common.chunks.regions.txt
        cut -f 4 chunks.txt > rare.chunks.regions.txt
    >>>

    output {
        File common_chunks = "common.chunks.regions.txt"
        File rare_chunks = "rare.chunks.regions.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task BcftoolsConcatNaive {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * ceil(size([vcfs], "GB"))

    command <<<
        set -euxo pipefail

        bcftools concat --no-version ~{sep=" " vcfs} --naive-force | \
            bcftools sort --no-version -Oz -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File concatenated_vcf = "~{output_prefix}.bcf"
        File concatenated_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

#task Shapeit4 {
#    input{
#        File vcf_input
#        File vcf_index
#        File mappingfile
#        String region
#        String output_prefix
#        Int cpu
#        Int memory
#        String extra_args = "--thread $(nproc)"

#        RuntimeAttr? runtime_attr_override
#        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
#    }
#    command <<<
#        set -euxo pipefail

#        shapeit4.2 --input ~{vcf_input} \
#                --map ~{mappingfile} \
#                --region ~{region} \
#                --sequencing \
#                --output ~{output_prefix}.bcf \
#                ~{extra_args}
#        bcftools index ~{output_prefix}.bcf
#    >>>

#    output{
#        # File resouce_monitor_log = "resources.log"
#        File phased_bcf = "~{output_prefix}.bcf"
#        File phased_bcf_index = "~{output_prefix}.bcf.csi"
#    }

# #########################
#    RuntimeAttr default_attr = object {
#        cpu_cores:          cpu,
#        mem_gb:             memory,
#        disk_gb:            100,
#        boot_disk_gb:       100,
#        preemptible_tries:  0,
#        max_retries:        0,
#        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/shapeit4:v1"
#    }
#    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
#    runtime {
#        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
#        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
#        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
#        zones: zones
#        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
#        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
#        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
#        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
#    }
#}

#task LigateVcfs {

#    input {
#        Array[File] vcfs
#        Array[File]? vcf_idxs
#        String output_prefix

#        RuntimeAttr? runtime_attr_override
#    }

#    Int disk_gb = 2*ceil(size(vcfs, "GB")) + 100

#    command <<<
#        set -euxo pipefail
#        if ! ~{defined(vcf_idxs)}; then
#            for ff in ~{sep=' ' vcfs}; do bcftools index $ff; done
#        fi

#        ligate_static --input ~{write_lines(vcfs)} --output ~{output_prefix}.vcf.gz
#        bcftools +fill-tags --no-version ~{output_prefix}.vcf.gz -- -t AF,AC,AN | \
#            bcftools annotate --no-version  -x 'FORMAT/PP' \
#                -Oz -o ~{output_prefix}.tagged.vcf.gz
#        
#        bcftools index -t ~{output_prefix}.tagged.vcf.gz
#    >>>

#    output {
#        File ligated_vcf_gz = "~{output_prefix}.tagged.vcf.gz"
#        File ligated_vcf_gz_tbi = "~{output_prefix}.tagged.vcf.gz.tbi"
#    }

#    #########################
#    RuntimeAttr default_attr = object {
#        cpu_cores:          2,
#        mem_gb:             8,
#        disk_gb:            disk_gb,
#        boot_disk_gb:       10,
#        preemptible_tries:  2,
#        max_retries:        1,
#        docker:"hangsuunc/shapeit5:v1"
#    }
#    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
#    runtime {
#        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
#        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
#        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
#        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
#        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
#        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
#        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
#    }
#}

#task Shapeit5PhaseRareNew{
#    input{
#        File vcf_input
#        File vcf_index
#        File scaffold_bcf
#        File scaffold_bcf_index
#        File mappingfile
#        String chunk_region
#        String scaffold_region
#        String output_prefix
#        Int chunknum
#        Int cpu
#        Int memory
#        String extra_args = "--thread $(nproc)"
#        
#        RuntimeAttr? runtime_attr_override
#        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
#    }
#    command <<<
#        set -euxo pipefail

#        bcftools +fill-tags ~{scaffold_bcf} --no-version -Ob -o tmp.scaffold.out.bcf -- -t AN,AC
#        bcftools index tmp.scaffold.out.bcf

#        bcftools +fill-tags ~{vcf_input} --no-version -Ob -o tmp.out.bcf -- -t AN,AC
#        bcftools index tmp.out.bcf
#        
#        # try to fix bugs in https://github.com/odelaneau/shapeit5/issues/33
#        # replace filtering with setGT to set missing genotypes to 0|0
#        bcftools +setGT tmp.out.bcf --no-version -Ob -o tmp.rare.out.bcf -- -t . -n 0p  
#        bcftools index tmp.rare.out.bcf
#        
#        /shapeit5/phase_rare --input tmp.rare.out.bcf \
#                    --scaffold tmp.scaffold.out.bcf \
#                    --map ~{mappingfile} \
#                    --input-region ~{chunk_region} \
#                    --scaffold-region ~{scaffold_region} \
#                    --output ~{output_prefix}.chunk.~{chunknum}.bcf \
#                    ~{extra_args}

#        bcftools +fill-tags ~{output_prefix}.chunk.~{chunknum}.bcf --no-version -Ob -o ~{output_prefix}.chunk.~{chunknum}.tagged.bcf -- -t AF,AC,AN
#        bcftools index ~{output_prefix}.chunk.~{chunknum}.tagged.bcf

#    >>>

#    output{
#        File chunk_vcf = "~{output_prefix}.chunk.~{chunknum}.tagged.bcf"
#        File chunk_vcf_index = "~{output_prefix}.chunk.~{chunknum}.tagged.bcf.csi"
#    }

#    Int disk_gb = 100 + ceil(2 * size(vcf_input, "GiB"))

#    #########################
#    RuntimeAttr default_attr = object {
#        cpu_cores:          cpu,
#        mem_gb:             memory,
#        disk_gb:            disk_gb,
#        boot_disk_gb:       100,
#        preemptible_tries:  0,
#        max_retries:        0,
#        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/shapeit5:develop"
#    }
#    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
#    runtime {
#        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
#        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
#        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
#        zones: zones
#        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
#        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
#        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
#        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
#    }
#}


#task FilterCommonandRareVariants {

#    input {
#        File vcf_gz        
#        File vcf_gz_tbi
#        String output_prefix
#        String region
#        String filter_common_args 
#        String filter_rare_args

#        RuntimeAttr? runtime_attr_override
#    }

#    command {
#        set -euxo pipefail

#        # filter common
#        bcftools +fill-tags --no-version -r ~{region} ~{vcf_gz} -- -t AF,AC,AN | \
#            bcftools view ~{filter_common_args} \
#                -Oz -o ~{output_prefix}.common.vcf.gz
#        bcftools index -t ~{output_prefix}.common.vcf.gz

#        # filter rare
#        bcftools +fill-tags --no-version -r ~{region} ~{vcf_gz} -- -t AF,AC,AN | \
#            bcftools view ~{filter_rare_args} \
#                --no-version -Oz -o ~{output_prefix}.rare.vcf.gz
#        bcftools index -t ~{output_prefix}.rare.vcf.gz

#    }

#    output {
#        File filter_common_vcf = "~{output_prefix}.common.vcf.gz"
#        File filter_common_vcf_tbi = "~{output_prefix}.common.vcf.gz.tbi"
#        File filter_rare_vcf = "~{output_prefix}.rare.vcf.gz"
#        File filter_rare_vcf_tbi = "~{output_prefix}.rare.vcf.gz.tbi"
#    }

#    #########################
#    RuntimeAttr default_attr = object {
#        cpu_cores:          2,
#        mem_gb:             8,
#        disk_gb:            1000,
#        boot_disk_gb:       10,
#        preemptible_tries:  0,
#        max_retries:        1,
#        docker:"hangsuunc/shapeit5:v1"
#    }
#    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
#    runtime {
#        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
#        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
#        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
#        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
#        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
#        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
#        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
#    }
#}
