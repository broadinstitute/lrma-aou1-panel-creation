version 1.0


workflow StatisticalPhasing {

    input {

        File joint_short_vcf
        File joint_short_vcf_idx
        File joint_sv_vcf
        File joint_sv_vcf_idx
        File reference_fasta
        File reference_fasta_fai
        File genetic_maps_tsv
        String chromosome
        String region
        String output_prefix

        Int filter_and_concat_shard_size = 2000000
        String? filter_and_concat_short_filter_args
        String filter_and_concat_short_view_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))<50'"
        String filter_and_concat_sv_view_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))>=50'"

        # inputs for FixVariantCollisions (see documentation for arguments in task)
        File fix_variant_collisions_java
        Int operation = 1
        String weight_tag = "SCORE"
        Int is_weight_format_field = 0
        Float default_weight = 0.5

        String chunk_extra_args = "--thread $(nproc) --window-size 2000000 --buffer-size 200000"

        Boolean do_shapeit5 = true
        String shapeit4_extra_args = "--thread $(nproc) --use-PS 0.0001"
        String shapeit5_extra_args =  "--thread $(nproc)"
        String filter_common_args = "-i 'MAF>=0.001'"
    }

    Map[String, String] genetic_maps_dict = read_map(genetic_maps_tsv)


    call CreateShards { input:
        region = region,
        bin_size = filter_and_concat_shard_size,
        pad_size = 0,
        output_prefix = output_prefix + ".shards"
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
            filter_and_concat_short_view_args = filter_and_concat_short_view_args,
            filter_and_concat_short_filter_args = filter_and_concat_short_filter_args,
            filter_and_concat_sv_view_args = filter_and_concat_sv_view_args
        }
        
        # added variant collision fix step
        call FixVariantCollisions { input:
            phased_vcf = FilterAndConcatVcfs.filter_and_concat_vcf,
            fix_variant_collisions_java = fix_variant_collisions_java,
            operation = operation,
            weight_tag = weight_tag,
            is_weight_format_field = is_weight_format_field,
            default_weight = default_weight,
            output_prefix = output_prefix + ".collisionless.shard-" + s,
        }
    }

    call BcftoolsConcatNaive as ConcatFixVariantCollisionsBeforeShapeit { input:
        vcfs = FixVariantCollisions.phased_collisionless_vcf,
        vcf_idxs = FixVariantCollisions.phased_collisionless_vcf_idx,
        output_prefix = output_prefix + ".collisionless"
    }

    call CreateShapeitChunks { input:
        vcf = ConcatFixVariantCollisionsBeforeShapeit.concatenated_vcf,
        vcf_idx = ConcatFixVariantCollisionsBeforeShapeit.concatenated_vcf_idx,
        region = region,
        extra_args = chunk_extra_args
    }

    Array[String] common_regions = read_lines(CreateShapeitChunks.common_chunks) # for shapeit4

    scatter (i in range(length(common_regions))) {
        if (!do_shapeit5) {
            # phase all using Shapeit4
            call Shapeit4 as Shapeit4All { input:
                vcf = ConcatFixVariantCollisionsBeforeShapeit.concatenated_vcf,
                vcf_idx = ConcatFixVariantCollisionsBeforeShapeit.concatenated_vcf_idx,
                genetic_map = genetic_maps_dict[chromosome],
                region = common_regions[i],
                output_prefix = output_prefix + ".phased",
                extra_args = shapeit4_extra_args
            }
        }
        if (do_shapeit5) {
            call FilterCommon { input:
                vcf = ConcatFixVariantCollisionsBeforeShapeit.concatenated_vcf,
                vcf_idx = ConcatFixVariantCollisionsBeforeShapeit.concatenated_vcf_idx,
                output_prefix = output_prefix + ".common",
                region = common_regions[i],
                filter_common_args = filter_common_args
            }
            # phase common scaffold using Shapeit4
            call Shapeit4 as Shapeit4Common { input:
                vcf = FilterCommon.common_vcf,
                vcf_idx = FilterCommon.common_vcf_idx,
                genetic_map = genetic_maps_dict[chromosome],
                region = common_regions[i],
                output_prefix = output_prefix + ".phased",
                extra_args = shapeit4_extra_args
            }
        }
    }

    call LigateVcfs as LigateScaffold { input:
        vcfs = select_all(flatten([Shapeit4All.phased_vcf, Shapeit4Common.phased_vcf])),
        vcf_idxs = select_all(flatten([Shapeit4All.phased_vcf_idx, Shapeit4Common.phased_vcf_idx])),
        output_prefix = output_prefix + ".phased.ligated"
    }

    if (do_shapeit5) {
        # phase rare using Shapeit5
        Array[String] rare_regions = read_lines(CreateShapeitChunks.rare_chunks) # for shapeit5

        scatter (i in range(length(rare_regions))) {
            call Shapeit5Rare { input:
                vcf = ConcatFixVariantCollisionsBeforeShapeit.concatenated_vcf,
                vcf_idx = ConcatFixVariantCollisionsBeforeShapeit.concatenated_vcf_idx,
                scaffold_vcf = LigateScaffold.ligated_vcf,
                scaffold_vcf_idx = LigateScaffold.ligated_vcf_idx,
                genetic_map = genetic_maps_dict[chromosome],
                region = rare_regions[i],
                scaffold_region = common_regions[i],
                output_prefix = output_prefix + ".phased." + "chunk-" + i,
                extra_args = shapeit5_extra_args
            }
        }

        call BcftoolsConcatNaive as ConcatShapeit5 { input:
            vcfs = flatten([Shapeit5Rare.phased_vcf]),
            vcf_idxs = flatten([Shapeit5Rare.phased_vcf_idx]),
            output_prefix = output_prefix + ".phased.concat"
        }
    }

    output {
        File phased_vcf = select_first([ConcatShapeit5.concatenated_vcf, LigateScaffold.ligated_vcf])
        File phased_vcf_idx = select_first([ConcatShapeit5.concatenated_vcf_idx, LigateScaffold.ligated_vcf_idx])
    }
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

    Int disk_gb = 3 * ceil(size([vcf, vcf_idx], "GiB"))

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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
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

        Int disk_gb = 10

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
        mem_gb:             3,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            false,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
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
        String? filter_and_concat_short_filter_args
        String filter_and_concat_short_view_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))<50'"
        String filter_and_concat_sv_view_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))>=50'"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 4 * ceil(size(short_vcf, "GiB") + size(sv_vcf, "GiB"))

    command <<<
        set -euxo pipefail

        # split to biallelic and filter short (re-fill tags when needed)
        bcftools norm --no-version -r ~{region} \
                -m-any -N -f ~{reference_fasta} ~{short_vcf} -Ou | \
            bcftools +fill-tags --no-version -- -t AF,AC,AN -Ou | \
            bcftools filter --no-version ~{filter_and_concat_short_filter_args} -Ou | \
            bcftools +fill-tags --no-version -- -t AF,AC,AN -Ou | \
            bcftools view --no-version ~{filter_and_concat_short_view_args} \
                -Ob -o ~{output_prefix}.short.bcf
        bcftools index ~{output_prefix}.short.bcf

        # populate missing with hom-ref and filter SV
        bcftools +setGT -r ~{region} ~{sv_vcf} --no-version -- -t . -n 0p -Ou | \
            bcftools +fill-tags --no-version -- -t AF,AC,AN -Ou | \
            bcftools view --no-version ~{filter_and_concat_sv_view_args} \
                -Ob -o ~{output_prefix}.SV.bcf
        bcftools index ~{output_prefix}.SV.bcf

        # concatenate with deduplication; providing SV VCF as first argument preferentially keeps those records
        bcftools concat --no-version \
            ~{output_prefix}.SV.bcf \
            ~{output_prefix}.short.bcf \
            --allow-overlaps --remove-duplicates -Ou | \
            bcftools sort -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    >>>

    output {
        File filter_and_concat_vcf = "~{output_prefix}.vcf.gz"
        File filter_and_concat_vcf_idx = "~{output_prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             6,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
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
        File phased_vcf                     # biallelic
        File fix_variant_collisions_java
        Int operation = 1                   # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "SCORE"         # ID of the weight field; weights are assumed to be non-negative; we set to SCORE to prefer kanpig records (and moreover, those with higher SCORE) over DeepVariant records (these should have no SCORE, and will be assigned the low default_weight below)
        Int is_weight_format_field = 0      # given a VCF record in a sample, assign it a weight encoded in the INFO field (0) or in the sample column (1)
        Float default_weight = 0.5          # default weight if the weight field is not found
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 50 + 4 * (ceil(size(phased_vcf, "GiB")))

    command <<<
        set -euxo pipefail

        java ~{fix_variant_collisions_java} \
            ~{phased_vcf} \
            ~{operation} \
            ~{weight_tag} \
            ~{is_weight_format_field} \
            ~{default_weight} \
            collisionless.vcf \
            windows.txt \
            histogram.txt \
            null                            # do not output figures

        # replace all missing alleles (correctly) emitted with reference alleles, since this is expected by PanGenie panel-creation script
        bcftools +setGT --no-version collisionless.vcf -- -t . -n 0p -Ou | \
            bcftools +fill-tags --no-version --threads $(($(nproc)-1)) -Oz -o ~{output_prefix}.phased.collisionless.vcf.gz -- -t AF,AC,AN
        # use vcf.gz to avoid errors from missing header lines
        bcftools index -t ~{output_prefix}.phased.collisionless.vcf.gz
    >>>

    output {
        File phased_collisionless_vcf = "~{output_prefix}.phased.collisionless.vcf.gz"
        File phased_collisionless_vcf_idx = "~{output_prefix}.phased.collisionless.vcf.gz.tbi"
        File windows = "windows.txt"
        File histogram = "histogram.txt"
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
        docker:             "us.gcr.io/broad-gatk/gatk:4.6.0.0"     # needs Java + bcftools
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
        String extra_args = "--thread $(nproc) --window-size 5000000 --buffer-size 500000"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size([vcf, vcf_idx], "GiB")) + 1

    command <<<
        set -euxo pipefail

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_chunk_static
        chmod +x GLIMPSE_chunk_static

        ./GLIMPSE_chunk_static \
            -I ~{vcf} \
            --region ~{region} \
            ~{extra_args} \
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

    Int disk_gb = 50 + 4 * ceil(size(vcfs, "GiB"))

    command <<<
        set -euxo pipefail

        bcftools concat --no-version ~{sep=" " vcfs} --naive -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    >>>

    output {
        File concatenated_vcf = "~{output_prefix}.vcf.gz"
        File concatenated_vcf_idx = "~{output_prefix}.vcf.gz.tbi"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
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

task FilterCommon {
    input {
        File vcf
        File vcf_idx
        String output_prefix
        String region
        String filter_common_args 

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 50 + 4 * ceil(size(vcf, "GiB"))

    command <<<
        set -euxo pipefail

        # filter to common
        bcftools +fill-tags --no-version -r ~{region} ~{vcf} -- -t AF,AC,AN -Ou | \
            bcftools view ~{filter_common_args} \
                -Ob -o ~{output_prefix}.common.bcf
        bcftools index ~{output_prefix}.common.bcf
    >>>

    output {
        File common_vcf = "~{output_prefix}.common.bcf"
        File common_vcf_idx = "~{output_prefix}.common.bcf.csi"
    }
    
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
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

task Shapeit4 {
    input {
        File vcf
        File vcf_idx
        File genetic_map
        String region
        String output_prefix
        String extra_args = "--thread $(nproc) --use-PS 0.0001"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 4 * ceil(size(vcf, "GiB"))
    
    command <<<
        set -euxo pipefail

        shapeit4.2 --input ~{vcf} \
                --map ~{genetic_map} \
                --region ~{region} \
                --sequencing \
                --output ~{output_prefix}.bcf \
                ~{extra_args}
        bcftools index ~{output_prefix}.bcf
    >>>

    output{
        File phased_vcf = "~{output_prefix}.bcf"
        File phased_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             16,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/shapeit4:v1"
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

task LigateVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 4 * ceil(size(vcfs, "GiB"))

    command <<<
        set -euxo pipefail

        ligate_static --input ~{write_lines(vcfs)} --output ~{output_prefix}.bcf
        bcftools +fill-tags --no-version ~{output_prefix}.bcf \
            -Oz -o ~{output_prefix}.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.vcf.gz
    >>>

    output {
        File ligated_vcf = "~{output_prefix}.vcf.gz"
        File ligated_vcf_idx = "~{output_prefix}.vcf.gz.tbi"
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
        docker:             "hangsuunc/shapeit5:v1"
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

task Shapeit5Rare {
    input {
        File vcf
        File vcf_idx
        File scaffold_vcf
        File scaffold_vcf_idx
        File genetic_map
        String region
        String scaffold_region
        String output_prefix
        String extra_args = "--thread $(nproc)"
        
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 4 * ceil(size(vcf, "GiB")) + ceil(size(scaffold_vcf, "GiB"))

    command <<<
        set -euxo pipefail

        # we only need to fill rare in input (common in scaffold should have been imputed or filled previously);
        # this also fills common in input, but those records will be ignored by Shapeit5
        bcftools +setGT --no-version ~{vcf} -- -t . -n 0p -Ou | \
            bcftools +fill-tags --no-version -Ob -o input.bcf -- -t AF,AC,AN
        bcftools index input.bcf

        /shapeit5/phase_rare --input input.bcf \
            --scaffold ~{scaffold_vcf} \
            --map ~{genetic_map} \
            --input-region ~{region} \
            --scaffold-region ~{scaffold_region} \
            --output phased.bcf \
            ~{extra_args}

        bcftools +fill-tags phased.bcf --no-version -Oz -o ~{output_prefix}.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.vcf.gz

    >>>

    output{
        File phased_vcf = "~{output_prefix}.vcf.gz"
        File phased_vcf_idx = "~{output_prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/shapeit5:develop"
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
