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

workflow AnnotateAndPop {
    input {
        Array[File] posteriors_vcf_gzs          # whole-genome, per-batch, biallelic GLIMPSE2
        Array[File] posteriors_vcf_gz_tbis
        Array[File] panel_split_vcf_gzs         # per-chromosome
        Array[File] panel_split_vcf_gz_tbis
        Array[File] panel_id_split_vcf_gzs      # per-chromosome
        Array[File] panel_id_split_vcf_gz_tbis
        Array[String]+ chromosomes

        Array[String] output_prefixes

        String docker
    }

    scatter (i in range(length(posteriors_vcf_gzs)))
    {
        scatter (j in range(length(chromosomes))) {
            call AnnotateAndPop as ChromosomeAnnotateAndPop {
                input:
                    posteriors_vcf_gz = posteriors_vcf_gzs[i],
                    posteriors_vcf_gz_tbi = posteriors_vcf_gz_tbis[i],
                    panel_split_vcf_gz = panel_split_vcf_gzs[j],
                    panel_split_vcf_gz_tbi = panel_split_vcf_gz_tbis[j],
                    panel_id_split_vcf_gz = panel_id_split_vcf_gzs[j],
                    panel_id_split_vcf_gz_tbi = panel_id_split_vcf_gz_tbis[j],
                    chromosome = chromosomes[j],
                    output_prefix = output_prefixes[i] + "." + chromosomes[j],
                    docker = docker
            }
        }

        # concat across chromosomes
        call ConcatVcfs {
            input:
                vcf_gzs = ChromosomeAnnotateAndPop.annotated_and_popped_vcf_gz,
                vcf_gz_tbis = ChromosomeAnnotateAndPop.annotated_and_popped_vcf_gz_tbi,
                output_prefix = output_prefixes[i],
                docker = docker
        }
    }

    output {
        Array[File] annotated_and_popped_vcf_gzs = ConcatVcfs.vcf_gz
        Array[File] annotated_and_popped_vcf_gz_tbis = ConcatVcfs.vcf_gz_tbi
    }
}

task AnnotateAndPop {
    input {
        # all VCFs should be split to biallelic
        File posteriors_vcf_gz
        File posteriors_vcf_gz_tbi
        File panel_split_vcf_gz
        File panel_split_vcf_gz_tbi
        File panel_id_split_vcf_gz
        File panel_id_split_vcf_gz_tbi
        String chromosome
        String output_prefix

        String docker
        File? monitoring_script
        RuntimeAttributes runtime_attributes = {"use_ssd": true}
    }

    Int disk_size_gb = 3 * ceil(size(posteriors_vcf_gz , "GB"))

    command <<<
        set -euox pipefail

        bcftools annotate -r ~{chromosome} -a ~{panel_split_vcf_gz} ~{posteriors_vcf_gz} \
            -c CHROM,POS,REF,ALT,ID:=INFO/ID,INFO/ID:=INFO/ID \
            -Oz -o ~{output_prefix}.annotated.vcf.gz
        bcftools index -t ~{output_prefix}.annotated.vcf.gz

        # modified version of convert-to-biallelic.py
        # DO NOT apply bcftools norm -m+ before using this, pass a biallelic VCF instead!
        python - --panel_id_split_vcf_gz ~{panel_id_split_vcf_gz} \
                 --input_vcf_gz ~{output_prefix}.annotated.vcf.gz \
                 --output_vcf /tmp/~{output_prefix}.annotated-popped.vcf \
                 <<-'EOF'
        import sys
        import argparse
        from collections import defaultdict
        import gzip
        import tqdm

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--panel_id_split_vcf_gz',
                                type=str)
            parser.add_argument('--input_vcf_gz',
                                type=str)
            parser.add_argument('--output_vcf',
                                type=str)

            args = parser.parse_args()

            # chromosome ->  ID -> [start, REF, ALT] per chromosome
            chrom_to_variants = defaultdict(lambda: defaultdict(list))

            # read the biallelic panel VCF containing REF/ALT for all variant IDs and store them in list of assigned IDs
            for line in gzip.open(args.panel_id_split_vcf_gz, 'rt'):
                if line.startswith('#'):
                    continue
                fields = line.split()
                info_field = { i.split('=')[0] : i.split('=')[1] for i in fields[7].split(';') if "=" in i}
                assert 'ID' in info_field
                ids = info_field['ID'].split(',')
                assert len(ids) == 1
                chrom_to_variants[fields[0]][ids[0]] = [fields[1], fields[3], fields[4]]

            with open(args.output_vcf, 'w') as f:
                for line in tqdm.tqdm(gzip.open(args.input_vcf_gz, 'rt')):
                    if line.startswith('#'):
                        # header line
                        #if any([i in line for i in ['INFO=<ID=AF', 'INFO=<ID=AK', 'FORMAT=<ID=GL', 'FORMAT=<ID=KC']]):
                        #    # these fields will not be contained in biallelic VCF
                        #    continue
                        f.write(line[:-1] + '\n')
                        continue
                    fields = line.split()
                    assert len(fields) > 7
                    # parse the INFO field
                    info_field = { i.split('=')[0] : i.split('=')[1] for i in fields[7].split(';') if "=" in i}
                    assert 'ID' in info_field
                    # determine ID string belonging to each allele (keep empty string for REF, as it does not have an ID)
                    allele_to_ids = [''] + info_field['ID'].split(',')
                    info_ids = info_field['ID'].split(',')
                    # allow bi-allelic records with unknown IDs (that are not in annotation VCF)
                    # if (len(info_ids) == 1) and any([x not in chrom_to_variants[fields[0]] for x in info_ids[0].split(':')]):
                    #    # unknown ID, leave record as is
                    #    print(line[:-1])
                    #    continue
                    # collect all variant IDs in this region
                    ids = set([])
                    if len(info_field['ID'].split(',')) > 1:
                        raise 'VCF should be biallelic'
                    bubble_id = info_field['ID']
                    for assigned_id in bubble_id.split(':'):
                        variants = chrom_to_variants[fields[0]][assigned_id]
                        if variants:
                            ids.add((assigned_id, int(variants[0])))
                    # sort the ids by the starting coordinate (to ensure the VCF is sorted)
                    ids = list(ids)
                    ids.sort(key=lambda x : x[1])
                    # create a single, biallelic VCF record for each ID
                    for (var_id, coord) in ids:
                        vcf_line = fields[:9]
                        # set start coordinate
                        vcf_line[1] = str(coord)
                        # also add assigned ID to ID column of the VCF
                        vcf_line[2] = var_id
                        # set REF
                        vcf_line[3] = chrom_to_variants[fields[0]][var_id][1]
                        # set ALT
                        vcf_line[4] = chrom_to_variants[fields[0]][var_id][2]
                        # set INFO/ID to bubble ID
                        vcf_line[7] = 'ID=' + bubble_id
                        # also add other INFO fields (except ID which was replaced)
                        for k,v in info_field.items():
                            if k == 'ID':
                                continue
                            if k in ['MA', 'UK']:
                                values = ';' + k + '=' + v
                                vcf_line[7] = vcf_line[7] + values
                        vcf_line[8] = fields[8]
                        # determine the genotype of each sample
                        for sample_field in fields[9:]:
                            # determine position of GT from FORMAT
                            assert 'GT' in fields[8]
                            format_field = fields[8].split(':')
                            index_of_gt = format_field.index('GT')
                            genotype = sample_field.split(':')
                            biallelic_genotype = []
                            for allele in genotype[index_of_gt].split('|'):
                                if allele == '.':
                                    # missing allele
                                    biallelic_genotype.append('.')
                                else:
                                    if var_id in allele_to_ids[int(allele)].split(':'):
                                        biallelic_genotype.append('1')
                                    else:
                                        biallelic_genotype.append('0')
                            # append other FORMAT fields
                            vcf_line.append('|'.join(biallelic_genotype) + ':' + ':'.join(genotype[1:]))
                        f.write('\t'.join(vcf_line) + '\n')

        if __name__ == '__main__':
            main()
        EOF

        bcftools sort /tmp/~{output_prefix}.annotated-popped.vcf -Oz -o ~{output_prefix}.annotated-popped.vcf.gz
        bcftools index -t ~{output_prefix}.annotated-popped.vcf.gz
    >>>

    output {
        File annotated_and_popped_vcf_gz = "~{output_prefix}.annotated-popped.vcf.gz"
        File annotated_and_popped_vcf_gz_tbi = "~{output_prefix}.annotated-popped.vcf.gz.tbi"
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
}

task ConcatVcfs {
    input {
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {"use_ssd": true}
    }

    Int disk_size_gb = 3 * ceil(size(vcf_gzs, "GB"))

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        mkdir inputs
        mv ~{sep=' ' vcf_gzs} inputs
        mv ~{sep=' ' vcf_gz_tbis} inputs

        # TODO NOTE USE OF SORT, ENSURE THIS GIVES DESIRED CHROMOSOME ORDER (COULD TAKE IN ARRAY OF CHROMOSOMES INSTEAD)
        if [ $(ls inputs/*.vcf.gz | wc -l) == 1 ]
        then
            cp $(ls inputs/*.vcf.gz) ~{output_prefix}.vcf.gz
            cp $(ls inputs/*.vcf.gz.tbi) ~{output_prefix}.vcf.gz.tbi
        else
            bcftools concat $(ls inputs/*.vcf.gz | sort -V -d) --naive -Oz -o ~{output_prefix}.vcf.gz
            bcftools index -t ~{output_prefix}.vcf.gz
        fi
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
