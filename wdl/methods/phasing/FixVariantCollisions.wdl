version 1.0

workflow FixVariantCollisions {
    input {
        File phased_vcf_gz                     # biallelic
        File phased_vcf_gz_tbi
        File fix_variant_collisions_java
        Int operation = 1                   # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "UNIT_WEIGHT"   # ID of the weight field; if this field is not found, all weights are set to one; weights are assumed to be non-negative
        Int is_weight_format_field = 0      # given a VCF record in a sample, assign it a weight encoded in the sample column (1) or in the INFO field (0)
        String output_prefix
    }

    call FixVariantCollisions {
        input:
            phased_vcf_gz = phased_vcf_gz,
            phased_vcf_gz_tbi = phased_vcf_gz_tbi,
            fix_variant_collisions_java = fix_variant_collisions_java,
            operation = operation,
            weight_tag = weight_tag,
            is_weight_format_field = is_weight_format_field,
            output_prefix = output_prefix
    }

    output {
        File phased_collisionless_vcf_gz = FixVariantCollisions.phased_collisionless_vcf_gz
        File phased_collisionless_vcf_gz_tbi = FixVariantCollisions.phased_collisionless_vcf_gz_tbi
        File windows = FixVariantCollisions.windows
        File histogram = FixVariantCollisions.histogram
    }
}

task FixVariantCollisions {

    input {
        File phased_vcf_gz                     # biallelic
        File phased_vcf_gz_tbi
        File fix_variant_collisions_java
        Int operation = 1                   # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "UNIT_WEIGHT"   # ID of the weight field; if this field is not found, all weights are set to one; weights are assumed to be non-negative
        Int is_weight_format_field = 0      # given a VCF record in a sample, assign it a weight encoded in the sample column (1) or in the INFO field (0)
        String output_prefix
    }

    command <<<
        set -euxo pipefail

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
        bcftools view collisionless.vcf | \
            sed -e 's/\.|0/0|0/g' | sed -e 's/0|\./0|0/g' | sed -e 's/\.|1/0|1/g' | sed -e 's/1|\./1|0/g' | sed -e 's/\.|\./0|0/g' | \
            bcftools view -Oz -o ~{output_prefix}.phased.collisionless.vcf.gz
        bcftools index -t ~{output_prefix}.phased.collisionless.vcf.gz
    >>>

    output {
        File phased_collisionless_vcf_gz = "~{output_prefix}.phased.collisionless.vcf.gz"
        File phased_collisionless_vcf_gz_tbi = "~{output_prefix}.phased.collisionless.vcf.gz.tbi"
        File windows = "windows.txt"
        File histogram = "histogram.txt"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "16 GiB"
        disks: "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-gatk/gatk:4.6.0.0"
    }
}
