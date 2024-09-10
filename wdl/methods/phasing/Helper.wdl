version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

struct DataTypeParameters {
    Int num_shards
    String map_preset
}

task HiPhase {

    meta {
        description: "Generates phased VCF. Note this runs fast so no need to parallize."
    }


    input {
        File bam
        File bai

        File unphased_snp_vcf
        File unphased_snp_tbi
        File unphased_sv_vcf
        File unphased_sv_tbi

        File ref_fasta
        File ref_fasta_fai
        String samplename

        Int memory
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        String extra_args

        RuntimeAttr? runtime_attr_override
    }

    # Int bam_sz = ceil(size(bam, "GB"))
	Int disk_size = 30 # if bam_sz > 200 then 2*bam_sz else bam_sz + 200
    Int thread_num = memory/2

    command <<<
        set -euxo pipefail

        touch ~{bai}

        hiphase \
        --threads ~{thread_num} \
        --bam ~{bam} \
        --reference ~{ref_fasta} \
        --global-realignment-cputime 300 \
        --vcf ~{unphased_snp_vcf} \
        --output-vcf ~{samplename}_phased_snp.vcf.gz \
        --vcf ~{unphased_sv_vcf} \
        --output-vcf ~{samplename}_phased_sv.vcf.gz \
        --haplotag-file ~{samplename}_phased_sv_haplotag.tsv \
        --stats-file ~{samplename}.stats.csv \
        --blocks-file ~{samplename}.blocks.tsv \
        --summary-file ~{samplename}.summary.tsv \
        --sample-name ~{samplename} \
        ~{extra_args}

        bcftools sort ~{samplename}_phased_snp.vcf.gz -O z -o ~{samplename}_phased_snp.sorted.vcf.gz
        tabix -p vcf ~{samplename}_phased_snp.sorted.vcf.gz

        bcftools sort ~{samplename}_phased_sv.vcf.gz -O z -o ~{samplename}_phased_sv.sorted.vcf.gz
        tabix -p vcf ~{samplename}_phased_sv.sorted.vcf.gz
        
    >>>

    output {
        File phased_snp_vcf = "~{samplename}_phased_snp.sorted.vcf.gz"
        File phased_snp_vcf_tbi = "~{samplename}_phased_snp.sorted.vcf.gz.tbi"
        File phased_sv_vcf   = "~{samplename}_phased_sv.sorted.vcf.gz"
        File phased_sv_vcf_tbi = "~{samplename}_phased_sv.sorted.vcf.gz.tbi"
        File haplotag_file = "~{samplename}_phased_sv_haplotag.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          thread_num,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "hangsuunc/hiphase:1.3.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task SubsetVCF {

    meta {
        description: "Subset a VCF file to a given locus"
    }

    parameter_meta {
        vcf_gz: "VCF file to be subsetted"
        vcf_tbi: "Tabix index for the VCF file"
        locus: "Locus to be subsetted"
        prefix: "Prefix for the output file"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File vcf_gz
        File vcf_tbi
        String locus
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf_gz, vcf_tbi], "GB")) + 1

    command <<<
        set -euxo pipefail

        bcftools view ~{vcf_gz} --regions ~{locus} | bgzip > ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longshot:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task SubsetBam {

    meta {
        description : "Subset a BAM file to a specified locus."
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        File bai
        String locus
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }



    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam
    >>>

    output {
        File subset_bam = "~{prefix}.bam"
        File subset_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task InferSampleName {
    meta {
        description: "Infer sample name encoded on the @RG line of the header section. Fails if multiple values found, or if SM ~= unnamedsample."
    }

    parameter_meta {
        bam: {
            localization_optional: true,
            description: "BAM file"
        }
    }

    input {
        File bam
        File bai
    }



    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} > header.txt
        if ! grep -q '^@RG' header.txt; then echo "No read group line found!" && exit 1; fi

        grep '^@RG' header.txt | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g' | sort | uniq > sample.names.txt
        if [[ $(wc -l sample.names.txt) -gt 1 ]]; then echo "Multiple sample names found!" && exit 1; fi
        if grep -iq "unnamedsample" sample.names.txt; then echo "Sample name found to be unnamedsample!" && exit 1; fi
    >>>

    output {
        String sample_name = read_string("sample.names.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker:         "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task SplitVCFbySample {
    input{       
        File joint_vcf
        String region
        String samplename
    }
    
    command <<<
        set -x pipefail

        bcftools index ~{joint_vcf}

        bcftools view -s ~{samplename} ~{joint_vcf} -o ~{samplename}.subset.g.vcf.gz

        tabix -p vcf ~{samplename}.subset.g.vcf.gz

    >>>
    
    output {
		File single_sample_vcf = "~{samplename}.subset.g.vcf.gz"
        File single_sample_vcf_tbi = "~{samplename}.subset.g.vcf.gz.tbi"
    }


    Int disk_size = 1 + ceil(2 * (size(joint_vcf, "GiB")))

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task MergePerChrVcfWithBcftools {
    parameter_meta {
        vcf_input: {localization_optional: true}
        tbi_input: {localization_optional: true}
    }
    input{
        Array[File] vcf_input
        Array[File] tbi_input
        String pref
        Int threads_num
    }

    command <<<
        set -eux

        # we do single-sample phased VCFs localization ourselves
        mkdir -p ssp_vcfs
        time \
        gcloud storage cp ~{sep=" " vcf_input} /cromwell_root/ssp_vcfs/

        time \
        gcloud storage cp ~{sep=" " tbi_input} /cromwell_root/ssp_vcfs/

        # then merge, and safely assume all ssp-VCFs are sorted in the same order, on one chr
        cd ssp_vcfs
        ls *.vcf.gz > my_vcfs.txt

        bcftools merge \
            --threads ~{threads_num} \
            --merge none \
            -l my_vcfs.txt \
            -O z \
            -o ~{pref}.AllSamples.vcf.gz

        tabix -@ ~{threads_num} -p vcf ~{pref}.AllSamples.vcf.gz

        # move result files to the correct location for cromwell to de-localize
        mv ~{pref}.AllSamples.vcf.gz ~{pref}.AllSamples.vcf.gz.tbi /cromwell_root/
    >>>

    output{
        File merged_vcf = "~{pref}.AllSamples.vcf.gz"
        File merged_tbi = "~{pref}.AllSamples.vcf.gz.tbi"
    }

    runtime {
        cpu: 16
        memory: "32 GiB"
        disks: "local-disk 375 LOCAL"
        preemptible: 1
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}

task Shapeit4 {
    input{
        File vcf_input
        File vcf_index
        File mappingfile
        String region
        String prefix
        Int num_threads
        Int memory
        String extra_args

        RuntimeAttr? runtime_attr_override
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    command <<<
        # add AN AC tag

        # export MONITOR_MOUNT_POINT="/cromwell_root/"
        # bash /opt/vm_local_monitoring_script.sh &> resources.log &
        # job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        shapeit4.2 --input ~{vcf_input} \
                --map ~{mappingfile} \
                --region ~{region} \
                --sequencing \
                --output ~{prefix}.bcf \
                --thread ~{num_threads} \
                ~{extra_args}

        # if ps -p "${job_id}" > /dev/null; then kill "${job_id}"; fi
    >>>

    output{
        # File resouce_monitor_log = "resources.log"
        File phased_bcf = "~{prefix}.bcf"
    }

    #Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

 #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_threads,
        mem_gb:             memory,
        disk_gb:            100,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "hangsuunc/shapeit4:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
