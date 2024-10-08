version 1.0

import "./Helper.wdl" as H

workflow SplitCohortVcf {

    input {
        File joint_short_vcf
        File joint_short_vcf_tbi
        File ref_dict
        Array[String] region_list
        Array[String] sample_list
        String outputdirectory
    }

    scatter (region in region_list)  {

        call H.SubsetVCF as SubVCF { input:
            vcf_gz = joint_short_vcf,
            vcf_tbi = joint_short_vcf_tbi,
            locus = region
        }

        call SplitVcf { input:
            joint_vcf = SubVCF.subset_vcf,
            joint_vcf_tbi = SubVCF.subset_tbi
        }

        call reorder_samples{ input:
            sample_names = sample_list,
            file_paths = SplitVcf.vcf_by_sample
        }

    }
    # output Array[Array[File]]
    Array[Int] indexes = range(length(sample_list))
    scatter (ind in indexes) {
        String sample_id = sample_list[ind]
        # flatten the array file
        scatter (array in reorder_samples.sample_map){
            File array_name = array[ind]
        }

        ### merge per_chromosome vcfs
        call MergePerChrCalls {input: 
            vcfs = array_name,
            ref_dict = ref_dict,
            prefix = sample_id
        }

    }

    call H.FinalizeToDir { input:
        files = MergePerChrCalls.vcf,
        outdir = outputdirectory
    }

    output {
        Array[File] splitted_vcf = MergePerChrCalls.vcf
    }
}


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


task SplitVcf {

    input {
        File joint_vcf
        File joint_vcf_tbi
        String chromo
        Int memory
    }

    command <<<
        set -euxo pipefail
        mkdir output
        bcftools +split -Oz -o output ~{joint_vcf}
        cd output
        for vcf in $(find . -name "*.vcf.gz"); do
            filename=$(basename "$vcf")
            mv "$vcf" "~{chromo}.$filename"
        done
        cd -

    >>>

    output {
        Array[File] vcf_by_sample = glob("output/*vcf.gz")
        # Array[File] vcf_by_sample_tbi = glob("output/*vcf.gz.tbi")

    }
    ###################
    runtime {
        cpu: 4
        memory: memory + " GiB"
        disks: "local-disk 1125 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}

task reorder_samples {
  input{
    # Command parameters
    Array[String] sample_names
    Array[String] file_paths
    String outfile 
    
    # Runtime parameters
    String docker = "python:latest"
    Int machine_mem_gb = 7
    Int disk_space_gb = 100
    Int preemptible_attempts = 0
  }
    command <<<
    set -oe pipefail
    
    python << CODE
    file_paths = ['~{sep="','" file_paths}']
    sample_names = ['~{sep="','" sample_names}']

    if len(file_paths)!= len(sample_names):
      print("Number of File Paths does not Equal Number of File Names")
      exit(1)
    
    reordered_filepath = []
    for sample in sample_names:
        for file in file_paths:
            if sample in file:
                reordered_filepath.append(file)
                break
    with open("reordered_file.txt", "w") as fi:
      for fp in reordered_filepath:
        fi.write(fp + "\n") 


    CODE

    mv reordered_file.txt ~{outfile}
    >>>

    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        disks: "local-disk " + disk_space_gb + " HDD"
        preemptible: 3
    }

    output {
        Array[File] sample_map = read_lines(outfile)
    }
}

task MergePerChrCalls {

    meta {
        description: "Merge per-chromosome calls into a single VCF"
    }

    parameter_meta {
        vcfs: "List of per-chromosome VCFs to merge"
        ref_dict: "Reference dictionary"
        prefix: "Prefix for output VCF"
    }

    input {
        Array[File] vcfs
        File ref_dict
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(vcfs, "GB")) + 1

    command <<<
        set -euxo pipefail

        VCF_WITH_HEADER=~{vcfs[0]}
        GREPCMD="grep"
        if [[ ~{vcfs[0]} =~ \.gz$ ]]; then
            GREPCMD="zgrep"
        fi

        $GREPCMD '^#' $VCF_WITH_HEADER | grep -v -e '^##contig' -e CHROM > header
        grep '^@SQ' ~{ref_dict} | awk '{ print "##contig=<ID=" $2 ",length=" $3 ">" }' | sed 's/[SL]N://g' >> header
        $GREPCMD -m1 CHROM $VCF_WITH_HEADER >> header

        ((cat header) && ($GREPCMD -h -v '^#' ~{sep=' ' vcfs})) | bcftools sort | bgzip > ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File vcf = "~{prefix}.vcf.gz"
        File tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             24,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
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