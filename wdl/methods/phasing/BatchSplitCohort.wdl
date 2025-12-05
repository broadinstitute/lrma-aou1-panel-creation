version 1.0

import "./Helper.wdl" as H

workflow SplitCohortVcf {

    input {
        File joint_short_vcf
        File joint_short_vcf_tbi
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
            joint_vcf_tbi = SubVCF.subset_tbi,
            chromo = region
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
        call ConcateAndSortVCFs {input: 
            vcfs = array_name,
            prefix = sample_id
        }

    }

    call H.FinalizeToDir { input:
        files = ConcateAndSortVCFs.vcf,
        outdir = outputdirectory
    }

    output {
        Array[File] splitted_vcf = ConcateAndSortVCFs.vcf
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
        disks: "local-disk 375 LOCAL"
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

task ConcateAndSortVCFs {

    meta {
        description: "Fast merging & sorting VCFs when the default sorting is expected to be slow"
    }

    input {
        Array[File] vcfs
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int sz = ceil(size(vcfs, 'GB'))
    Int disk_sz = if sz > 100 then 5 * sz else 375  # it's rare to see such large gVCFs, for now

    Int cores = 8

    # pending a bug fix (bcftools github issue 1576) in official bcftools release,
    # bcftools sort can be more efficient in using memory
    Int machine_memory = 48 # 96
    Int work_memory = ceil(machine_memory * 0.8)

    command <<<
        set -euxo pipefail

        echo ~{sep=' ' vcfs} | sed 's/ /\n/g' > all_raw_vcfs.txt

        echo "==========================================================="
        echo "starting concatenation" && date
        echo "==========================================================="
        bcftools \
            concat \
            --naive \
            --threads ~{cores-1} \
            -f all_raw_vcfs.txt \
            --output-type z \
            -o concatedated_raw.vcf.gz  # fast, at the expense of disk space
        for vcf in ~{sep=' ' vcfs}; do rm $vcf ; done

        # this is another bug in bcftools that's hot fixed but not in official release yet
        # (see bcftools github issue 1591)
        echo "==========================================================="
        echo "done concatenation, fixing header of naively concatenated VCF" && date
        echo "==========================================================="

        echo "==========================================================="
        echo "starting sort operation" && date
        echo "==========================================================="
        bcftools \
            sort \
            --temp-dir tm_sort \
            --output-type z \
            -o ~{prefix}.vcf.gz \
            concatedated_raw.vcf.gz
        bcftools index --tbi --force ~{prefix}.vcf.gz
        echo "==========================================================="
        echo "done sorting" && date
        echo "==========================================================="
    >>>

    output {
        File vcf = "~{prefix}.vcf.gz"
        File tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cores,
        mem_gb:             "~{machine_memory}",
        disk_gb:            disk_sz,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}