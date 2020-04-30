workflow CC_Crams{
    Array[File] crams
    Array[File] crais
    Array[String] samples
    File hg38_reference
    File hg38_reference_fai
    File hg38_reference_dict
    File intervals_exons
    scatter (scatter_index in range(length(crams))){
        call cram2bam {
            input:
                intervals_exons = intervals_exons,
                cram = crams[scatter_index],
                crai = crais[scatter_index],
                hg38_reference = hg38_reference,
                hg38_reference_fai = hg38_reference_fai,
                hg38_reference_dict = hg38_reference_dict,
                sample = samples[scatter_index]
        }
    }
    output {
        Array [String] sample = cram2bam.entity_id
        Array [File] bams = cram2bam.bam
        Array [File] bais = cram2bam.bai
    }

    meta {
        author: "Jack Fu"
        email: "jmfu@mgh.harvard.edu"
    }
}
task cram2bam {
    File intervals_exons
    File cram
    File crai
    File hg38_reference
    File hg38_reference_fai
    File hg38_reference_dict
    String sample
    # Runtime parameters
    Int? disk_space_gb
    Boolean use_ssd = false
    Int cpu=1
    Int? preemptible_attempts
    Int machine_mem_mb = 600
    Int command_mem_mb = machine_mem_mb - 100
    String base_filename = "${sample}"
    String counts_exons_filename = "${base_filename}.exons.counts.tsv"
    command <<<
        set -e
        set -o pipefail
        ln -vs ${hg38_reference} reference.fasta
        ln -vs ${hg38_reference_fai} reference.fasta.fai
        ln -vs ${hg38_reference_dict} reference.dict
        samtools view -h -T reference.fasta ${cram} |
        samtools view -b -o ${sample}.bam -
        samtools index -b ${sample}.bam
        mv ${sample}.bam.bai ${sample}.bai
    >>>
    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(cram, "GB")) + 200]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
        maxRetries: 3
    }
    output {
        String entity_id = base_filename
        File bam = "${sample}.bam"
        File bai = "${sample}.bai"
    }
}