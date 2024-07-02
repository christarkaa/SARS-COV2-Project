/*
   VARIANT CALLING NEXTFLOW WORKFLOW
   SARS-COV2 PROJECT
*/

nextflow.enable.dsl=2

// Pipeline Input Parameters
params.outdir = 'results'
params.genome = "/Users/christophertarkaa/sars-cov2-project/Reference/MN908947.3.fasta"
params.reads = "/Users/christophertarkaa/sars-cov2-project/covid_samples/*_{1,2}.fastq.gz"
params.minQuality = 20
params.minLength = 50

// CHANNELS
ref_Ch = Channel.fromPath(params.genome, checkIfExists: true)
reads_ch = Channel.fromFilePairs(params.reads, checkIfExist: true)

// WORKFLOW
workflow {

    // Quality control
    reads_ch | FASTQC

    // Trimming of reads
    trimmed_reads_ch = reads_ch | TRIMMING

    // Index the reference genome
    indexed_genome_ch = BWA_INDEX(ref_Ch)

    // Align reads to reference genome
    aligned_bams_ch = BWA_ALIGN(indexed_genome_ch.out.bwa_index, trimmed_reads_ch)

    // Sort and index BAM files
    sorted_bams_ch = SORT_AND_INDEX_BAM(aligned_bams_ch)

    // Calculate read coverage (mpileup)
    pileup_ch = BCFTOOLS_MPILEUP(sorted_bams_ch, ref_Ch)

    // Call variants
    raw_vcf_ch = BCFTOOLS_CALL(pileup_ch, ref_Ch)

    // Filter and report SNVs
    filtered_vcf_ch = VCFUTILS(raw_vcf_ch)

    // Additional steps like extracting SNPs and Indels
    SNPs = EXTRACT_SNPS(filtered_vcf_ch)
    INDELs = EXTRACT_INDELS(filtered_vcf_ch)
}

// PROCESSES

// Quality Control
process FASTQC {
    tag { "FASTQC ${sample_id}" }
    label 'process_low'
    cpus 6

    publishDir("${params.outdir}/QC", mode: 'copy')

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*_fastqc*"), emit: fastqc_out

    script:
    """
    fastqc ${reads}
    """
}

// Trimming
process TRIMMING {
    tag { "TRIMMING ${sample_id}" }
    label 'process_low'
    cpus 4

    publishDir("${params.outdir}/Trimming", mode: 'copy')

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1.trimmed.fastq.gz"), path("${sample_id}_2.trimmed.fastq.gz"), emit: trimmed_fastq

    script:
    """
    sickle pe -f ${reads[0]} -r ${reads[1]} -t sanger -o ${sample_id}_1.trimmed.fastq.gz -p ${sample_id}_2.trimmed.fastq.gz -s ${sample_id}_singletons.fastq.gz -q ${params.minQuality} -l ${params.minLength}
    """
}

// Index the reference genome
process BWA_INDEX {
    tag { "BWA_INDEX ${genome}" }
    label 'process_low'
    cpus 6

    publishDir("${params.outdir}/bwa_index", mode: 'copy')

    input:
    path(genome)

    output:
    path("*"), emit: bwa_index

    script:
    """
    bwa index ${genome}
    """
}

// Align reads to reference genome & create BAM file.
process BWA_ALIGN {
    tag { "BWA_ALIGN ${sample_id}" }
    label 'process_medium'
    cpus 6

    publishDir("${params.outdir}/bwa_align", mode: 'copy')

    input:
    tuple path(index), path("*"), emit: bwa_index
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), emit: aligned_bam

    script:
    """
    bwa mem ${index[0]} ${reads[0]} ${reads[1]} | samtools view -Sb - > ${sample_id}.aligned.bam
    """
}

// Sort and index BAM files
process SORT_AND_INDEX_BAM {
    tag { "SORT_AND_INDEX_BAM ${sample_id}" }
    label 'process_medium'
    cpus 4

    publishDir("${params.outdir}/sorted_bam", mode: 'copy')

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: sorted_bam

    script:
    """
    samtools sort -o ${sample_id}.sorted.bam ${bam}
    samtools index ${sample_id}.sorted.bam
    """
}

// Calculate read coverage of positions in the genome
process BCFTOOLS_MPILEUP {
    tag { "BCFTOOLS_MPILEUP ${sample_id}" }
    label 'process_high'
    cpus 8

    publishDir("${params.outdir}/mpileup", mode: 'copy')

    input:
    tuple val(sample_id), path(bam)
    path(reference)

    output:
    path("${sample_id}.mpileup"), emit: pileup

    script:
    """
    bcftools mpileup -f ${reference} ${bam} -o ${sample_id}.mpileup
    """
}

// Detect the single nucleotide variants (SNVs)
process BCFTOOLS_CALL {
    tag { "BCFTOOLS_CALL ${sample_id}" }
    label 'process_high'
    cpus 8

    publishDir("${params.outdir}/variants", mode: 'copy')

    input:
    path(pileup)
    path(reference)

    output:
    path("${sample_id}.raw.vcf"), emit: raw_vcf

    script:
    """
    bcftools call -m -v -Ov -o ${sample_id}.raw.vcf ${pileup}
    """
}

// Filter and report the SNVs in VCF (variant calling format)
process VCFUTILS {
    tag { "VCFUTILS ${sample_id}" }
    label 'process_high'
    cpus 4

    publishDir("${params.outdir}/filtered_vcf", mode: 'copy')

    input:
    path(raw_vcf)

    output:
    path("${sample_id}.filtered.vcf"), emit: filtered_vcf

    script:
    """
    bcftools filter -s LOWQUAL -e '%QUAL<20 || DP<10' ${raw_vcf} -o ${sample_id}.filtered.vcf
    """
}

// Extract SNPs
process EXTRACT_SNPS {
    tag { "EXTRACT_SNPS ${sample_id}" }
    label 'process_low'
    cpus 2

    publishDir("${params.outdir}/snps", mode: 'copy')

    input:
    path(vcf)

    output:
    path("${sample_id}.snps.vcf"), emit: snps

    script:
    """
    bcftools view -v snps ${vcf} -o ${sample_id}.snps.vcf
    """
}

// Extract Indels
process EXTRACT_INDELS {
    tag { "EXTRACT_INDELS ${sample_id}" }
    label 'process_low'
    cpus 2

    publishDir("${params.outdir}/indels", mode: 'copy')

    input:
    path(vcf)

    output:
    path("${sample_id}.indels.vcf"), emit: indels

    script:
    """
    bcftools view -v indels ${vcf} -o ${sample_id}.indels.vcf
    """
}

