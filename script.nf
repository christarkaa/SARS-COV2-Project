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
ref_ch = Channel.fromPath(params.genome, checkIfExists: true)
reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

// WORKFLOW
workflow {
    FASTQC(reads_ch)
    BWA_INDEX(ref_ch)
    BWA_ALIGN(BWA_INDEX.out.bwa_index.combine(reads_ch))
    SAMTOOLS_SORT(BWA_ALIGN.out.aligned_bam)
    BCFTOOLS_MPILEUP(SAMTOOLS_SORT.out.sorted_bam.combine(ref_ch))
    BCFTOOLS_CALL(BCFTOOLS_MPILEUP.out.pileup)
    VCFUTILS(BCFTOOLS_CALL.out.raw_vcf)
    EXTRACT_SNPS(VCFUTILS.out.filtered_vcf)
    EXTRACT_INDELS(VCFUTILS.out.filtered_vcf)
}

// PROCESSES

// Quality Control
process FASTQC {
    tag { "FASTQC ${sample_id}" }
    label 'process_low'
    cpus 4 

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

// Index the reference genome
process BWA_INDEX {
    tag { "BWA_INDEX ${genome}" }
    label 'process_low'
    cpus 4 

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

// Align reads to reference genome & create BAM file
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

// Sort BAM files
process SAMTOOLS_SORT {
    tag { "SORT_BAM ${sample_id}" }
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
    cpus 7 

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
    cpus 7 
    publishDir("${params.outdir}/variants", mode: 'copy')

    input:
    path(pileup)

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
    cpus 2 // Adjusted to fit within the CPU limit

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
