/*
   VARIANT CALLING NEXTFLOW WORKFLOW
   SARS-COV2 PROJECT
*/

nextflow.enable.dsl=2

// Pipeline Input Parameters
params.outdir = 'results'
params.genome = "/Users/christophertarkaa/sars-cov2-project/Reference/MN908947.3.fasta"
params.reads = "/Users/christophertarkaa/sars-cov2-project/covid_samples/*_{1,2}.fastq"
params.minQuality = 20
params.minLength = 50

// CHANNELS
ref_ch = Channel.fromPath(params.genome, checkIfExists: true)
reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

workflow {
    BWA_INDEX(ref_ch)
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
