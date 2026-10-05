#!/usr/bin/env nextflow

/*
 * Nextflow pipeline to test Immcantation Docker container and pipelines
 * Recreates test_pipelines.sh tests
 */

// Parameters
params.image = "immcantation/suite:devel"
params.pythonwarnings = "default"
params.nproc = 2
params.ext = "tsv"

// Get current date and version for output directory
def version = params.image.split(':')[1]
def date = new Date().format('yyyy.MM.dd')

// Set up directories
def data_dir = "${workflow.launchDir}/data"
def run_dir = "${workflow.launchDir}/run/nextflow-${date}-${version}"

// Update params
params.data_dir = data_dir
params.run_dir = run_dir

// Create output directory
file(params.run_dir).mkdirs()

log.info """
    ╔═══════════════════════════════════════════════════════════════
    ║ Immcantation Pipeline Tests (Nextflow)
    ╠═══════════════════════════════════════════════════════════════
    ║ Image          : ${params.image}
    ║ Python warnings: ${params.pythonwarnings}
    ║ Processors     : ${params.nproc}
    ║ Data directory : ${params.data_dir}
    ║ Run directory  : ${params.run_dir}
    ╚═══════════════════════════════════════════════════════════════
    """.stripIndent()

process versions_report {
    container params.image
    publishDir params.run_dir, mode: 'copy'
    
    output:
    path 'versions-report.txt'
    
    script:
    """
    versions report > versions-report.txt
    """
}

process preprocess_phix {
    container params.image
    publishDir "${params.run_dir}/phix", mode: 'copy'
    
    input:
    val ready
    
    output:
    path '**'
    
    script:
    """
    preprocess-phix \
        -s /data/sequences/HD13M_10K_R1.fastq \
        -o . \
        -p ${params.nproc}
    """
}

process presto_abseq {
    container params.image
    publishDir "${params.run_dir}/abseq", mode: 'copy'
    
    input:
    val ready
    
    output:
    path '**'
    
    script:
    """
    presto-abseq \
        -1 /data/sequences/HD13M_10K_R1.fastq \
        -2 /data/sequences/HD13M_10K_R2.fastq \
        -y /data/report.yaml \
        -n HD13M \
        -o . \
        -p ${params.nproc}
    """
}

process presto_clontech {
    container params.image
    publishDir "${params.run_dir}/clontech", mode: 'copy'
    
    input:
    val ready
    
    output:
    path '**'
    
    script:
    """
    presto-clontech \
        -1 /data/sequences/HD13M_10K_R1.fastq \
        -2 /data/sequences/HD13M_10K_R2.fastq \
        -j /usr/local/share/protocols/Universal/Human_IG_CRegion_RC.fasta \
        -r /usr/local/share/igblast/fasta/imgt_human_ig_v.fasta \
        -y /data/report.yaml \
        -n HD13M \
        -o . \
        -p ${params.nproc}
    """
}

process presto_clontech_umi {
    container params.image
    publishDir "${params.run_dir}/clontech-umi", mode: 'copy'
    
    input:
    val ready
    
    output:
    path '**'
    
    script:
    """
    presto-clontech-umi \
        -1 /data/sequences/HD13M_10K_R1.fastq \
        -2 /data/sequences/HD13M_10K_R2.fastq \
        -j /usr/local/share/protocols/Universal/Human_IG_CRegion_RC.fasta \
        -r /usr/local/share/igblast/fasta/imgt_human_ig_v.fasta \
        -n HD13M \
        -o . \
        -p ${params.nproc}
    """
}

process changeo_10x {
    container params.image
    publishDir "${params.run_dir}/10x", mode: 'copy'
    
    input:
    val ready
    
    output:
    path '**'
    
    script:
    """
    bash changeo-10x \
        -s /data/sequences/PBMC2B.fasta \
        -a /data/sequences/PBMC2B_annotations.csv \
        -x 0.15 \
        -n PBMC2B \
        -o . \
        -p ${params.nproc} \
        -z
    """
}

process changeo_igblast {
    container params.image
    publishDir "${params.run_dir}/changeo", mode: 'copy'
    
    input:
    val ready
    
    output:
    path '**', emit: results
    path 'combined_db-pass.tsv', emit: db_pass
    
    script:
    """
    changeo-igblast \
        -s /data/mixed-data/test_combined_repertoire.fasta \
        -n combined \
        -o . \
        -p ${params.nproc} \
        -z
    """
}

process changeo_clone {
    container params.image
    publishDir "${params.run_dir}/changeo", mode: 'copy'
    
    input:
    path db_pass
    
    output:
    path '*.tsv'
    
    script:
    """
    changeo-clone \
        -d ${db_pass} \
        -x 0.15 \
        -n combined \
        -o . \
        -p ${params.nproc} \
        -z
    """
}

process tigger_genotype {
    container params.image
    publishDir "${params.run_dir}/changeo", mode: 'copy'
    
    input:
    path db_pass
    
    output:
    path 'combined_igh_db-pass_parse-select.tsv'
    path 'combined_geno*'
    
    script:
    """
    # Subset to heavy chain
    ParseDb.py select \
        -d ${db_pass} \
        -f v_call \
        -u IGH \
        --regex \
        --outname combined_igh_db-pass \
        --outdir .
    
    # Run TIgGER genotyping
    tigger-genotype \
        -d combined_igh_db-pass_parse-select.tsv \
        -v v_call_genotyped \
        -x 2 \
        -y 2 \
        -u F \
        -n combined \
        -o . \
        -p ${params.nproc}
    """
}

process shazam_threshold {
    container params.image
    publishDir "${params.run_dir}/changeo", mode: 'copy'
    
    input:
    path db_pass
    
    output:
    path 'combined-*'
    
    script:
    """
    # Test density method
    shazam-threshold \
        -d ${db_pass} \
        -m density \
        -n combined-1 \
        -o . \
        -p ${params.nproc}
    
    # Test GMM method
    shazam-threshold \
        -d ${db_pass} \
        -m gmm \
        -n combined-2 \
        -o . \
        -p ${params.nproc}
    
    # Test with subsampling
    shazam-threshold \
        -d ${db_pass} \
        --subsample 100 \
        --repeats 2 \
        -n combined-3 \
        -o . \
        -p ${params.nproc}
    
    # Test none method
    shazam-threshold \
        -d ${db_pass} \
        -m none \
        -n combined-4 \
        -o . \
        -p ${params.nproc}
    """
}

workflow {
    // Run version report first
    versions_report()
    
    // Create a ready channel to trigger parallel processes
    ready = versions_report.out.collectFile()
    
    // Run all preprocessing and presto pipelines in parallel
    preprocess_phix(ready)
    presto_abseq(ready)
    presto_clontech(ready)
    presto_clontech_umi(ready)
    changeo_10x(ready)
    
    // Run changeo-igblast
    changeo_igblast(ready)
    
    // Run dependent processes that need the db_pass file
    changeo_clone(changeo_igblast.out.db_pass)
    tigger_genotype(changeo_igblast.out.db_pass)
    shazam_threshold(changeo_igblast.out.db_pass)
}

workflow.onComplete {
    println ""
    println "Pipeline execution summary"
    println "---------------------------"
    println "Completed at : ${workflow.complete}"
    println "Duration     : ${workflow.duration}"
    println "Success      : ${workflow.success}"
    println "Work directory: ${workflow.workDir}"
    println "Results saved to: ${params.run_dir}"
    println ""
}
