#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- Default Parameters ---
params.input_dir = "data/fastq_wgbs" // Directory where FASTQ files matching globs are located
params.samplesheet = "samples_wgbs.csv"
params.outdir = "results_wgbs"

// Trim Galore parameters
params.clip_r1 = 15
params.clip_r2 = 15
params.trim_galore_cores = 2 // Cores for Trim Galore itself (Cutadapt within it)

// Bismark parameters
params.bismark_genome_folder = "genome/" // Path to Bismark prepared genome
params.bismark_non_directional = true
params.bismark_aligner = "bowtie2" // or "hisat2"
params.bismark_aligner_cores = 4 // Cores for Bowtie2/HISAT2 (Bismark uses 1 for itself)

// Methylation extractor parameters
params.meth_extractor_no_overlap = true
params.meth_extractor_paired_end = true // For WGBS this is usually true
params.meth_extractor_cytosine_report = true
params.meth_extractor_bedgraph = true
params.meth_extractor_zero_based = true
params.meth_extractor_multicore = 2 // Cores for methylation extractor sub-processes

// QC parameters
params.skip_fastqc_raw = false
params.skip_multiqc = false

log.info """
         W G B S - P I P E L I N E (Bismark)
         =====================================
         Input Samplesheet: ${params.samplesheet}
         Input Directory  : ${params.input_dir}
         Output Directory : ${params.outdir}
         ---
         Bismark Genome   : ${params.bismark_genome_folder}
         Bismark Aligner  : ${params.bismark_aligner}
         Non-Directional  : ${params.bismark_non_directional}
         Trim R1 clip     : ${params.clip_r1}
         Trim R2 clip     : ${params.clip_r2}
         ---
         Run Raw FastQC   : ${!params.skip_fastqc_raw}
         Run MultiQC      : ${!params.skip_multiqc}
         """.stripIndent()

// --- Channel for input samples ---
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row ->
        def meta = [:]
        meta.id = row.sample
        def r1_glob = row.fastq_r1_glob
        def r2_glob = row.fastq_r2_glob

        def r1_files = file("${params.input_dir}/${r1_glob}").toList()
        def r2_files = file("${params.input_dir}/${r2_glob}").toList()

        if (r1_files.isEmpty() || r2_files.isEmpty()) {
            log.warn "WARNING: No R1 or R2 files found for sample ${meta.id} with globs:\n  R1: ${params.input_dir}/${r1_glob}\n  R2: ${params.input_dir}/${r2_glob}"
            return null
        }
        // Sort files to ensure deterministic cat order
        r1_files.sort()
        r2_files.sort()
        return tuple(meta, r1_files, r2_files)
    }
    .filter { it != null }
    .set { ch_raw_fastqs_pe }


// --- Workflow Definition ---
workflow {
    MERGE_FASTQS_PE(ch_raw_fastqs_pe)
    ch_merged_reads = MERGE_FASTQS_PE.out.merged_reads // tuple(meta, r1_merged, r2_merged)

    if (!params.skip_fastqc_raw) {
        FASTQC_RAW_PE(ch_merged_reads)
    }

    TRIM_GALORE_PE(ch_merged_reads)
    ch_trimmed_reads = TRIM_GALORE_PE.out.trimmed_reads // tuple(meta, r1_val, r2_val)

    BISMARK_ALIGN(
        ch_trimmed_reads,
        file(params.bismark_genome_folder)
    )
    ch_aligned_bam = BISMARK_ALIGN.out.bam // tuple(meta, bam, report)

    DEDUPLICATE_BISMARK(ch_aligned_bam.map{ meta, bam, report -> tuple(meta, bam) }) // only need bam
    ch_dedup_bam = DEDUPLICATE_BISMARK.out.dedup_bam // tuple(meta, dedup_bam, report)

    BISMARK_METHYLATION_EXTRACTOR(
        ch_dedup_bam.map{ meta, bam, report -> tuple(meta, bam) }, // only need bam
        file(params.bismark_genome_folder)
    )

    if (!params.skip_multiqc) {
        ch_multiqc_files = Channel.empty()
        if (!params.skip_fastqc_raw) {
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW_PE.out.zip.flatMap()) // flatMap if process emits multiple zips per item
        }
        ch_multiqc_files = ch_multiqc_files.mix(TRIM_GALORE_PE.out.fastqc_zip.flatMap())
        ch_multiqc_files = ch_multiqc_files.mix(TRIM_GALORE_PE.out.trim_report.map{ meta, r1_rep, r2_rep -> [r1_rep, r2_rep] }.flatMap())
        ch_multiqc_files = ch_multiqc_files.mix(BISMARK_ALIGN.out.report)
        ch_multiqc_files = ch_multiqc_files.mix(DEDUPLICATE_BISMARK.out.report)
        ch_multiqc_files = ch_multiqc_files.mix(BISMARK_METHYLATION_EXTRACTOR.out.reports.map { meta, mbias1, mbias2, splitting_report, nucl_report -> [mbias1, mbias2, splitting_report, nucl_report] }.flatMap())

        MULTIQC(ch_multiqc_files.collect().map{ it.unique() }) // Ensure unique files before collect
    }
}

// --- Process Definitions ---

process MERGE_FASTQS_PE {
    tag "${meta.id}"
    publishDir "${params.outdir}/01_merged_fastq", mode: 'copy', pattern: "*_merged.fastq.gz"

    input:
    tuple val(meta), path(r1_files), path(r2_files)

    output:
    tuple val(meta), path("${meta.id}_R1_merged.fastq.gz"), path("${meta.id}_R2_merged.fastq.gz"), emit: merged_reads

    script:
    """
    cat ${r1_files.join(' ')} > ${meta.id}_R1_merged.fastq.gz
    cat ${r2_files.join(' ')} > ${meta.id}_R2_merged.fastq.gz
    """
}

process FASTQC_RAW_PE {
    tag "${meta.id}"
    publishDir "${params.outdir}/02_fastqc_raw", mode: 'copy'
    conda "bioconda::fastqc=0.12.1"

    when: !params.skip_fastqc_raw

    input:
    tuple val(meta), path(r1_merged), path(r2_merged)

    output:
    tuple val(meta), path("*_fastqc.{zip,html}"), emit: zip_html // Using a glob to catch both
    path "*.zip", emit: zip
    path "*.html", emit: html

    script:
    """
    fastqc -o . $r1_merged $r2_merged
    """
}

process TRIM_GALORE_PE {
    tag "${meta.id}"
    publishDir "${params.outdir}/03_trimmed_fastq", mode: 'copy', pattern: "*_val_?.fq.gz"
    publishDir "${params.outdir}/03_trimmed_fastq/reports", mode: 'copy', pattern: ["*_trimming_report.txt", "*_fastqc.{zip,html}"]
    conda "bioconda::trim-galore=0.6.10 bioconda::cutadapt=4.1 bioconda::fastqc=0.12.1" // FastQC needed if --fastqc is used

    input:
    tuple val(meta), path(r1_merged), path(r2_merged)

    output:
    tuple val(meta), path("${r1_merged.baseName}_val_1.fq.gz"), path("${r2_merged.baseName}_val_2.fq.gz"), emit: trimmed_reads
    tuple val(meta), path("*_R1_merged_trimming_report.txt"), path("*_R2_merged_trimming_report.txt"), emit: trim_report
    tuple val(meta), path(["*_R1_merged_val_1_fastqc.zip", "*_R2_merged_val_2_fastqc.zip"]), emit: fastqc_zip // Correctly named outputs from trim_galore --fastqc
    tuple val(meta), path(["*_R1_merged_val_1_fastqc.html", "*_R2_merged_val_2_fastqc.html"]), emit: fastqc_html

    script:
    def clip_r1_opt = params.clip_r1 > 0 ? "--clip_R1 ${params.clip_r1}" : ""
    def clip_r2_opt = params.clip_r2 > 0 ? "--clip_R2 ${params.clip_r2}" : ""

    """
    trim_galore \\
        --paired \\
        --gzip \\
        $clip_r1_opt \\
        $clip_r2_opt \\
        --cores ${params.trim_galore_cores} \\
        --fastqc \\
        -o . \\
        $r1_merged \\
        $r2_merged
    """
}

process BISMARK_ALIGN {
    tag "${meta.id}"
    publishDir "${params.outdir}/04_bismark_alignment", mode: 'copy', pattern: ["*.bam", "*.{txt,html}"]
    // Bismark needs its chosen aligner (bowtie2 or hisat2) and samtools in the environment
    conda "bioconda::bismark=0.24.2 bioconda::samtools=1.17 ${params.bismark_aligner == 'bowtie2' ? 'bioconda::bowtie2=2.5.1' : 'bioconda::hisat2=2.2.1'}"

    input:
    tuple val(meta), path(r1_trimmed), path(r2_trimmed)
    path bismark_genome_dir

    output:
    tuple val(meta), path("${r1_trimmed.baseName}_pe.bam"), emit: bam
    tuple val(meta), path("${r1_trimmed.baseName}_PE_report.txt"), emit: report
    path "*_SE_report.txt", emit: se_reports, optional: true // if any single reads align
    path "*.html", emit: html_report, optional: true // Bismark summary report

    script:
    def non_directional_opt = params.bismark_non_directional ? "--non_directional" : ""
    int bismark_parallel_instances = Math.max(1, params.bismark_aligner_cores / 2 ) // e.g. 4 cores -> 2 parallel instances
    int threads_per_instance = Math.max(1, params.bismark_aligner_cores / bismark_parallel_instances)
    def aligner_core_opts = "--${params.bismark_aligner}_options \\\"--threads $threads_per_instance\\\" --parallel $bismark_parallel_instances"
    if (params.bismark_aligner_cores < 2) aligner_core_opts = "" // single threaded if not enough cores

    """
    bismark \\
        --genome_folder $bismark_genome_dir \\
        --${params.bismark_aligner} \\
        $aligner_core_opts \\
        $non_directional_opt \\
        --bam \\
        -o . \\
        -1 $r1_trimmed \\
        -2 $r2_trimmed

    # Bismark might create an HTML summary report with a different name
    if ls *.html 1> /dev/null 2>&1; then
      mv *.html ${meta.id}_bismark_summary.html
    fi
    """
}

process DEDUPLICATE_BISMARK {
    tag "${meta.id}"
    publishDir "${params.outdir}/05_deduplicated_bam", mode: 'copy', pattern: ["*.bam", "*.txt"]
    conda "bioconda::bismark=0.24.2" // Samtools is usually bundled or a dependency

    input:
    tuple val(meta), path(aligned_bam)

    output:
    tuple val(meta), path("${aligned_bam.baseName}.deduplicated.bam"), emit: dedup_bam
    tuple val(meta), path("${aligned_bam.baseName}.deduplication_report.txt"), emit: report

    script:
    """
    deduplicate_bismark \\
        --paired \\
        --bam \\
        -o . \\
        $aligned_bam
    """
}

process BISMARK_METHYLATION_EXTRACTOR {
    tag "${meta.id}"
    publishDir "${params.outdir}/06_methylation_calls", mode: 'copy'
    conda "bioconda::bismark=0.24.2" // Samtools is usually bundled or a dependency

    input:
    tuple val(meta), path(deduplicated_bam)
    path bismark_genome_dir

    output:
    tuple val(meta), path("*.bedGraph.gz"), optional: true, emit: bedgraph
    tuple val(meta), path("*.cov.gz"), optional: true, emit: covgz
    tuple val(meta), path("*M-bias_R1.txt"), path("*M-bias_R2.txt"), path("*PE_report.txt"), path("*splitting_report.txt"), emit: reports // Multiple report files
    path "*", emit: all_files // Catch all other output files

    script:
    def paired_opt = params.meth_extractor_paired_end ? "--paired-end" : "--single-end"
    def no_overlap_opt = params.meth_extractor_no_overlap && params.meth_extractor_paired_end ? "--no_overlap" : ""
    def cytosine_report_opt = params.meth_extractor_cytosine_report ? "--cytosine_report" : ""
    def bedgraph_opt = params.meth_extractor_bedgraph ? "--bedGraph" : ""
    def zero_based_opt = params.meth_extractor_zero_based && params.meth_extractor_bedgraph ? "--zero_based" : "" // only with bedGraph
    def multicore_opt = params.meth_extractor_multicore > 1 ? "--multicore ${params.meth_extractor_multicore}" : ""

    """
    bismark_methylation_extractor \\
        $paired_opt \\
        $no_overlap_opt \\
        $cytosine_report_opt \\
        $bedgraph_opt \\
        $zero_based_opt \\
        $multicore_opt \\
        --gzip \\
        --report \\
        --genome_folder $bismark_genome_dir \\
        -o . \\
        $deduplicated_bam
    """
}


process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    conda "bioconda::multiqc=1.19"

    when: !params.skip_multiqc

    input:
    path ('*') // Collects all files passed to it

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc -f .
    """
}
