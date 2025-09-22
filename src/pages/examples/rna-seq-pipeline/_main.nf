#!/usr/bin/env nextflow

/*
 * Proof of concept of a RNAseq pipeline implemented with Nextflow
 */

params.reads = "$baseDir/data/ggal/ggal_gut_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = "results"
params.multiqc = "$baseDir/multiqc"

// import modules
include { RNASEQ } from './modules/rnaseq'
include { MULTIQC } from './modules/multiqc'

workflow {
log.info """\
  R N A S E Q - N F   P I P E L I N E
  ===================================
  transcriptome: ${params.transcriptome}
  reads        : ${params.reads}
  outdir       : ${params.outdir}
  """

  read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
  RNASEQ( params.transcriptome, read_pairs_ch )
  MULTIQC( RNASEQ.out, params.multiqc )
}