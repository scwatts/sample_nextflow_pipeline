#!/usr/bin/env nextflow
// Enable DSL2
nextflow.enable.dsl=2


// Create channel from reads on disk
ch_reads = Channel.fromFilePairs('data/*_{1,2}.fastq.gz', flat: true)
// Set output directory
output_dir = 'output/'


process read_assembly {
  // Directive to copy FASTA assembly to the output directory
  publishDir path: {"${output_dir}"}, mode: 'copy'

  input:
  // Input is ch_reads
  // The isolate_id variable is provided so we can informatively name output files
  tuple val(isolate_id), path(reads_fwd), path(reads_rev)

  output:
  // Here we include isolate_id in the output so that it is also available to downstream processes
  // The '*fasta' glob captures skesa output and exposes it to the publishDir directive
  tuple val(isolate_id), path('*.fasta')

  script:
  // The skesa output is named using the input isolate_id variable
  """
  skesa --reads ${reads_fwd},${reads_rev} --cores 1 > ${isolate_id}.fasta
  """
}


process assembly_stats {
  // Directive to copy assembly stats output to the output directory
  publishDir path:{"${output_dir}"}, mode: 'copy'

  input:
  // Input is output of the read_assembly process
  // Again, the isolate_id is provided to name the output file
  tuple val(isolate_id), path(assembly_fasta)

  output:
  // The output file is defined here so that it is captured by publishDir
  path('*_stats.tsv')

  script:
  // The command output is uniquely named with the isolate_id
  """
  assembly_stats.py -a ${assembly_fasta} > ${isolate_id}_stats.tsv
  """
}


workflow {
  // Here we define how channels connect processes
  // ch_reads (channel) -> read_assembly (process) -> ch_assemblies (channel) -> assembly_stats (process)
  ch_assemblies = read_assembly(ch_reads)
  assembly_stats(ch_assemblies)
}
