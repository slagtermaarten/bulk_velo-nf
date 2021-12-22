#!/usr/bin/env nextflow

/*
===============================================================================
                           BULK VELOCYTO PIPELINE
===============================================================================
*/


params.bam_dir = "${baseDir}/test_data"
params.input_bam_glob = '*.bam'
params.outdir = 'results'
params.cpu = 8
params.chroms = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X'
params.mask = "${baseDir}/test_data/mask.gtf"
params.gtf = "${baseDir}/test_data/genes.gtf"
params.N_pc = 5

bam_glob = params.bam_dir + '/' + params.input_bam_glob
i_bams = Channel.fromPath(bam_glob)
  .map { tuple( it.baseName, it ) }

// Add faux UMI information to reads to emulate SC data
process add_UMI {
  maxForks params.cpu
  conda "${baseDir}/env/samtools.yml"
  input:
    tuple val(sample_id), file('i.bam') from i_bams
  output:
    tuple val(sample_id), file('i.sam') into i_sams
  shell:
    """
    #!/usr/bin/env python

    import simplesam
    from simplesam import Reader, Writer

    with Reader(open('i.bam')) as in_bam:
      with Writer(open('i.sam', 'w'), in_bam.header) as out_sam:
        for read in in_bam:
          read["UB"] = read.qname.split(":")[2] # add the umi tag
          read["CB"] = read.qname.split(":")[1] # add the barcode tag
          out_sam.write(read)
    """
}

// Convert SAM to BAM and sort
process binarize_add_UMI {
  conda "${baseDir}/env/samtools.yml"
  maxForks params.cpu
  input:
    tuple val(sample_id), file('i.sam') from i_sams
  output:
    tuple val(sample_id), file('umi.bam') into umi_bams
  shell:
    """
    samtools view -b -o umi.bam i.sam
    samtools sort umi.bam -o umi.bam
    """
}

// Subset to chromosomes of interest
process subset_bam {
  conda "${baseDir}/env/samtools.yml"
  maxForks params.cpu
  input:
    tuple val(sample_id), file('umi.bam') from umi_bams
  output:
    tuple val(sample_id), file('subset.bam') into subset_bams
  shell:
    """
    samtools index umi.bam
    samtools view -F 4 -b umi.bam !{params.chroms} > subset.bam
    samtools sort subset.bam -o subset.bam
    """
}

// Run Velocyto
process run_velocyto {
  conda "${baseDir}/env/velocyto.yml"
  maxForks params.cpu
  input:
    tuple val(sample_id), file("subset.bam") from subset_bams
  output:
    tuple val(sample_id), file("${sample_id}.loom") into velo_ch
  shell:
    """
    samtools index subset.bam
    velocyto run -U -m !{params.mask} subset.bam !{params.gtf} \
      --sampleid !{sample_id} --outputfolder ./
    """
}

all_looms = velo_ch
  .toSortedList( { a, b -> a[0] <=> b[0] } ) // Sort by sample names
  .get()
  .collect({ it[1] }) // Restrict to file names (second item of each tuple)
  .join(' ')


process combine_velocyto {
  conda "${baseDir}/env/velocyto.yml"
  input:
    val(all_looms)
  output:
    file('combined.loom') into comb_ch, comb_ch_dup
  shell:
    """
    #!/usr/bin/env python
    import loompy, os

    file_list = '!{all_looms}'.split(' ')
    comb = loompy.combine(file_list, 'combined.loom',
      key="Accession")
    """
}

// loomR isn't working for me and I circumvent its use this way
process loom_to_scv {
  conda "${baseDir}/env/velocyto.yml"
  publishDir "${params.outdir}", mode: 'copy'

  input:
    file('combined.loom') from comb_ch_dup
  output:
    file('ambiguous.csv')
    file('matrix.csv')
    file('spliced.csv')
    file('unspliced.csv')
  shell:
    """
    python ${baseDir}/scripts/loom_to_csv.py combined.loom
    touch ambiguous.csv
    """
}

process sc_velo {
  conda "${baseDir}/env/velocyto.yml"
  publishDir "${params.outdir}", mode: 'copy'

  input:
    file('combined.loom') from comb_ch
  output:
    file('scVelo_analysis.loom')
  shell:
    """
    #!/usr/bin/env python

    import scvelo as scv
    import numpy as np
    import scanpy as sc

    adata = scv.read('combined.loom', cache=True)
    adata.var_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=0)
    sc.pp.filter_genes(adata, min_cells=0)

    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    N_pc = min(!{params.N_pc}, adata.obs.shape[0], adata.var.shape[0])
    scv.pp.moments(adata, n_pcs=N_pc)
    scv.tl.velocity(adata, mode='stochastic')
    scv.tl.velocity_graph(adata)

    adata.write('scVelo_analysis.loom')
    """
}

// vim: filetype=nextflow et ts=2 sw=2
