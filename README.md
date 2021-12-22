# Velocyto/scVelo.py for a set of bulk RNASeq samples

A Nextflow pipeline based on an existing 
[tutorial](https://github.com/praneet1988/Inferring-and-Visualizing-RNA-Velocity-in-Bulk-RNA-SEQ) 
by Praneet Chaturvedi and Abhijit Badve.

# Instructions

There's no need to download this repository, it can be downloaded and run directly from 
the command line as such:

   nextflow run slagtermaarten/bulk_velo-nf \
     --bam_dir /path/to/dir/containing/bam_files \
     --bam_glob '*.bam' \
     --cpu 8 \
     --outdir 'results' \
     --chroms '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y' \
     --mask /path/to/intron_location_mask \
     --gtf /path/to/gene_locations \
     --N_pc 10

The parameter names should mostly be self-explanatory. `bam_glob` allows you to sub-select 
bam files in the `bam_dir` using regex. `N_pc` sets the number of principal components on 
which scVelo is run.

Of course, the intron masking file and gtf files should be compatible with the reference 
genome to which the RNASeq data was aligned. For more information, please see the original 
[tutorial](https://github.com/praneet1988/Inferring-and-Visualizing-RNA-Velocity-in-Bulk-RNA-SEQ) 
on which this was based.
