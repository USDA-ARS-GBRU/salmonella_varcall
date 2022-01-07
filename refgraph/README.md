## Pangenome Reference Graphs

This folder contains the custom *S. enterica* reference graph files and should eventually also include the script used to create them.

The graph files with "pggb_100" in their names were constructed using 100 complete, contiguous *S. enterica* genomes available in the NCBI assembly database. The genomes were selected from nearly 600 complete, circular, contiguous, high-quality genomes whose similarity was calculated using MASH. The MASH similaries were subsequently used for clustering into 100 clusters, from which 100 representative genomes were selected, 1 from each cluster.

Since the graph files are large (i.e. approaching 1 Gbyte in size), they are stored using the Git Large File Storage (LFS) system.

The graph files are used in the graph-based variant calling approach that uses vg tools.
