# Variant calling pipeline

basic usage:

`vg_giraffe_pipeline.sh [name_of_ref_graph] [fastq_1] [fastq_2]`

Any line with [tk] needs to be filled in with the appropriate hard-coded directory on the final server/compute node.
It should read `name_of_ref_graph` and find the indices we've already built from that graph (with `vg gbwt` etc) that start with that name.
For example, "pggb_100.mod" (we had to "mod" the graph to have no nodes >1000 bp with `vg mod -X`).
