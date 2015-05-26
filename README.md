# graphsample

Subsample FASTQ by sampling connected components of a de-Bruijn graph

Taking a subsample from a FASTQ can lead to poor coverage of some regions in the subsample, causing the subsample to have informational properties that are not representative of the full read set.

`graphsample` addresses this problem by building a de-Bruijn graph from the reads, identifying all the components (connected sub-graphs), and randomly sampling those components. It outputs the reads that belong to the sampled components.

What's the point? This allows you to take small subsamples from large sets of reads, and use the subsample to optimise the parameters of any tools and algorithms you want to run on the full set.
