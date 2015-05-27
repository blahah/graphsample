# graphsample

Subsample FASTQ by sampling connected components of a de-Bruijn graph

Taking a subsample from a FASTQ can lead to poor coverage of some regions in the subsample, causing the subsample to have informational properties that are not representative of the full read set.

`graphsample` addresses this problem by building a [de-Bruijn graph](http://en.wikipedia.org/wiki/De_Bruijn_graph) from the reads, identifying all the [connected components](http://en.wikipedia.org/wiki/Connected_component_%28graph_theory%29), and randomly sampling those components. It outputs the reads that belong to the sampled components.

What's the point? Glad you asked! `graphsample` allows you to take a small subsample from a large set of reads, and use the subsample to optimise the parameters of any tools and algorithms you want to run on the full set.

## Compiling

```bash
$ git clone --recursive https://github.com/Blahah/graphsample.git
$ cd graphsample
$ make
```

## Running

```bash
$ bin/graphsample --help
```
