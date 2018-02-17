# Requirements

**NOTE: DMFinder can only analyze PAIRED-END data produced with Illumina platforms (i.e. read pairs with forward-reverse orientation).**

Please make sure that chromosome names in the input BAM file are in the following format: chr1, chr2, chr3...

DMFinder has three external dependencies: SamTools, Graph, and GraphEasy. The link to each is provided here:

* SamTools - http://samtools.sourceforge.net. **Be sure to add the SamTools directory to your PATH variable**.
* Graph - http://search.cpan.org/~jhi/Graph-0.9704/lib/Graph.pod
* GraphEasy - http://search.cpan.org/~tels/Graph-Easy/lib/Graph/Easy.pm	

Requiremens for AmpliconFinder see in tools/ampiconfinder/README.md.

The Graph package is provides the infrastructure for implementing the graph algorithms used by DMFinder. The GraphEasy package allows for visualization of the predicted double minutes (written to a file called "graph.dot").


# Installation

```
git clone --recursive https://github.com/rmarduga/DMFinder.git
```

Once the dependencies have been installed and the files have been extracted, DMFinder can be executed with the following command:

perl dm_find.pl
	
NOTE: The dm_find_core.pl program contains most of the logic to run the program, but it is not recommended that you execute this file directly. 




# Running DMFinder

```
dm_find.pl [OPTIONS] --input_bam INPUT_BAM --sv SV_FILE --cn CN_FILE
                     --report_file REPORT_FILE  --graph_file GRAPH_FILE

    INPUT_BAM           indexed and sorted bam file
    SV_FILE             structural variant coordinates in VCF format
    CN_FILE             copy number variant coordinates in BED format
    REPORT_FILE         path to report file with found double minutes in csv
                        format
    GRAPH_FILE          path to graph file with found double minutes in dot
                        format
OPTIONS

  --min_cyclic INT      Minimum number of amplicons in cyclic chain to predict
                        as a double minute [2]
  --min_non_cyclic INT  Minimum number of amplicons in a non-cyclic chain to
                        predict as a DM [2]
  --cutoff INT          Cutoff coefficient of mapped read pair standard
                        deviation [4]
  --mean INT            Mean mapped distance between read pairs [400]
  --stdev INT           Standard deviation of mapped distances between
                        reads [80]
  --minqual INT         Minimum read mapping quality to estimate mapping
                        coverage [30]
  --window INT          Maximum number of base pairs to define proximity
                        between copy number breakpoint and structural variant
                        breakpoint [2*(cutoff*stdev+mean)]
  --split_amplicons     Allow splitting amplicons from copy number file if one
                        of the ends of a structural variant edge falls in
                        between of the amplicon region. (experimental)
  --no_avg_cov_check    Predict DMs regardless of average mapping coverage
                        among its amplicons.
  --verbose             Make the operation more talkative
```

* **input_bam** -	This should be a coordinate-sorted bam file with an associated index file. The index filename should be your bam filename with ".bai" affixed to the end.
* **sv** - This is a VCF format file that contains the structural variant breakpoints. **You must convert your structural variant breakpoint prediction files into VCF format to be read by DMFinder**.
* **cn** - This is a BED format file that contains coordinates corresponding to copy number amplified segments. **DMFinder assumes these regions correspond to segments of copy number gain**.
* **min_cyclic** - For circular (cyclic) sub-graphs, this is the minimum number of amplicons (vertices) that must exist for predicting a double minute. 
* **min_non_cyclic** - For NON-circular sub-graphs, this is the minimum number of amplicons (vertices) that must exist for predicting double minute. Per algorithm 3, this is the smallest allowable weakly connected component that is allowed to predict a double minute.
* **cutoff** - This is the coefficient of standard deviations of the mean for measuring the window size. Not used if the user specifies a window length.
* **mean** - This is the mean mapped distance between read pairs (insert size) of the input bam file. Not used if the user specifies a window length.
* **stdev** - Standard deviation of mapped distances of read pairs of the input bam file. Not used if the user specifies a window length.
* **minqual** - Minimum read mapping quality for computing amplicon read coverage.
* **window** - Maximum distance between amplicons (vertices) to connect by an edge. 	 

# Output format

DMFinder 


DMFinder produces two files:
1. GRAPH_FILE which contains the graphs of each predicted double minute. It can be visualized using GraphViz (http://www.graphviz.org).

2. REPORT_FILE which contains the report in CSV format of predicted double minutes. The first row in the file contains the columns' headers. 
   * **DM_index** index of form DM<INDEX> (for instance, DM1). All amplicons with the same index belong to the same DM. 
   * **chromosome**
   * **amplicon start**
   * **amplicon end**
   * **type** one of 2 types: *complete* double minute when the method correctly identified all DM amplicons connected in a cycle of structural variant (SV) breakpoints (a hallmark of DMs). Else, if the DM was partially identified, meaning that at least *m* (specified with --min_non_cyclic parameter) amplicons were discovered that were linearly connected by SV breakpoints , then we say that an *incomplete* double minute was discovered.

All amplicons of *complete* DM, whose rows appeared adjacently in REPORT_FILE are connected to each other. Additionally, the first and the last amplicons of the double minutes of *compleate* type are also connected to each other.


# Publication
Matthew Hayes and Jing Li. "An integrative framework for the identification of double minute chromosomes using next generation sequencing data". BMC Genetics, 2015, 16(Suppl 2):S1.
