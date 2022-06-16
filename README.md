# DNA clustering project

This is an implementation of DNA clustering that improves on previous work done by Raschtian et al. in their paper:
["Clustering billions of reads for DNA data storage"](https://papers.nips.cc/paper/6928-clustering-billions-of-reads-for-dna-data-storage)

Please see the research report, for full details and expected performance: in `research_report.pdf`

## Usage
To run the clustering algorithm on your own data , bring an input file of simulated DNA strands,  in the form ${A,C,G,T}$.  
Place the input into the folder `files/minion_idt`, and then call the clustering algorithm from `main.py` with the correct input file path.

Placeholder input files exist with files of 3000 and 6000 strands, for example:            
`files\minion_idt\3000 strands in size 150 with x2 errors and cluster avg of 40\evyat files\evyat00.txt`.     
Make sure that your input is in the same format as this file.   
There are input files both with and without indices in the strands metadata.

## Important Functions
__`algo_clustering_to_file`__ :(in Clustering.py), takes the input file, runs the Microsoft clustering algorithm, and saves the result to a file.


__`test_stats`__: Takes the output file form the Microsoft algorithm, and runs the `handle_unions` and `handle_singletons` methods from `clustering.py`.
It then produces statistics from the result files and logs the results.  


__`test_times`__: Receives functions to apply to the input files and runs them, logging the runtimes of the different functions. Helpful for evaluating the runtime complexity of our algorithm.



 
