This repository contains my implementations of various algorithms related to the [Generalized Time-Reversible (GTR) model](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#GTR:_Generalised_time-reversible_.28Tavar.C3.A9_1986.29.5B9.5D).

* **[generateSeq](generateSeq.py): Generate sequences on a phylogenetic tree using the GTR model**
    * This program performs a preorder traversal on the input tree and uses the GTR model to simulate sequences down the tree
    * Usage: `generateSeq.py [-h] -t TREE -p GTRPARAMS [-i ROOTSEQ] [-r ROOTSEQLEN] -o OUT`
	    * `-h`: Show help message
    	* `-t TREE`: Specify `TREE` (a Newick file) to be the input phylogenetic tree
	    * `-p GTRPARAMS`: Specify `GTRPARAMS` (a text file in the format specified in the homework instructions) to be the input GTR parameters
    	* `-i ROOTSEQ`: Specify `ROOTSEQ` (a FASTA file) to be the sequence of the root node of `TREE`
	    * `-r ROOTSEQLEN`: If `-i ROOTSEQ` is not specified, generate a random sequence of length `ROOTSEQLEN` using the stationary vector in `GTRPARAMS` to be the sequence of the root
    	    * If both `-i ROOTSEQ` and `-r ROOTSEQLEN` are specified, `-i ROOTSEQ` will have preference (i.e., `-r ROOTSEQLEN` will be ignored)
	    * `-o OUT`: Specify `OUT` to be the output file
    	    * Note that the file `OUT` will be overwritten if it already exists
* **[estimateGTRPairs](estimateGTRPairs.py): Estimate GTR parameters from pairs of DNA sequences**
	* Given a pair of DNA sequences and their distance, this program uses the "Sequential Least Squares Programming" (SLSQP) method to estimate maximum-likelihood GTR parameters
    * Usage: `estimateGTRPairs.py [-h] -f BRANCHLENGTH -d SEQS -o OUT`
	    * `-h`: Show help message
    	* `-f BRANCHLENGTH`: Specify `BRANCHLENGTH` to be the branch length (i.e., distance) between the two sequences
	    * `-d SEQS`: Specify `SEQS` (a FASTA file) to be the pair of sequences from which to estimate GTR parameters
    	* `-o OUT`: Specify `OUT` to be the output file
        	* Note that the file `OUT` will be overwritten if it already exists
* **[computeML](computeML.py): Compute the likelihood of a tree given GTR parameters**
    * Given a phylogenetic tree, sequences, and GTR parameters, this program uses the Felsenstein tree-pruning algorithm to compute the likelihood of the tree
    * Usage: `computeML.py [-h] -t TREE -p GTRPARAMS -d SEQS -o OUT`
	    * `-h`: Show help message
    	* `-t TREE`: Specify `TREE` (a Newick file) to be the input phylogenetic tree
	    * `-p GTRPARAMS`: Specify `GTRPARAMS` (a text file in the format specified in the homework instructions) to be the input GTR parameters
    	* `-d SEQS`: Specify `SEQS` (a FASTA file) to be the set of sequences for the leaves of `TREE`
	    * `-o OUT`: Specify `OUT` to be the output file
    	    * Note that the file `OUT` will be overwritten if it already exists
* **[estimateGTR](estimateGTR.py): Estimate GTR parameters on multiple DNA sequences**
    * Given a phylogenetic tree and sequences, this program uses the Felsenstein tree-pruning algorithm and the "Sequential Least Squares Programming" (SLSQP) method to estimate maximum-likelihood GTR parameters
    * Usage: `estimateGTR.py [-h] -t TREE -d SEQS [-i MAXIT] -o OUT`
	    * `-h`: Show help message
    	* `-t TREE`: Specify `TREE` (a Newick file) to be the input phylogenetic tree
	    * `-d SEQS`: Specify `SEQS` (a FASTA file) to be the set of sequences for the leaves of `TREE`
    	* `-i MAXIT`: Specify `MAXIT` (a positive integer) to be the maximum number of iterations for the optimization
	    * `-o OUT`: Specify `OUT` to be the output file
    	    * Note that the file `OUT` will be overwritten if it already exists

REQUIREMENTS
===
* [Python 3](https://www.python.org/downloads/) (NOT PYTHON 2!)
* [DendroPy](http://www.dendropy.org/) `pip3 install dendropy`
* [NumPy](http://www.numpy.org/) `pip3 install numpy`
* [SciPy](http://scipy.org/) `pip3 install scipy`
