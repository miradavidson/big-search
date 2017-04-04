# Big Search Case Study

## Background

Biology is a data intensive science. The simplest representations of the
human genome require roughly 3 GB of data. Fully characterizing a single
cell requires over a terabyte. Thus before a Data Scientist can apply the
power of modern machine learning techniques, one must first learn to write
code capable of handling such volumes of data.

While Python is generally a strong language choice for Data Science, it
demands a higher level of skill with the language, its idioms, and its
libraries is required to work with real-world biological data sets.

This case study was inspired by a 2008 blog post by Python's very own
benevolent dictator for life, Guido van Rossum [available
here](http://neopythonic.blogspot.co.uk/2008/10/sorting-million-32-bit-integers-in-2mb.html)
along with our own work on large genomic datasets.


## Objective

Your objective is to solve the k-mismatch problem in CPython and prepare a 30 min talk on your solution.
A formal statement of the k-mismatch problem problem is:

```
Given a set of strings selected from the alphabet ACGT which together form a Text, T, and a query string, P, find *all* occurrences of P
 in T (substrings of T) with up to k mismatches.
```

In this case, T is the reference primary assembly of the human genome (~3GB) and P is a query pattern between 4 and 150 characters in length.
In real world scenarios this pattern may be a Feature one wishes to detect, a sequencing read, or the target sequence of a DNA modifying enzyme. The internet abounds with O(m log n) solutions to this problem.
Knuth-Morris-Pratt and Rabin-Karp are both popular if indexing (pre-processing) is not permitted.
Inverted Indices, Suffix Arrays, and FM Indices are all popular data structures where the text is pre-processed.

The underlying goal is to demonstrate your problem solving skills and knowledge of Python, its idioms, and its libraries, *not* to develop a wholly novel algorithm.
A naive solution in C/C++ can solve this problem in about 1 second without any preprocessing of the text, but for this case study you should make your own implementation in Python.
This is a deliberate choice as this problem forces the developer to contend with several quirks of CPython's memory management model, IO functionality, List and String data structures, and sort and find standard library functions.
You may write your own CPython C extensions, use external libraries such as numpy, or tools such as Cython or CFFI, but you must the extend the setup.py script to handle compilation of any such extensions.
However, you may *not* simply call a 3rd party executable such as bwa, bowtie, nblast, grep, etc. but you may statically link to them if you automate a portable build process.

The problem is deliberately presented in an open ended way.
It may be solved either through profiling and optimization, algorithmic sophistication, or effective use of 3rd party libraries.
It is deliberately designed to expose some of Python's most interesting quirks and we hope you learn from it.
Since the problem permits both very simple and very complex solutions and your creativity is encouraged.
Here is one breakdown of the problem from which you may want to pick an strategy depending on time constraints, interest, and ability level.

### Hints

1. First solve the problem with a slow and simple searcher for one record
   and without preprocessing.
1. Solve only for k = 0, one relatively long p, and a single member string of the text (good if you're stuck and for checking).
1. Assume the text and pattern are restricted to the four letter alphabet 'ACGT'
1. Don't consider insertions, deletions, gaps, or wildcard patterns (simplifies from O(len(P)^2*len(T)) to O(len(P)*log(len(T)))).
1. Search the text using built in library functions (find, re.match).
1. For a given k, a substring of the text of length r must match the pattern exactly, where = (len(P)/(k+1)).


### Example Strategy 1
1. Pre-process the text into an inverted index (keyword, n-gram, k-mer index) constructed out of built-in Python structures, numpy arrays, or third-party libraries
1. Process both the forward and reverse strands of the text, but write the "postings list" to disk, using one file per word.
4. Use a "seed and extend" search strategy: find an exactly matching subsequence of the pattern and then extend it
1. For each chuck of the pattern, load the relevant posting list as the seeds. Numpy arrays may be helpful here.
1. For each seed, slice out the relevant chunk of the text and see if the rest matches.
5. Use multiprocessing to process multiple seed matches in parallel.
NB: This is similar to to how Elastic Search and Google work. Its very fast, but very memory intensive if the postings lists are stored in RAM. Work well for small values of k and len(P).

### Example Strategy 2
1. Use a fast implementation of an index-free string search algorithm  KMP, RK, or Smith-Waterman
1. Memory map the text file
1. Use multiprocessing to search each DNA strand in parallel
1. Merge the results together.
NB: This approach is used by cli tools like [the silver searcher](https://github.com/ggreer/the_silver_searcher) which takes a few seconds. Work in production for large values of k.

### Example Strategy 3
1. Enumerate the suffices of the Text as a set of (record.id, sequence_offset) tuples. Numpy record arrays are helpful here.
1. Use the numpy sort routines to sort these suffices lexographically to make a suffix array (or use a 3rd party library).
1. Initalize a list of match tuples marking the interval of the suffix array that correspond to matching segments of the array and the total number of mismatches.
1. For each character in the pattern, for each candidate match interval, perform a binary search to find the intervals which contain the next matching character.
1. Alternatively increment the mismatch counter for the candidate match..
1. Remove candidate matches that exceed k mismatches.
1. Stop when the pattern is fully traversed or when no candidate matches remain.
NB: This is the standard binary search method for a suffix array and example code is easily found online. The key idea here is to use numpy instead of built in Python lists to reduce the memory footprint.


### Advanced Suggestions
1. Consider a high-performance index data-structure such as a suffix array, FM-index, or efficient key-value store.
1. Consider impact of time spent writing results to output, OS level IO caching, and converting data types.
1. Consider start up time spent loading the text and index vs. run-time spent finding and reporting results
1. Consider reducing the  memory footprint by storing the text and index on disk, memory-mapping them, and only loading relevant subsections
1. Consider the non-random, highly-repetitive nature of the pattern and text, alternative encodings, and compression opportunities.
1. Consider a 2bit binary encoding: A => 00, C => 01, G => 10, T => 11 and use of bitwise operations for comparison
1. Does the performance vary for different values of k and P?


## Evaluation Criteria
In evaluating solutions the following criteria will be considered in order of importance:

1. Does the solution solve the problem and do the technical decisions make sense?
1. Does a green test suite demonstrate the solution is correct?
1. Does the solution build properly? Does setup.py or an install.sh script cover all steps in building the code?
1. Does the code use a consistent coding style and has it been linted and follow PEP8 conventions? Are Python idioms followed?
1. Does the profiler demonstrate the code is fast? A production-grade solution would < 1 second; naive python solutions take a few minutes. Is the cost of this speed in system resources, complexity, dependency burden and readability justified?
1. Does the code communicate clearly through its style, comments, docstrings, test suite, and/or accompanying docs?

## Guidelines

1. Your work should be your own, but you may use external libraries, resources, and friends. You may ask questions.
1. You may use as much time as you like, but bear in mind we consider applications on a rolling basis.
1. You should solve the problem for a least one example pattern (see test suite) and one value of k > 0. The optimal solution is a function of k, len(P), len(T), and len(alphabet), but the depth to which you explore this is discretionary.
1. You should consider the speed, memory footprint, and complexity of your solution and justify your decisions.
1. You should assume T is immutable and known in advance (reference genomes are updated only a few times a year). The pattern P is only known at run time. You may preprocess T but not P, and this preprocessing can be done "offline".
1. Profiling should include loading T (or its preprocessed version), reading P from stdin, executing the search, and writing both the location and substring results to a file. Beware of time spent on serialization and marshalling!
1. Your solution should find *all* matches meeting the criteria. There is a very similar problem of finding the *optimal* matching substring or *first* matching substring that have different solutions, but don't solve these problems.
1. The PySAM (pysam) library is included in the requirements.txt file. This a Python wrapper around the htslib C library used by the "pros" and includes a parser for the "fasta" file format used to store the genome sequence along with a wide range of other tools.
1. Indexing the genome AND searching the genome are two separate but interrelated problems; focus on the latter first.

## Data files

Human reference genome:
ftp://ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

Example output for two sequences in .sam format (generated using bowtie aligner).

