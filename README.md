# big-search
Case Study in Python searching datasets larger than RAM. This describes my implementation.

## Structure
### Input
The program takes genome input in FASTA file format using the `pysam` library, which helps speed-up the throughput of the sequencing data with the use of the C [hstlib](https://github.com/samtools/htslib) library.

At this point, I considered using Python's `mmap` function, which memory-maps files. However, it was unclear how this could interface with `pysam`. Further research indicated that although memory-mapping would reduce the physical RAM usage, the memory _address_ overhead was still significant and potentially a problem in a 32-bit addressing environment.

The genome is then split into an arbitrary number of (overlapping) sections, each handled by a separate process using Python's `multiprocessing` module. In my initial testing, using 4 separate processes showed a speedup of approximately 2x when run on a quad-core machine. Any more than 4 processes had little added effect.

The query string and maximum _k_ mismatches are taken as input from the user. This could easily be modified to be a command-line argument.

### Searching
After some research into string-search algorithms, I decided to use the [Knuth-Morris-Pratt](https://en.wikipedia.org/wiki/Knuth%E2%80%93Morris%E2%80%93Pratt_algorithm) (KMP) algorithm for its excellent worst-case _O(n)_ complexity for data of length _n_. This worst-case performance is particularly important given the limited base alphabet and the fact that the data is not necessarily completely random. The sacrifice in generating the LPS (Longest Prefix Suffix) array (which is _O(l)_ for a query of length _l_) is negligible given that the query length is capped at 150.

KMP doesn't easily allow for the introduction of mismatches given the fixed LPS array which is precalculated. I modified the algorithm to allow for any number of mismatches, which resulted in obvious but not unreasonable slowing of the run time. However, the modification does not take into account the fact that the LPS isn't the same when a mismatch occurs. This could potentially lead to the occasional missed match, but avoids the enormous slowdown that would result from either recalculating the LPS array for each mismatch, or returning to the point where the first mismatch occurred (and the LPS array is guaranteed to be correct).

When a match is found, it's added to the list `all_alignments` which is a global list managed by the `multiprocessing.Manager`. This saves having to combine the lists produced by each process at the end, and allows each process to complete without having to wait for every process to finish.

### Output
Once all processes have completed their search, the all_alignments list is saved to a file with one index per line, along with the corresponding match at that index.

## Performance
Preliminary testing was performed on a very small dataset. Once the initial functionality had been implemented, I started testing on the first chromosome, to give some idea of how long it would take on the full genome. The table below shows some experimental evidence of the performance of the program running on a MacBook with a quad-core processor and X GB of RAM.

| Query String      | Mismatches           | Runtime  |
| ------------- |:-------------:| -----:|
| insert-here      | insert-here | insert-here |
| insert-here      | insert-here      |   insert-here |
| insert-here | insert-here      |    insert-here |

The following table shows the runtimes recorded when the entire human genome was considered, searching each chromosome one by one.