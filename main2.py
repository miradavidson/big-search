import pysam
import timeit
from multiprocessing import Process, Manager
from time import localtime, strftime

sequence_filename = "Homo_sapiens.fa"

def find_alignment(query, dat, sec_len, proc_offset, rest, all_alignments, k_max):
    """
    Finds all substrings of dat that match query up to k_max mismatches.
    Input: query string, full dataset, desired sequence length, process number, maximum k mismatch.
    Output: all alignments as list, managed by process manager.
    """
    start_point = sec_len * proc_offset - len(query)
    end_point = sec_len * (proc_offset + 1) + rest
    if start_point < 0:
        start_point = 0
    data = dat[start_point:end_point] # define the data section per process
    m = len(query)
    n = len(data)

    lps = [0] * m # create a vector with length of query
    j = 0 # index for query[]

    lps_array(query, m, lps) # preprocess the query

    i = 0 # index for data[]
    k = 0
    while i < n:
        # find substring for k = 0
        if query[j] == data[i]:
            i += 1
            j += 1
        if j == m:
            all_alignments.append(start_point + i - j)
            j = lps[j - 1]
            k = 0
        # find substring with mismatches at end of query
        elif j == m - (k_max - k):
            all_alignments.append(start_point + i - j)
            print(start_point+i-j)
            j = lps[j - 1]
            i += k_max - k
            k = 0
        # find substring with k mismatches
        elif i < n and query[j] != data[i]:
            if k < k_max:
                match = False
                for z in range(1, k_max-k+1):
                    if j < m-z and i < n-z and query[j+z] == data[i+z]:
                        i += z
                        j += z
                        k += z
                        match = True
                        break
                if not match:
                    k = 0
                    if j != 0:
                        j = lps[j - 1]
                    else:
                        i += 1
            else:
                k = 0
                if j != 0:
                    j = lps[j - 1]
                else:
                    i += 1


def lps_array(query, m, lps):
    """
    Calculates LPS array of query for KMP string search.
    Input: query string, length of query and array of the same length.
    Output: LPS array
    """
    length = 0
    ind = 1

    while ind < m:
        if query[ind] == query[length]:
            length += 1
            lps[ind] = length
            ind += 1
        else:
            if length != 0:
                length = lps[length - 1]
            else:
                lps[ind] = 0
                ind += 1


if __name__ == '__main__':
    """
    Takes FASTA formatted genome and, given a query, finds substring that match with up to k mismatches.
    """
    samfile = pysam.FastaFile(sequence_filename)

    chromosome = list(range(1, 23))
    chromosome.extend(['X', 'Y'])

    query = input("Enter your query sequence: ")
    while not 3 < len(query) < 151:
        query = input("Please enter a query sequence between 4 and 150 bps: ")
    k_max = int(input("Enter the value of k: "))

    for c in chromosome: # look at chromosomes individually
        genome = samfile.fetch(str(c))
        num_processes = 4  # define number of processes
        section_length = int(len(genome) / num_processes) # divide the data into number of processes
        extra = len(genome) % num_processes

        manager = Manager()
        all_alignments = manager.list() # create a list of indices

        procs = []
        start = timeit.default_timer()
        leftovers = 0
        for i in range(num_processes):  # start each of the processes and add to procs list
            if i == num_processes - 1:
                leftovers = extra
            p = Process(target=find_alignment, args=(query, genome, section_length, i, leftovers, all_alignments, k_max))
            p.start()
            procs.append(p)
        print("Started all processes.")
        while any([p.is_alive() for p in procs]):  # check if any processes are still running
            pass
        stop = timeit.default_timer()
        print("All matches written to file. Time taken: {} seconds.".format(stop - start))

        # write output to file
        with open('results.txt', 'a') as f:
            if all_alignments != []:
                f.write("Found on chromosome " + str(c) + '\n')
                for i in sorted(all_alignments, key=int):
                    f.write(str(i) + ": " + str(genome[i:(i + len(query))]) + '\n')
        f.close()