import pysam
import mmap
import timeit
from multiprocessing import Process, Manager

sequence_filename = "sampleseq.fasta"
sequence_name = "sp|P37840|SYUA_HUMAN"


def find_alignment(query, dat, sec_len, proc_offset, rest, all_alignments, k_max):
    start_point = sec_len * proc_offset - len(query)
    end_point = sec_len * (proc_offset + 1) + rest
    if start_point < 0:
        start_point = 0
    data = dat[start_point:end_point]
    m = len(query)
    n = len(data)

    lps = [0] * m
    j = 0


    lps_array(query, m, lps)

    i = 0
    k = 0
    while i < n:
        if query[j] == data[i]:
            i += 1
            j += 1

        if j == m:
            all_alignments.append(start_point + i - j)
            j = lps[j - 1]

        elif i < n and query[j] != data[i]:
            if k < k_max and query[j+1] == data[i+1]:
                i += 1
                j += 1
                k += 1
            else:
                if j != 0:
                    j = lps[j - 1]
                else:
                    i += 1

def lps_array(query, m, lps):
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
    samfile = pysam.FastaFile(sequence_filename)
    genome = samfile.fetch(sequence_name)

    num_processes = 4  # define number of processes
    section_length = int(len(genome) / num_processes)
    extra = len(genome) % num_processes

    manager = Manager()
    all_alignments = manager.list()

    query = input("Enter your query sequence: ")
    k_max = int(input("Enter the value of k: "))
    procs = []
    start = timeit.default_timer()
    for i in range(num_processes):  # start each of the processes and add to procs list
        p = Process(target=find_alignment, args=(query, genome, section_length, i, 0, all_alignments, k_max))
        p.start()
        procs.append(p)
    print("Started all processes.")
    while any([p.is_alive() for p in procs]):  # check if any processes are still running
        pass
    stop = timeit.default_timer()
    print("Time taken: {} seconds.".format(stop - start))
    print(all_alignments)