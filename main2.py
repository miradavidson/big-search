import multiprocessing

import numpy as np

from Bio import SeqIO

record = SeqIO.read("sampleseq.fasta", "fasta")
gene = record.format("fasta")

def findAlignment(query, data):
    M = len(query)
    N = len(data)

    lps = [0] * M
    j = 0

    LPSArray(query, M, lps)

    i = 0
    while i < N:
        if query[j] == data[i]:
            i += 1
            j += 1

        if j == M:
            print("Found pattern at index " + str(i - j))
            j = lps[j - 1]
        elif i < N and query[j] != data[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1

def LPSArray(query, M, lps):
    len = 0
    i = 1

    while i < M:
        if query[i] == query[len]:
            len += 1
            lps[i] = len
            i += 1
        else:
            if len != 0:
                len = lps[len - 1]
            else:
                lps[i] = 0
                i += 1

listFormat = list(gene)
endOfDscp = list(gene).index('\n')
reformat = (listFormat[(endOfDscp+1):])
Bases = [s for s in reformat if s != '\n']
print(Bases)