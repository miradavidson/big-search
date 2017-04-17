import mmap
import multiprocessing

from Bio import SeqIO
record = SeqIO.read("sampleseq.fasta", "fasta")
gene = record.format("fasta")

endOfDscp = gene.index('\n')
reformat = gene[(endOfDscp+1):]
Bases = reformat.replace("\n", "")
a = int(len(Bases)/4)
extra = len(Bases)%4

def findAlignment(query, dat, pno, rest):
    data = dat[(a*pno):(a*(pno+1)+rest)]
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
            print("Found pattern at index " + str((a*pno) + i - j) + ".")
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

query_here = input("Enter your query sequence: ")

if __name__ == '__main__':
    for i in range(3):
        p = multiprocessing.Process(target=findAlignment, args=(query_here, Bases, (1*i),0))
        p.start()

def multiprocess():
    Process(target=findAlignment, args = (query_here, Bases, 0,0)).start()
    Process(target=findAlignment, args = (query_here, Bases, 1,0)).start()
    Process(target=findAlignment, args = (query_here, Bases, 2,0)).start()
    Process(target=findAlignment, args = (query_here, Bases, 3,extra)).start()
    print("End of search.")