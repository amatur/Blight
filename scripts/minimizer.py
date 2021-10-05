# seq="ATGCGATATCGTAGGCGTCGATGGAGAGCTAGATCGATCGATCTAAATCCCGATCGATTCCGAGCGCGATCAAAGCGCGATAGGCTAGCTAAAGCTAGCA"


import sys
filepath=sys.argv[1]
Kmer=int(sys.argv[2])
M=int(sys.argv[3])

# file1=open("/Users/Sherlock/amaturWS/mini")
file1=open(filepath,"r")
Lines = file1.readlines()

count = 0
# Strips the newline character
for seq in Lines:
    seq=seq.strip()
    rev=seq[::-1]

    rev=rev.replace("A","X")
    rev=rev.replace("T","A")
    rev=rev.replace("X","T")
    rev=rev.replace("C","X")
    rev=rev.replace("G","C")
    rev=rev.replace("X","G")

    count += 1
    #print("Seq: ", seq, ": \n")

    L=len(seq)
    for i in range(0, L-Kmer+1):
        sub_f=seq[i:i+Kmer]
        sub_r=rev[L-Kmer-i:L-i]

        min="ZZZZZZZZZZZZZ"
        min_fonly="ZZZZZZZZZZZZZ"
        for j in range(0, Kmer-M+1):
            sub2=sub_f[j:j+M]
            if sub2 < min_fonly:
                min_fonly=sub2

            if sub2 < min:
                min=sub2
            sub2=sub_r[j:j+M]
            
            if sub2 < min:
                min=sub2

                
        print(sub_f,min,min_fonly)

file1.close()
