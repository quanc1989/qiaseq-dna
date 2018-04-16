import os
import sys
import collections
import itertools
import edlib
import cutadapt.adapters
# our modules
from umi_filter import reverseComplement


    
def create_primer_search_datastruct(primer_file):
    '''
    '''
    primer_kmer = collections.defaultdict(list)    
    primer_8bp = collections.defaultdict(list)
    primers = collections.defaultdict(list)
    primers_cutadapt = {}

    with open(primer_file,'r') as IN:
        for line in IN:
            contents = line.strip('\n').split('\t')
            primer = contents[-1]
            primer_8bp[primer[0:8]].append(primer)
            if primer in primers:
                raise Exception("Duplicate Primer Encountered : %s"%line)
            contents.append(len(primer))
            primers[primer] = contents

            # create k-mer index
            k=8
            kmers = set(''.join(itertools.islice(primer,i,i+k)) for i in range(len(primer)+1-k))
            for oligo in kmers:
                primer_kmer[oligo].append(primer)

            # create object for cutadapt
            revcomp = reverseComplement(primer)
            primers_cutadapt[primer] = cutadapt.adapters.Adapter(sequence=revcomp,where=cutadapt.adapters.BACK,max_error_rate=0.1,min_overlap=3)
            
            
    return (primer_8bp,primer_kmer,primers,primers_cutadapt)

def trim_primer(primer_datastruct,R1,R2):
    '''
    '''
    k = 8
    best_score = None
    best_primer = None
    primer_8bp, primer_kmer , primers, primers_cutadapt = primer_datastruct
    
    #candidates = primer_8bp[R1[0:8]]
    candidates = []
    if len(candidates) == 0: # could not look up first 8bp of the read in the index ; search the k-mer index
        candidates = set()
        for oligo in set(''.join(itertools.islice(R1[0:20],i,i+8)) for i in range(len(R1[0:20])+1-k)): # all 8-mers of the first 20 bp of the read
            if oligo in primer_kmer:
                for c in primer_kmer[oligo]:
                    candidates.add(c)
        if len(candidates) == 0: # exhaustive search over all primers
            candidates = primers        
       
    for p in candidates:

        p_len = primers[p][-1]

        if R1[0:p_len] == p: # check for exact match first
            return (p_len,p)
        
        else:
            alignment = edlib.align(R1[0:p_len],p,mode="NW")
            score =  float(alignment['editDistance'])/p_len
            if score <= 0.05:
                return (alignment['locations'][0][1],p)
            else:
                if best_score == None or score < best_score:
                    temp = alignment
                    best_primer = p
                    best_score = score

                    
    assert best_score != None, "Primer could not be scored correctly !"
    
    if best_score <= 0.10:
        return (temp['locations'][0][1],best_primer)
    else:
        return (-1,None)

def iterate_fastq(R1_fastq,R2_fastq):
    '''
    '''
    with open(R1_fastq,'r') as IN1 , open(R2_fastq,'r') as IN2:
        while True:
            R1_info =  (IN1.next(),IN1.next(),IN1.next(),IN1.next())
            R2_info =  (IN2.next(),IN2.next(),IN2.next(),IN2.next())
            yield (R1_info,R2_info)

    
def main(R1_fastq,R2_fastq,primer_file):
    '''
    '''
    # counters
    num_R1=0
    trimmed_R1=0
    trimmed_R2= 0

    # primer search datastruct
    primer_datastruct = create_primer_search_datastruct(primer_file)

    # output files
    R1_fastq_trimmed = R1_fastq + ".trimmed"
    R2_fastq_trimmed = R2_fastq + ".trimmed"

    batch = os.path.basename(R1_fastq).split(".")[2]    
    primer_trimming_info = os.path.join(os.path.dirname(R1_fastq_trimmed),"trimming_info.{}.txt".format(batch))

    primer_8bp,primer_kmer,primers,primers_cutadapt = primer_datastruct
    
    with open(R1_fastq_trimmed,'w') as OUT1, open(R2_fastq_trimmed,'w') as OUT2, open(primer_trimming_info,'w') as OUT3:
        for R1_info, R2_info in iterate_fastq(R1_fastq,R2_fastq):
            R1_id,R1_seq,R1_t,R1_qual = R1_info
            R2_id,R2_seq,R2_t,R2_qual = R2_info
            temp_id = R1_id.split(" ")[0]
            trim_pos,primer = trim_primer(primer_datastruct,R1_seq,R2_seq)
            if trim_pos == -1:
                OUT1.write(R1_id+R1_seq+R1_t+R1_qual)
                OUT2.write(R2_id+R2_seq+R2_t+R2_qual)
                OUT3.write(temp_id +"\tN/A\tN/A\n")                
            else:
                chrom,pos,strand,seq,p_len = primers[primer]
                primer_info = chrom+"-"+strand+"-"+pos+"\t"
                OUT3.write(temp_id+"\t"+primer_info)
                trimmed_R1+=1
                # trimed R1
                OUT1.write(R1_id+R1_seq[trim_pos+1:]+R1_t+R1_qual[trim_pos+1:])
                
                # trimmed R2; note: stripped new line here to get the match coordinates correctly
                read_seq = cutadapt.seqio.Sequence(name=R2_id,sequence=R2_seq.rstrip('\n'))
                match = primers_cutadapt[primer].match_to(read=read_seq)
                if match == None: # could not match primer
                    OUT2.write(R2_id+R2_seq+R2_t+R2_qual)
                    OUT3.write("N/A"+"\n")
                else:
                    trimmed_R2+=1
                    trimmed_seq = R2_seq[0:match.rstart]
                    trimmed_qual = R2_qual[0:match.rstart]                
                    OUT2.write(R2_id+trimmed_seq+"\n"+R2_t+trimmed_qual)
                    OUT3.write(primer_info+"\n")
                    
            num_R1+=1
            
    print "Total R1 Reads: {}".format(num_R1)
    print "R1 Reads Trimmed: {}".format(trimmed_R1)
    print "Total R2 Reads: {}".format(num_R1) 
    print "R2 Reads Trimmed: {}".format(trimmed_R2)
    

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3])    
