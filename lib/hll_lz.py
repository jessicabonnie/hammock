#!/usr/bin/env python
import Digest
import numpy as np
import sys
import argparse
from math import log2, pow

## Question: make sure this works for canonicalized kmers... is there a way to keep ordering information without losing this?
''' can we insure we select the same minimizer for both directions.... canonicalize substring.  Try same assembly same data but reverse complement. BED intervals, forward version and reverse version.  '''


def split_hash(stringh, bits=32,registerb=10):
    ''' This function splits a hash into a register and a remainder. It returns a tuple of the register and the leading zero count of the remainder.'''
    binh=format(stringh,f'0{bits}b')
    # if binh[0] == '1':
    #     print(f'I have a 1 and I am {len(binh)} long')
    register=int(binh[:registerb],2)
    remainder=binh[registerb:]
    # do we want numbers with no leading zeros to still receive a value of 1?
    lzc=remainder.find('1') #+ 1
    # print(stringh,binh,register,remainder,lzc)
    return (register,lzc)

def kmer_counts(line, klen, kmer_dict):
    ''' This function counts the number of times each kmer appears in a given string. It takes the string, kmer length, and a dictionary of kmers as optional parameters. It returns a dictionary of kmers and their counts.
    '''
    for i in range(len(line)-klen+1):
        kmer=line[i:i+klen]
        if kmer in kmer_dict:
            kmer_dict[kmer]+=1
        else:
            kmer_dict[kmer]=1
    return kmer_dict

def minimizer_lzcs(file_loc, window_size=40, klen=8, registerb=10, bits=32, actual=False, verbose=False):
    ''' This function calculates the LZCs (Leading Zero Count) for each register in a given file. It takes the file location, window size, kmer length, register size, and number of bits as optional parameters. It returns a list of LZCs for each register.
    '''
    counter=0
    kmer_dict= {}
    # names=[]
    lzcs = [0] * (2**registerb)
    howmany=0
    curmax=0
    registers=set()
    hashv=set()
    for line in open(file_loc,"r"):
        line = line.strip()
        if counter % 2 == 0:
            pass
        elif line=="":
            pass
            # names.append(line[1:])
        else:
            # if actual:
                # kmer_dict=kmer_counts(line,klen, kmer_dict)
            mnzrs= Digest.window_minimizer(line,w=window_size,k=klen, include_hash=True)
            howmany+=len(mnzrs)
            # print(len(mnzrs))
            # print(type(mnzrs))
            # x=len(line)
            # if counter % 501 == 0:
            #     print(counter)
            # print(line)
            if mnzrs == []:
                raise ValueError(f"There is a sequence of length {len(line)}, which cannot be analyzed using a window of {window_size} and a kmer length of {klen}.")
            ## analyze what the breakdown is for sequence lengths
            # can use the flanking sequence still though -- kmer that is part first kmer and last kmer, round up or down for klen/2 ... mimic fof stategy c when we change klen and winlen larger and smaller
            hashlist=[h for (m,h) in mnzrs]
  
            hashv.update(hashlist)
            curmax = max([max(hashlist),curmax])
            # kmers = [line[m:m+klen] for (m,h) in mnzrs]
            reg_lzc = [split_hash(stringh=h,bits=bits,registerb=registerb) for (m,h) in mnzrs]
            for (reg,lzc) in reg_lzc:
                registers.add(reg)
                if lzcs[reg] < lzc:
                    lzcs[reg] = lzc
        counter+=1
        # if counter > 20:
        #     exit()
    # print(kmer_dict)
    if verbose:
        print("Total Minimizer Values So Far: ",howmany)
        print(len(hashv))
        print("curmax: ", curmax)
        print("registers: ", registers)
        print("max register: ", max(registers))
        # print("Hash Set: ", hashv)
    return lzcs, len(hashv)

def hll_cardinality(lzcs, registerb, bits=32, verbose=False):
    #count number of registers with zeros
    # numb0=lzcs.count(0)
    if verbose: print('lzcs passed to hll_cardinality' ,lzcs)
    numb0 = 0
    z=0
    for lzc in lzcs:
        if lzc == 0:
            numb0 += 1
        z += pow(2,-lzc)

    #calculate the harmonic mean
    # z= sum([2**(-x) for x in lzcs])
    if verbose: print("z: ",z)
    #determine the bias correction
    m = pow(2,registerb)#2**registerb
    if (m <= 16): alpha = 0.673
    elif (m == 32): alpha = 0.697
    elif (m == 64): alpha = 0.709
    else: alpha = 0.7213/(1+1.079/m)
    if verbose: print("alpha: ", alpha)
    cardinality = alpha * m * m * (1/z); # Calculate the estimated cardinality
    if verbose: print("initial cardinality: ",cardinality)
    # Adjust for small cardinalities
    if cardinality <= 2.5 * m:
        if numb0 > 0:
            cardinality = m * log2(m / numb0)
    if verbose: print("cardinality after small cardinality adjustment: ",cardinality)
    # Adjust for large cardinalities using range correction
    
    if cardinality > (1/30) * pow(2,32): 
        cardinality = (-1*pow(2,32)) * log2(1 - (cardinality/pow(2,32)))
    if verbose: print("cardinality after large cardinality adjustment: ",cardinality)
    return cardinality

# def lzcs_union(lzcs1, lzcs2):
#     lu=[max(a,b) for a,b in zip(lzcs1,lzcs2)]
#     print(lu)
#     return lu

def lzcs_jaccard(lzcs1, lzcs2, registerb=8, bits=32):
    union = hll_cardinality([max(a,b) for a,b in zip(lzcs1,lzcs2)], registerb=registerb, bits=bits)
    intersection = hll_cardinality(lzcs1, registerb=registerb, bits=bits) + hll_cardinality(lzcs2, registerb=registerb, bits=bits) - union
    return intersection/union


def parse_clargs(argv):
    parser = argparse.ArgumentParser(description='Compute LZCs from a fasta file')
    parser.add_argument('fasta1', type=str, help='The path to the first fasta file')
    parser.add_argument('fasta2', type=str, help='The path to the second fasta file')
    # parser.add_argument('--table','-t', type=str, default=None, help='Table of kmer and window length options to try')
    parser.add_argument('--output', '-o', type=str, default=None, help='The output file path')
    parser.add_argument('--window_size','-w', type=str, default='40', help='A comma separated list of window sizes to use for minimizers')
    parser.add_argument('--klen', '-k', type=str, default='8', help='A comma separated list of kmer lengths to use for minimizers')
    parser.add_argument('--registerb', '-r', type=str, default='8', help='Comma separated list of values of the number of bits to use for the register')
    parser.add_argument('--bits', '-b', type=int, default=32, help='Number of bits of the number output by the hash function')
    args = parser.parse_args(args=argv)
    print(args)
    return args


def main(args):
    verbose=True
    args=parse_clargs(args)
    args.window_size = [int(x) for x in args.window_size.split(",")]
    args.klen = [int(x) for x in args.klen.split(",")]
    args.registerb = [int(x) for x in args.registerb.split(",")]
    if args.output is None:
        out=sys.stdout
    else:
        out=open(args.output, "w")
    out.write("Window\tKmer\tRBits\tJaccard\tCard1\tTruth1\tCard2\tTruth2\n")
    for i in args.window_size:
        for j in args.klen:
            for r in args.registerb:
                lzcs1, real_count1 = minimizer_lzcs(file_loc=args.fasta1, window_size=i, klen=j, registerb=r, actual=False, verbose=verbose)
                #99750935  
                lzcs2,real_count2 = minimizer_lzcs(file_loc=args.fasta2, window_size=i, klen=j, registerb=r, verbose=verbose)
                card1 = hll_cardinality(lzcs1,registerb=r)
                card2 = hll_cardinality(lzcs2,registerb=r)
                out.write(f"{i}\t{j}\t{r}\t{lzcs_jaccard(lzcs1, lzcs2, registerb=r, bits=args.bits)}\t{card1}\t{real_count1}\t{card2}\t{real_count2}\n")
    out.close()

    # lzcs1 = minimizer_lzcs(file_loc=args.fasta1, window_size=args.window_size, klen=args.klen)
    # #, window_size=args.window_size, klen=args.klen, registerb=args.registerb)

    # lzcs2 = minimizer_lzcs(file_loc=args.fasta2, window_size=args.window_size, klen=args.klen)
    #, window_size=args.window_size, klen=args.klen, registerb=args.registerb)

    # jaccard = lczs_jaccard(lzcs1,lzcs2)
    # print(f"Jaccard: {lczs_jaccard(lzcs1,lzcs2)}")
    # print(f"Hamming: {lczs_hamming(lzcs1,lzcs2)}")

    '''
    stress test:
    random sample of 10e9 numbers from 2^32-1. feed to split hash and then hll and see if cardinality is close to 10e9. Try with 4-8 register bits. Does it fill in all the registers? Then do some repeats with different sample sizes.... see if this affects the cardinality accuracy.
    
    
    
    '''

if __name__ == "__main__":
    # print(split_hash(4294967295,registerb=4))
    
    args = parse_clargs(sys.argv[1:])
    args.window_size = [int(x) for x in args.window_size.split(",")]
    args.klen = [int(x) for x in args.klen.split(",")]
    args.registerb = [int(x) for x in args.registerb.split(",")]
    if args.output is None:
        out=sys.stdout
    else:
        out=open(args.output, "w")
    out.write("Window\tKmer\tRBits\tJaccard\tCard1\tTruth1\tCard2\tTruth2\n")
    for i in args.window_size:
        for j in args.klen:
            for r in args.registerb:
                lzcs1, real_count1 = minimizer_lzcs(file_loc=args.fasta1, window_size=i, klen=j, registerb=r, actual=False)
                #99750935  
                lzcs2,real_count2 = minimizer_lzcs(file_loc=args.fasta2, window_size=i, klen=j, registerb=r)
                card1 = hll_cardinality(lzcs1,registerb=r)
                card2 = hll_cardinality(lzcs2,registerb=r)
                out.write(f"{i}\t{j}\t{r}\t{lzcs_jaccard(lzcs1, lzcs2, registerb=r, bits=args.bits)}\t{card1}\t{real_count1}\t{card2}\t{real_count2}\n")
    out.close()

    # lzcs1 = minimizer_lzcs(file_loc=args.fasta1, window_size=args.window_size, klen=args.klen)
    # #, window_size=args.window_size, klen=args.klen, registerb=args.registerb)

    # lzcs2 = minimizer_lzcs(file_loc=args.fasta2, window_size=args.window_size, klen=args.klen)
    #, window_size=args.window_size, klen=args.klen, registerb=args.registerb)

    # jaccard = lczs_jaccard(lzcs1,lzcs2)
    # print(f"Jaccard: {lczs_jaccard(lzcs1,lzcs2)}")
    # print(f"Hamming: {lczs_hamming(lzcs1,lzcs2)}")

    '''
    stress test:
    random sample of 10e9 numbers from 2^32-1. feed to split hash and then hll and see if cardinality is close to 10e9. Try with 4-8 register bits. Does it fill in all the registers? Then do some repeats with different sample sizes.... see if this affects the cardinality accuracy.
    
    
    
    '''