#!/usr/bin/env python3

import argparse
import sys
from Bio import SeqIO

"""
Authors: Alexandre CORMIER
        Bioinformatics engineer
        SeBiMER, Ifremer

Creation Date: 2020-02-26
Modified on: 2020-02-27

Copyright (c): SeBiMER, february-2020

Modified from https://github.com/stajichlab/AAFTF/tree/master/AAFTF/vecscreen.py
Original authors : Jason Stajich and Jon Palmer
Licence: MIT
"""

def getArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i',dest="blastn",type=argparse.FileType('r'),required=True,help='Blast result against UniVecDB')
    parser.add_argument('-f',dest="fasta",type=str,required=True,help='Fasta to clean')
    parser.add_argument('-m',dest="minSize",type=int,default=500,help='Minimal contig size')
    parser.add_argument('-o',dest="output",type=str,required=True,help='Output name')

    arg = parser.parse_args()

    return arg

def main(args):

    ### Step - 1 - Blastn parsing results and attribution of hits
    VecHits, Excludes = parse_blastn(args.blastn)

    ### Step - 2 - Clean the fasta file and create a new record
    clean_fasta(VecHits, Excludes, args.fasta, args.output, args.minSize)

def parse_blastn(blastn):

    '''
    Blast header rows:
    qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue score qlen
    '''

    Excludes = {}
    VecHits = {}
    found_vector_seq = 0
    terminalDist = 200

    for row in blastn:
        qaccver,saccver,pid,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,score,qlen = row.split()
        scoverage = int(length) * 100 / int(qlen)
        #vecscreen https://www.ncbi.nlm.nih.gov/tools/vecscreen/about/#Moderate
        #says to use score here (I'm interpret as score not bitscore)
        #need to determine if match is terminal or if internal
        loc = [int(qstart), int(qend)]
        if loc[0] > loc[1]:
            loc = [loc[1],loc[0]]
        #check for location
        terminal = False
        position = None
        #considere hit as terminal if it's starts in the last 50bp
        #otherwise, considere it as internal
        if loc[0] <= terminalDist:
            terminal = True
            position = '5'
        if (int(qlen) - loc[1]) <= terminalDist:
            terminal = True
            position = '3'

        # Take account of the blast score
        # Default = weak
        Match = 0 # weak=0, moderate=1, strong=2
        score = int(score)
        if terminal:
            if score >= 19:
                Match = 1
            if score >= 24:
                Match = 2
        else:
            if score >= 25:
                Match = 1
            if score >= 30:
                Match = 2
        if Match == 0:
            continue

        if Match > 0:
            found_vector_seq += 1
            # remove the contig with more than half contamination
            if scoverage >= 50:
                Excludes[qaccver] = True
                # VecHits.pop(qaccver, None)
            # elif not qaccver in VecHits:
            elif not qaccver in VecHits:
                VecHits[qaccver] = [(saccver, int(qlen), loc, int(score), terminal, position)]
            else:
                VecHits[qaccver].append((saccver, int(qlen), loc, int(score), terminal, position))

    # delete key in VecHits that are in Excludes regardless of whether it is in the dictionary
    for qaccver in Excludes:
        VecHits.pop(qaccver, None)

    return(VecHits, Excludes)

def clean_fasta(VecHits, Excludes, fasta, output, minSize):

    trimTerminal = 0
    splitContig = 0
    f = open(output, "w")
    for record in SeqIO.parse(fasta, "fasta"):
        FiveEnd = 0
        ThreeEnd = len(record.seq)
        internals = []
        slicer = []
        sInt = []
        Seq = str(record.seq)
        if not record.id in VecHits and record.id not in Excludes:
            f.write('>'+record.description+'\n')
            while len(Seq) > 0:
                f.write(Seq[:70]+'\n')
                Seq = Seq[70:]
        elif record.id in Excludes:
            pass
        else:
            #VecHits contains list of tuples of information, if terminal, then just truncate
            #off the closest side. Also, need to check if multiple intervals are within 50
            #bp of each other, that whole interval is removed.
            #should be able to accomplish above with the several rounds that it runs with,
            #so split on internal and trim terminal. done.
            for hit in VecHits[record.id]:
                ID,length,loc,score,terminal,pos = hit
                if terminal and pos == '5':
                    if loc[1] > FiveEnd:
                        FiveEnd = loc[1]
                elif terminal and pos == '3':
                    if loc[0] < ThreeEnd:
                        ThreeEnd = loc[0]
                else: #internal hits to add to list
                    if not loc in internals:
                        internals.append(loc)

            #now sort intervals
            sInt = sorted(internals, key=lambda x: int(x[0]))
            #now construct slicing list
            if len(sInt) < 1:
                slicer = [FiveEnd,ThreeEnd]
            else:
                slicer = [FiveEnd]
                for x in sInt:
                    slicer = slicer + x
                slicer.append(ThreeEnd)
            paired_slicer = list(group(slicer, 2))

            if len(paired_slicer) < 2:
                print('Terminal trimming {:} to {:}'.format(record.id, paired_slicer))
                newSeq = Seq[paired_slicer[0][0]:paired_slicer[0][1]]
                if len(newSeq) >= minSize:
                    f.write('>'+record.description+'\n')
                    while len(newSeq) > 0:
                        f.write(newSeq[:70]+'\n')
                        newSeq = newSeq[70:]

            else:
                print('Spliting contig {:} into {:}'.format(record.id, paired_slicer))
                for num,y in enumerate(paired_slicer):
                    newSeq = Seq[y[0]:y[1]]
                    if len(newSeq) >= minSize:
                        f.write('>'+record.id+'.'+str(num+1)+'\n')
                        while len(newSeq) > 0:
                            f.write(newSeq[:70]+'\n')
                            newSeq = newSeq[70:]

    f.close()

def group(lst, n):
    for i in range(0, len(lst), n):
        val = lst[i:i+n]
        if len(val) == n:
            yield tuple(val)

if __name__ == '__main__':
    args = getArgs()
    main(args)
