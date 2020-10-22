#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "KLEIN Dylan"
__copyright__ = "Universite de Paris"
__credits__ = ["KLEIN Dylan"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "KLEIN Dylan"
__email__ = "klein.dylan@outlook.com"
__status__ = "Developpement"


def isfile(path):
    """
    Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """
    Retrieves the arguments of the program.
    Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file : str, minseqlen : int):
    """
    Take a file and return sequences yield 
    with a length greater or equal to minseqlen
    """
    
    with gzip.open(amplicon_file, "rb") as filin:
        for sequence in filin:
            if len(sequence) >= minseqlen:
                yield sequence.strip()
            continue


def dereplication_fulllength(amplicon_file : str, minseqlen : int, mincount : int):
    """
    Take fasta file, min length of sequences, and their min count
    return a yield of the sequence and its count
    """
    
    seq_occ_dict = {}
    
    for sequence in read_fasta(amplicon_file, minseqlen):
        if sequence not in seq_occ_dict.keys():    
            seq_occ_dict[sequence] = 0
        seq_occ_dict[sequence] += 1

    seq_occ_dict = dict(sorted(seq_occ_dict.items(), key = lambda x : x[1], reverse = True))

    for sequence, count in seq_occ_dict.items():
        if count >= mincount:
            yield [sequence, count]


def get_chunks(sequence : str, chunk_size : int):
    """
    Take a sequence and a chunk size
    return a chunk yield
    """
    
    chunk = []
    for i in range(0, len(sequence), chunk_size):
        if chunk:    
            if len(sequence[i:i + chunk_size]) != len(chunk[0]):
                break
        chunk.append(sequence[i:i + chunk_size])

    if len(chunk) >= 4:
        yield chunk
    else:    
        return None


def cut_kmer(sequence : str, kmer_size : int):
    """
    Take a sequence and a kmer_size
    return a kmer yield
    """
    
    for i in range(len(sequence) - kmer_size + 1):
        yield sequence[i : i + kmer_size]


def get_unique(liste : list):
    """
    take a list
    return the list with unique items
    """
    
    return set(liste)


def common(liste1 : list, liste2 : list):
    """
    Take two lists
    return the common items between them
    """
    
    return list(set(liste1).intersection(liste2))


def get_unique_kmer(kmer_dict : dict, sequence : str, id_seq : int, kmer_size : int):
    """
    Take a kmer dictionnary, a sequence, a seq ID, a kmer_size
    return unique kmer
    """
    
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer not in kmer_dict:
            kmer_dict[kmer] = []
        kmer_dict[kmer].append(id_seq)

    return kmer_dict



def search_mates(kmer_dict : dict, sequence : str, kmer_size : int):
    """
    Take a kmer dictionnary, a sequence and a kmer size
    return 8 sequences with the most affinity for the input sequence
    """

    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size) 
        if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]



def get_identity(alignment_list : list):
    """
    Take an alignment list of two aligned sequences
    return the identity percentage of them
    """
    
    tmp = len(common(alignment_list[0], alignment_list[1]))
    identity = (tmp / len(alignement_list[0])) * 100
    
    return identity


def detect_chimera(perc_identity_matrix : matrix):
    """
    Take an identity matrix
    """

    

    pass


def chimera_removal(amplicon_file : str, minseqlen : int, mincount : int, chunk_size : int, kmer_size : int):
    """
    
    """
    pass


def abundance_greedy_clustering(amplicon_file : str, minseqlen : int, mincount : int, chunk_size : int, kmer_size : int):
    """
    
    """
    pass


def write_OTU(OTU_list : list, output_file : str):
    """
    
    """
    pass


def fill(text : str, width = 80):
    """
    Split text with a line return to respect fasta format
    """

    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))




#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    
    args = get_arguments()

    dereplication_fulllength(args.amplicon_file, args.minseqlen, args.mincount)

    kmer_dict = get_unique_kmer({}, "TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAG", 0, 8)
    print(kmer_dict)
    kmer_dict = get_unique_kmer(kmer_dict, "GGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGC", 1, 8)
    print(kmer_dict)
    kmer_dict = get_unique_kmer(kmer_dict, "GGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCC", 2, 8)
    print(kmer_dict)

if __name__ == '__main__':
    main()