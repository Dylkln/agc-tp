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
    Take compressed fasta file and return sequences yield
    with a length greater or equal to minseqlen
    """

    with gzip.open(amplicon_file, "rb") as filin:
        for sequence in filin:
            if len(sequence) >= minseqlen:
                yield sequence.decode("UTF-8").strip()
            continue


def dereplication_fulllength(amplicon_file : str, minseqlen : int, mincount : int):
    """
    Take a compressed fasta file, min length of sequences, and their min count
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
        return chunk
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

    tmp = 0
    for nt in zip(alignment_list[0], alignment_list[1]):
        if nt[0] == nt[1]:
            tmp += 1

    identity = (tmp / len(alignment_list[0])) * 100

    return identity


def std(value_list):
    """
    Return the standard deviation of two values
    """

    return statistics.stdev(value_list)

def mean(value, liste):
    """
    Return the mean.
    """
    return (value/len(liste))


def detect_chimera(perc_identity_matrix):
    """
    Take an identity matrix
    return boolean : True if the target sequence is a chimera,
    False if the target sequence is not a chimera
    """

    standard_deviation = 0
    similarity_chunk1 = set()
    similarity_chunk2 = set()
    similarity = 0

    for sim in perc_identity_matrix:
        standard_deviation += std(sim)
        similarity_chunk1.add(sim[0])
        similarity_chunk2.add(sim[1])

    standard_deviation_mean = mean(standard_deviation, perc_identity_matrix)

    if len(similarity_chunk2) >= 2 or len(similarity_chunk1) >= 2:
        similarity = 1

    if standard_deviation_mean > 5 and similarity == 1:
        return True
    return False


def chimera_removal(amplicon_file : str, minseqlen : int,
    mincount : int, chunk_size : int, kmer_size : int):
    """
    Take a compressed fasta file, the minimum length for sequences, their min occurences,
    the chunk size and the kmer size
    return a generator of non chimera sequences
    """

    sequences = []
    occurences = []

    for seq_occ in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        sequences.append(seq_occ[0])
        occurences.append(seq_occ[1])

    segments = []
    kmer_dict = {}

    for index, seq in enumerate(sequences):
        segments.append(get_chunks(seq, chunk_size))
        kmer_dict = get_unique_kmer(kmer_dict, seq, index, kmer_size)

    mates = []

    for chunks in segments:
        for chunk in chunks:
            mates.append(search_mates(kmer_dict, chunk, kmer_size))

    parent_seq = common(mates[0], mates[1])

    chim_id = []

    chunk_list = [get_chunks(sequences[parent_seq[0]], chunk_size)]
    chunk_list += [get_chunks(sequences[parent_seq[1]], chunk_size)]

    for index, seq in enumerate(sequences):
        if seq not in parent_seq:
            chimera = get_chunks(seq, chunk_size)
            matrix = [[] for chunk in range(len(chimera))]

            for chunk in range(len(chunk_list)):
                for index2, chunk2 in enumerate(chimera):
                    matrix[index2].append(get_identity(nw.global_align(chunk2,
                        chunk_list[chunk][index2], gap_open = -1,
                        gap_extend = -1,
                        matrix = os.path.abspath(os.path.join(
                            os.path.dirname(__file__), "../agc")) + "/MATCH")))

            if detect_chimera(matrix):
                chim_id.append(index)

    for index, seq in enumerate(sequences):
        if index not in chim_id:
            yield [seq, occurences[index]]



def abundance_greedy_clustering(amplicon_file : str, minseqlen : int,
    mincount : int, chunk_size : int, kmer_size : int):
    """
    Take a compressed fasta file, the minimum length for sequences, their min occurences,
    the chunk size and the kmer size
    return the coccurence of each non chimera sequences
    """

    seq_count = chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)
    OTU = []

    for seq, count in seq_count:
        OTU.append(tuple(seq, count))

    return OTU


def write_OTU(OTU_list : list, output_file : str):
    """
    Take OTU list and write it in a fasta file
    """

    with open(output_file, "w") as filout:
        for index, OTU in enumerate(OTU_list):
            filout.write(f">OTU_{index + 1}, occurence : {OTU[1]}" + "\n")
            filout.write(fill(OTU[0]))
            filout.write("\n")


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

    OTU = abundance_greedy_clustering(args.amplicon_file, args.minseqlen,
        args.mincount, args.chunk_size, args.kmer_size)

    write_OTU(OTU, args.output_file)


if __name__ == '__main__':
    main()