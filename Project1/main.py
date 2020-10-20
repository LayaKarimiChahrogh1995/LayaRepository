#This file is a main file for project of Protein/Nucleotide Sequence Alignment.
#It contains interface to the four parts which make a whole, those are the following:
#1. Global and local alignment of nucleotide sequences using given matrix as a score matrix
#2. Global and local alignments of protein sequences using PAM and BLOSUM score matrices
#3. Alignment of protein and nucleotide sequences, depending on the type
#4. Alignment of multiple sequences using CLUSTALW algorithm

import first, second, third, utilities
from screens import screen

def main():
        screen()

if __name__ == "__main__":
        main()
