"""
Description: Reads the dna sequence from an input file and returns six reading frames in a new output file
Functions within are: reading the input file, reading string in dna strand format, finding the complement strand of dna, finding the reading frame, writing the reading frame into an output file and writing the amino acids reading frame to the output file.
Author: Khushbu Desai
"""

from sys import argv
import re
from os import path

def read_file(filename):
  """
  If opens and reads the input file contents.
  :param : filename: input file name
  :return: entire string taht is read from the input file.
  """
  with open(filename,'r') as infile:
    read_string = infile.read()

  return read_string

def read_dna(seq):
  """
  Converts a string input and removes characters invalid for a dna sequence
  param: seq: string input with any characters
  return: returns string with no other characters except acgtACGT and converts them to uppercase
  """
  seq_dna = re.sub('[^acgtACGT]', '', seq).upper()

  return seq_dna

def complement(seq,comp):
  """
  Find the complement of characters in a string based on key values from a input dictionary
  :param: seq: string of characters
  :param : comp: dictionary to find compliments of any character that must be substituted with a given key value
  :return: string with characters reveresed than seq string and new values from the comp dictionary keys
  """
  rev_seq = seq[::-1]
  comp_rev_seq = ""
  for i in rev_seq:
    comp_rev_seq += comp[i]

  return comp_rev_seq

def find_reading_frame(seq,pos,codon_table):
  """
  To find a set of three characters or three nucletides (codon) of a string
  param: seq: input string of nucleotides/ dna sequence
  param: pos: the index number from which the reading of codons should commence
  param: codon_table: dictionary with keys of three nucleotides (codon) with a value corresponding amino acid character
  return: returns a string of amino acids that are corresponding to the codon in the reading frame
  """
  amino_acid = ""
  for i in range(pos, len(seq)-2, 3):
    if i+2 < len(seq):
      codon = seq[i:i+3]
      aa = codon_table[codon]
      amino_acid += aa
 
  return amino_acid

def write_frames(gene, filename):
  """
  Appends values into an output file. 
  param: gene: values that need to be appended into the output file
  param: filename: output file that is appended
  """
  with open(filename, 'a') as outfile:
    outfile.write(gene)


def write_amino_acids(seq_dna, out_file):
  """
  It writes the sequence of amino acids into different reading frames to the output file
  param: seq_dna: input string of amino acids
  param: out_file: the output file that needs to be appended with reading frames
  """
  codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'} 
  # gene dictionary reference: https://pythonforbiologists.com/dictionaries

  input_pos = [0,1,2]
  for pos in input_pos:
    amino_acids = find_reading_frame(seq=seq_dna, pos=pos, codon_table=codon_table)
    # print(amino_acids)
    comment = "5'3' Frame " + str(pos+1) + "\n"
    write_frames(gene=comment, filename=out_file)
    write_frames(gene=amino_acids +  "\n\n", filename=out_file)

  comp_dict = {'A': 'T', 'T':'A', 'C':'G', 'G':'C'}
  comp_rev_seq = complement(seq=seq_dna,comp=comp_dict)
  for pos in input_pos:
    rev_amino_acids = find_reading_frame(seq=comp_rev_seq, pos=pos, codon_table=codon_table)
    # print(rev_amino_acids)
    comment = "3'5' Frame " + str(pos+1) + "\n"
    write_frames(gene=comment, filename=out_file)
    write_frames(gene=rev_amino_acids + "\n\n",filename=out_file)


def main():
  my_filename = argv[1]
  if path.exists(my_filename):
    filename, _ = path.splitext(my_filename)
    out_file = filename + ".out.txt"
    read_string = read_file(my_filename)
    seq_dna = read_dna(read_string)
    write_amino_acids(seq_dna=seq_dna, out_file=out_file)   
  else:
    print("Error: Input file does not exist")
  

if __name__ == "__main__":
  main()