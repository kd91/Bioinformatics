# https://pythonforbiologists.com/dictionaries
# replace dashes with stars
# save output to a file. and check your output for correctness with ex-pasy output
# DOCSTRINGS in pycharm automatic :param filename: 'string' & :return: & :raises: explain what situation raises an exception

# add exceptions of file doesnt exist/ filename doesnt exist. Use try & except code
# dict for codon table

from sys import argv
import re

def read_file(filename):
  """
  :param filename:
  :return:
  """
  with open(filename,'r') as infile:
    read_string = infile.read()

  # print(read_seq)

  return read_string

def read_dna(seq):
  seq_dna = re.sub('[^acgtACGT]', '', seq).upper()

  return seq_dna

def complement(seq,comp):
  """
  :param seq:
  :param comp: cpmpliment
  :return:
  """
  rev_seq = seq[::-1]
  comp_rev_seq = ""
  for i in rev_seq:
    comp_rev_seq += comp[i]

  return comp_rev_seq

def find_reading_frame(seq,pos,codon_table):
  amino_acid = ""
  for i in range(pos, len(seq)-2, 3):
    if i+2 < len(seq):
      codon = seq[i:i+3]
      # print(codon)
      aa = codon_table[codon]
      amino_acid += aa
 
  return amino_acid

def write_frames(gene, filename):
  with open(filename, 'a') as outfile:
    outfile.write(gene)
    outfile.write("\n")


def test(seq_dna, out_file):
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

  # out_file = 'readingFrames.out.txt'
  input_pos = [0,1,2]
  for pos in input_pos:
    amino_acids = find_reading_frame(seq=seq_dna, pos=pos, codon_table=codon_table)
    print(amino_acids)
    comment = "5'3' Frame " + str(pos+1)
    write_frames(gene=comment, filename=out_file)
    write_frames(gene=amino_acids, filename=out_file)

  comp_dict = {'A': 'T', 'T':'A', 'C':'G', 'G':'C'}
  comp_rev_seq = complement(seq=seq_dna,comp=comp_dict)
  for pos in input_pos:
    rev_amino_acids = find_reading_frame(seq=comp_rev_seq, pos=pos, codon_table=codon_table)
    print(rev_amino_acids)
    comment = "3'5' Frame " + str(pos+1)
    write_frames(gene=comment, filename=out_file)
    write_frames(gene=rev_amino_acids,filename=out_file)


def main():
  my_filename = argv[1]
  out_file = argv[2]
  read_string = read_file(my_filename)
  seq_dna = read_dna(read_string)
  # seq_dna = 'CCTGCCCGTGCAGCA'
  print(seq_dna)
  test(seq_dna=seq_dna, out_file=out_file)
  

if __name__ == "__main__":
  main()