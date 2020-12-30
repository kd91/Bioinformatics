# Reading FASTA sequence into amino acid sequence reading frames

Description of the project:

Write a Python3 script named readingFrames.py that will translate a FASTA sequence of nucleotides into the amino acid sequences of all 6 possible reading frames.Your program should be able to read the sequence from a file whose filename is input as a command line argument from the command prompt:

    python readingFrames.py seq.fasta.txt

An example of the translation into the 6 reading frames is the Expasy tool (http://us.expasy.org/tools/dna.html (Links to an external site.)), which also finds open reading frames (ORFs) and color-codes them with a reddish highlight. 

Your code should output in the amino acid sequence to a plain text file named readingFrames.out.txt in the format output by Expasy  (except there would be no highlighting, no color and no bold font in the output). 

Ex:
5'3' Frame 1
P A R A A R G S S T S L G A S P S S A L G S T D L F S R R K V H W N W K K W I Y L P S K F K K Y K Met S F Met L C R K S Stop .....

An ORF is the part of a reading frame that has the potential to code for a protein or peptide. It is a continuous stretch of codons beginning with a start codon and ending with a stop codon. You are NOT required to find the ORFs, but only to provide the 6 reading frames in which an ORF could be picked out.

You may assume that the gene is in proper FASTA format (ie- it contains only letters encoding nucleotides, numbers and whitespace). You may assume the letters are restricted to capital and lowercase A,C,T and G.

Your script should be able to translate the gene found here in the file - seq.fasta.txt 

 
