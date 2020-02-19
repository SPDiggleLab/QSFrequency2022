#Python pipeline for generating protein files matching PAO1 or PAK, depending on the protein. Takes excel input files generated by NCBI-BLASTn and compares them to nucleotide FASTA files for the gene in question. Nucleotide sequences pulled out by that function are then translated. These files are manipulated using Pandas. Before translation, to ensure we don't have false truncations, we will delete all duplicates - aka any strains that are sequenced twice on different contigs, and may have been "truncated" onto two contigs. The protein sequence is then outputted into a txt file which will be used to generate dissimilarity graphs using R. All code in in UsefulFunctions (UF) python file.
import UsefulFunction_forlabs_kao_update as UF
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description = 'This program takes extracted lab excel files and translates them, then outputs the translation')
parser.add_argument("Excel", type=str, help="First enter excel file")
parser.add_argument("LabExcel", type=str, help="Enter lab excel file with gene IDs")
args = parser.parse_args()


nuc_seqs = UF.ReadExcel(args.Excel, args.LabExcel)

#Gets rid of duplicates, nucleotides before ATG, etc.
#list1CU = UF.Cleanup(list1[0]) #all accessions
nuc_seqs_clean = UF.Cleanup(nuc_seqs) #only accessions w/ no truncations

#translates DNA to AA, throws out anything after a stop codon
#list1AA = UF.Translate(list1CU)
prot_seqs = UF.Translate(nuc_seqs_clean)

name = str(args.LabExcel)[:-5]

df = pd.DataFrame(prot_seqs)

#export to CSV
df.to_csv(name + "-LasR.txt", sep='\t', index=False, header=False)