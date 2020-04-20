import pandas as pd
import matplotlib.pyplot as plt

def SplitExcel(excel):

	df = pd.read_excel(excel, header=None)

	df.columns = ['Title', 'Sequence', 'Mismatches', 'Gaps', 'Percent Identity to Query', 'Query Coverage', 'Bitscore', 'E-Value']
	
	df['Accessions'] = df['Title'].str[:4]

	return df
    
def Cleanup(excel):

	df = excel.iloc[:,:2] # only keep 'title' and 'sequence' columns

	df['Sequencefix'] = df['Sequence'].str.replace('-','')

	searchfor = ['Y','W','M','R','S','K']
	df = df[~df['Sequencefix'].str.contains('|'.join(searchfor))]

	df['Sequencefix'].replace(to_replace='.*ATG', value='ATG', inplace=True)
	
	return df['Sequencefix']

def Translate(sequence):
    # This function translates DNA sequence into amino acid sequence
	codontable = {
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
	'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
	}

	list2 = list(sequence)
	outList = []

	for row in list2:
		amino_acid = ''
		for i in range(0, len(row), 3):
			codon = row[i:i + 3]
			if len(codon) == 3:
				if codontable[codon] == '*':
					break
				else:
					amino_acid += codontable[codon]
		outList.append(amino_acid)

	return outList






