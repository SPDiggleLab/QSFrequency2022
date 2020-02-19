import pandas as pd

def ReadExcel(excel, labexcel):

	#read in the gene excel file generated from BLASTn
	#call on translation function to add protein column
	#get rid of duplicates to ensure no false duplicates

	df = pd.read_excel(excel, header=None)

	df.columns = ['Title', 'Sequence', 'Mismatches', 'Gaps', 'Percent Identity to Query', 'Query Coverage', 'Bitscore', 'E-Value']

	df['StrainID'] = df['Title'].str[45:-41]

	return df

	#now lab excel

	dflist = pd.read_excel(labexcel, index=False, header=None)

	dflist.columns = ['Number','StrainID', 'Blank', 'Country', 'IsolationAnimal', 'IsolationSource', 'Environment?', 'Investigator', 'Details', 'Details']

	list1 = list(dflist['StrainID'])

	df1 = pd.concat([df[df['StrainID'].str.contains(i+" ")] for i in list1])
		#df2 = df1[df1['StrainID'].isin(list1)]
		#df2 = df1[df1['StrainID'].isin(['ATCC BAA-2109 IP','AUS165'])]

	df1 = pd.DataFrame(df1)

	df2 = df1.iloc[:,1:]

	return df2
	    
def Cleanup(nuc_seqs):

	df = nuc_seqs.iloc[:,:2] # only keep 'title' and 'sequence' columns

	df['Sequencefix'] = df['Sequence'].str.replace('-','')

	searchfor = ['Y','W','M','R','S','K']
	df = df[~df['Sequencefix'].str.contains('|'.join(searchfor))]

	df['Sequencefix'].replace(to_replace='.*ATG', value='ATG', inplace=True)
	
	return df['Sequencefix']
	#df.to_csv("sequencefix.txt", sep='\t', index=False)

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


#df['Protein'] = (Translate(df4['Sequencefix']))



