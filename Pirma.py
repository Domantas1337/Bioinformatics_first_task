import re

START_CODON = 0
STOP_CODON = 1


def transform_to_a_protein_sequence(sequence):
	protein_table = {
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

	protein = ""
	for i in range(0, len(sequence), 3):
		codon = sequence[i:i + 3]
		protein += protein_table[codon]

	return protein

def remove_short_fragments(codon_pairs):
	codon_pair_index = 0
	while True:
		if(codon_pair_index == len(codon_pairs)):
			break
		if(codon_pairs[codon_pair_index][1] - codon_pairs[codon_pair_index][0] < 97):
			del codon_pairs[codon_pair_index]
			continue
		if(codon_pair_index == len(codon_pairs) - 1):
			break
		codon_pair_index = codon_pair_index + 1



def find_longest_codon_pairs(start_stop_pairs):
	
	current_stop_codon = -1
	longest_codon_pairs = []


	for codon_pair in start_stop_pairs:
		if (current_stop_codon != codon_pair[1]):
			longest_codon_pairs.append(codon_pair)
			current_stop_codon = codon_pair[1]
		else:
			continue

	return longest_codon_pairs

def find_start_stop_pairts(start_positions, stop_positions):
	
	start_stop_pairs = []
	last_codon_pair = (-1, -1)

	current_stop_position_index = 0

	while start_positions and stop_positions:

		if(start_positions[0] > stop_positions[current_stop_position_index]):
			del stop_positions[current_stop_position_index]
		else:
			if(stop_positions[current_stop_position_index] - start_positions[0] >= 4
				and  (stop_positions[current_stop_position_index] - start_positions[0]) % 3 == 0):
				start_stop_pairs.append((start_positions[0], stop_positions[current_stop_position_index]))
			del start_positions[0]


	return start_stop_pairs



def find_stop_codon(sequence, position):
	stop_codons = ['TAA', 'TAG', 'TGA']

	for codon in stop_codons:
		result = sequence.find(codon, position, position + 3)
		
		if(result != -1):
			return result
	
	return -1 

def find_codons(sequence):
	
	start_codon = 'ATG'

	start_positions = [m.start() for m in re.finditer(start_codon, sequence)]
	print(start_positions)
	stop_positions = []

	position = 0
	result = -1
	
	while position != len(sequence) - 1:
		
		result = find_stop_codon(sequence, position)
		
		if(result != -1):
			stop_positions.append(result)

		position = position + 1

	return start_positions, stop_positions

def find_reverse_compliment(sequence):
	reverse_sequence = reversed(sequence)

	reverse_compiment = ""

	for char in reverse_sequence:
		if char == 'A':
			reverse_compiment += 'T'
		elif char == 'T':
			reverse_compiment += 'A'
		elif char == 'G':
			reverse_compiment += 'C'
		else:
			reverse_compiment += 'G'
	return reverse_compiment


def parse_fasta(fileName):

	sequence = []

	with open(fileName, 'r') as file:

		for line in file:
			if line.startswith(">"):
				sequence_name = line.strip().lstrip(">")
				sequence = ""
			else:
				sequence += line.strip()
	return sequence



def find_codon_rate(protein_sequence):
	codons = {'A':0, 'R':0, 'N':0, 'D':0, 'C':0	
			 ,'Q':0, 'E':0, 'G':0, 'H':0, 'I':0
			 ,'L':0, 'K':0, 'M':0, 'F':0, 'P':0	
	         ,'O':0, 'S':0, 'U':0, 'T':0, 'W':0	
			 ,'Y':0, 'V':0, 'B':0, 'Z':0, 'X':0	
			 ,'J':0}
	
	dicodons = {}
	for codon in codons:
		for second_codon in codons:
			dicodons[codon + second_codon] = 0
	
	for j in range(len(protein_sequence)):
		
		codons[protein_sequence[j][0]] = 1
		
		for i in range(0, len(protein_sequence[j]) -2):
			codons[protein_sequence[j][i + 1]] += 1
			dicodons[protein_sequence[j][i] + protein_sequence[j][i+1]] += 1

	for i in codons:
		codons[i] = codons[i] / (len(protein_sequence[j]) - 1)

def process_a_sequence(sequence):
	start_positions, stop_positions = find_codons(sequence)
	pairs = find_start_stop_pairts(start_positions, stop_positions)
	longest_pairs = find_longest_codon_pairs(pairs)

	remove_short_fragments(longest_pairs)
	print(longest_pairs)


	protein_sequence = []
	for i in range(0, len(longest_pairs)):
		protein_sequence.append(transform_to_a_protein_sequence(sequence[longest_pairs[i][0]:longest_pairs[i][1] + 3]))

	print(protein_sequence)
	find_codon_rate(protein_sequence)

	return (longest_pairs, protein_sequence)

# lines = fasta_string.strip().split("\n")
# description = lines[0]

# sequence = "".join(lines[1:])
	return sequence


fileName = "mamalian2.fasta"


sequence = parse_fasta(fileName)
reverse_compiment = find_reverse_compliment(sequence)


result = process_a_sequence(reverse_compiment)

# Output results
# print(f"Start codon positions: {start_positions}")
# print(f"Stop codon positions: {stop_positions}")

