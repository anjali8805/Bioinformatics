# project header 

print("=" * 50)
print("DNA SEQUENCE ANALYZER")
print("=" * 50)

# function to read a FASTA file and return the DNA sequence as a string
def read_fasta(filename):
    sequence = ""

 # Open the FASTA file in read mode, Read the file line by line, Skip the header line (starts with >) and Convert sequence to uppercase for consistency

    with open(filename, "r") as file:
        for line in file:
            if not line.startswith(">"):
                sequence += line.strip()

    return sequence.upper()

#  Read DNA sequence from FASTA file
dna = read_fasta("NM_007294.fasta")

# to display the length of the DNA sequence.
print("\n===== BASIC SEQUENCE ANALYSIS =====")
print("Length:", len(dna))

# Find the count of each nucleotide (A, T, C, G) in the DNA sequence
a_count = dna.count("A")
t_count = dna.count("T")
c_count = dna.count("C")
g_count = dna.count("G")

# to display the count of each nucleotide in the DNA sequence.
print("A:", a_count)
print("T:", t_count)
print("C:", c_count)
print("G:", g_count)

# to calculate the GC content of the DNA sequence and display it as a percentage.
gc_content = (g_count + c_count) / len(dna) * 100
print("GC Content:", round(gc_content, 2), "%")

# to find the reverse DNA
reverse_dna = dna[::-1]
print("Reverse DNA:", reverse_dna)

# Generate the reverse complement of the DNA sequence
complement = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

reverse_complement = ""

for base in dna[::-1]:
    reverse_complement += complement[base]

print("Reverse Complement:", reverse_complement)

#ORF finding
print("\n===== ORF ANALYSIS =====")
stop_codons = ["TAA","TGA","TAG"]

longest_orf = ""
longest_length = 0
longest_start = 0
longest_stop = 0

for i in range (len(dna)-2):
    codon = dna[i:i+3]
    if codon == "ATG":
        print("Position of start codon is found", i)
        # Loop through the sequence after the start codon to find the stop codon
        for j in range(i+3, len(dna)-2, 3):
            current_codon = dna[j:j+3]
            if current_codon in stop_codons:
                print("Position of stop codon is found", j)

                # Extract the ORF
                orf = dna[i:j+3]
                orf_length = len(orf)

                print("Start Position:", i)
                print("Stop Position:", j)  
                print("ORF Length:", orf_length)
                print("ORF:", orf)
                
                # find the longest ORF
               
                if orf_length > longest_length:
                    longest_length = orf_length
                    longest_orf = orf
                    longest_start = i
                    longest_stop = j
                break

               
print("\n===== LONGEST ORF =====")
print("Start Position:", longest_start)
print("Stop Position:", longest_stop)
print("Length:", longest_length)
print("ORF:", longest_orf)

print("\n===== TRANSCRIPTION =====")
#DNA to RNA transcription
rna = dna.replace("T", "U")
print("RNA Sequence:", rna)

## Standard genetic code: DNA codons mapped to amino acids
print("\n===== TRANSLATION =====")
genetic_code = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",

    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",

    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",

    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

protein = ""
for i in range(0, len(rna), 3):
    codon = rna[i:i+3]
    if len(codon) == 3:
        protein += genetic_code[codon]

# Print complete protein sequence
print("Protein Sequence:", protein)

# Print protein fragments separated by stop codons
print("Protein Fragments:")

for fragment in protein.split("*"):
    if fragment:
        print(fragment)


# Dinucleotide Frequencies
dinucleotide_counts = {}
for i in range(len(dna) -1):
    dinucleotide = dna[i: i+2]
    if dinucleotide in dinucleotide_counts:
        dinucleotide_counts[dinucleotide] += 1
    else:
        dinucleotide_counts[dinucleotide] = 1
print("Dinucleotide Frequencies:", dinucleotide_counts)

# Protein Properties - Protein Length
print("\n===== PROTEIN ANALYSIS =====")
protein_length = len(protein)
print ("Protein Length:", protein_length)

# Amino Acid Frequencies
amino_acid_counts = {}
for amino_acid in protein:
    if amino_acid in amino_acid_counts:
        amino_acid_counts[amino_acid] += 1
    else:
        amino_acid_counts[amino_acid] = 1
print("Amino Acid Frequencies:", amino_acid_counts)

# Protein Molecular Weight

weights = {
    'A': 89.1,
    'R': 174.2,
    'N': 132.1,
    'D': 133.1,
    'C': 121.2,
    'Q': 146.2,
    'E': 147.1,
    'G': 75.1,
    'H': 155.2,
    'I': 131.2,
    'L': 131.2,
    'K': 146.2,
    'M': 149.2,
    'F': 165.2,
    'P': 115.1,
    'S': 105.1,
    'T': 119.1,
    'W': 204.2,
    'Y': 181.2,
    'V': 117.1
}

molecular_weight = 0
for amino_acid in protein:
    if amino_acid in weights:
        molecular_weight += weights[amino_acid]
print("Protein Molecular Weight:", molecular_weight, "Da")

# ORFs in All 6 Reading Frames
# Forward reading frames
frame1 = dna
frame2 = dna[1:]
frame3 = dna[2:]

# Reverse reading frames
frame4 = reverse_complement
frame5 = reverse_complement[1:]
frame6 = reverse_complement[2:]

print("\n===== ORFs in All 6 Reading Frames =====")
print(
    "Frame +1:", frame1,
    "\nFrame +2:", frame2,
    "\nFrame +3:", frame3,
    "\nFrame -1:", frame4,
    "\nFrame -2:", frame5,
    "\nFrame -3:", frame6
)