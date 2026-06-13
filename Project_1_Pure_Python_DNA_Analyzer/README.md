# DNA Sequence Analyzer (Pure Python)

## Overview

This project is a beginner bioinformatics tool developed in pure Python.

The program reads a DNA sequence from a FASTA file and performs several common sequence analysis tasks, including nucleotide composition, GC content calculation, transcription, translation, protein analysis, ORF detection, and six reading frame generation.

The project was developed using the BRCA1 mRNA sequence (NM_007294) obtained from NCBI.

---

## Features

* Read DNA sequences from FASTA files
* Calculate sequence length
* Count nucleotide frequencies (A, T, C, G)
* Calculate GC content
* Generate reverse DNA sequence
* Generate reverse complement sequence
* Transcribe DNA into RNA
* Translate RNA into protein
* Detect Open Reading Frames (ORFs)
* Find the longest ORF
* Calculate dinucleotide frequencies
* Calculate protein length
* Calculate amino acid frequencies
* Calculate protein molecular weight
* Generate all six reading frames

---

## Project Structure

DNA_Sequence_Analyzer/

├── DNA_Sequence_Analyzer_Using_Python.py

├── NM_007294.fasta

└── README.md

---

## Requirements

* Python 3

No external libraries are required.

---

## Usage

Run the script:

```bash
python DNA_Sequence_Analyzer_Using_Python.py
```

Make sure the FASTA file is present in the same directory as the Python script.

---

## Dataset

BRCA1 mRNA sequence

Accession: NM_007294

Source: NCBI

---

## Concepts Learned

* FASTA file parsing
* DNA sequence manipulation
* Reverse complements
* Transcription
* Translation
* Genetic code and codons
* Open Reading Frames (ORFs)
* Reading frames
* Dinucleotide analysis
* Protein analysis

---

## Future Improvements

* Biopython implementation
* argparse command-line interface
* ORF histogram visualization
* GenBank support
* Codon usage analysis
* Automated report generation

---

## Author

Anjali Sinha

Bioinformatics Learning Project
