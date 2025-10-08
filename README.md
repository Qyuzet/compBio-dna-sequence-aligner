# DNA Sequence Aligner - Assignment 3

## Computational Biology (SCIE6062001)
**Odd Semester 2025–2026 AY**

This repository contains the implementation of Assignment 3: Pairwise and Multiple Sequence Alignment.

## Overview

This assignment implements three fundamental bioinformatics algorithms for sequence alignment:

1. **Global Alignment (Needleman-Wunsch)** - For aligning homologous genes
2. **Local Alignment (Smith-Waterman)** - For finding similar protein regions
3. **Multiple Sequence Alignment (MSA)** - For comparing multiple homologous proteins

## Project Structure

```
compBio-dna-sequence-aligner-dp/
├── sequence-aligner.ipynb    # Main Jupyter notebook with all solutions
├── task.md                   # Assignment requirements
├── test_alignment.py         # Test script to verify implementations
├── example/
│   ├── lab-5.ipynb          # Reference implementation from lab
│   └── blosum62.mat         # BLOSUM62 substitution matrix
└── datasets/
    ├── NC_000913.3-*.fasta  # E. coli genome sequences
    ├── NC_003198.-*.fasta   # Salmonella genome sequences
    ├── P00395.fasta         # Human COX1 protein
    ├── P00405.fasta         # Mouse COX1 protein
    ├── P68871.fasta         # Human hemoglobin beta
    ├── P02088.fasta         # Mouse hemoglobin beta
    ├── P01994.fasta         # Chicken hemoglobin beta
    ├── P02109.fasta         # Whale hemoglobin beta
    └── P02125.fasta         # Frog hemoglobin beta
```

## Problems Solved

### Problem 1: Global Alignment of Homologous Genes (Needleman–Wunsch)

**Objective:** Measure overall similarity between homologous recA genes from different bacterial species.

**Dataset:**
- *E. coli* recA gene (GenBank: NC_000913.3, coordinates 2,634,390–2,635,806)
- *Salmonella enterica* recA gene (GenBank: NC_003198.1, coordinates 2,711,800–2,713,200)

**Implementation:**
- Uses Needleman-Wunsch algorithm for global alignment
- Simple DNA scoring matrix (match=1, mismatch=-1, gap=-2)
- Reports alignment score, number of mismatches, and number of gaps

### Problem 2: Local Alignment of Similar Protein Regions (Smith–Waterman)

**Objective:** Identify regions of strong local similarity between related proteins.

**Dataset:**
- Human cytochrome c oxidase subunit I (COX1) - UniProt P00395
- Mouse cytochrome c oxidase subunit I (COX1) - UniProt P00405

**Implementation:**
- Uses Smith-Waterman algorithm for local alignment
- BLOSUM62 substitution matrix for protein scoring
- Gap penalty of -8
- Reports local alignment region, score, and gaps in aligned segment

### Problem 3: Multiple Sequence Alignment of Hemoglobin Proteins (MSA)

**Objective:** Study evolutionary conservation across homologous hemoglobin beta subunit proteins.

**Dataset:**
- Human HBB (UniProt P68871)
- Mouse HBB (UniProt P02088)
- Chicken HBB (UniProt P01994)
- Whale HBB (UniProt P02109)
- Frog HBB (UniProt P02125)

**Implementation:**
- Progressive alignment approach using pairwise Needleman-Wunsch
- BLOSUM62 substitution matrix
- Identifies conserved amino acid positions across all species
- Reports alignment scores and statistics

## Requirements

### Python Dependencies

```bash
pip install biopython jupyter
```

### Required Files

- Python 3.8 or higher
- Jupyter Notebook
- Biopython library
- Dataset files (included in `datasets/` folder)
- BLOSUM62 matrix (included in `example/` folder)

## Usage

### Running the Jupyter Notebook

1. Open Jupyter Notebook:
   ```bash
   jupyter notebook sequence-aligner.ipynb
   ```

2. Run all cells in order (Cell → Run All)

3. The notebook will:
   - Load all required sequences
   - Perform alignments for all three problems
   - Display results with statistics

### Running the Test Script

To verify the implementation works correctly:

```bash
python test_alignment.py
```

This will run basic tests on the alignment functions.

## Algorithms Implemented

### Needleman-Wunsch (Global Alignment)

- **Time Complexity:** O(mn) where m and n are sequence lengths
- **Space Complexity:** O(mn)
- **Use Case:** Aligning sequences of similar length end-to-end

### Smith-Waterman (Local Alignment)

- **Time Complexity:** O(mn)
- **Space Complexity:** O(mn)
- **Use Case:** Finding conserved regions within longer sequences

### Key Functions

- `Needleman_Wunsch(seq1, seq2, sm, g)` - Global alignment
- `Smith_Waterman(seq1, seq2, sm, g)` - Local alignment
- `recover_align(T, seq1, seq2)` - Traceback for global alignment
- `recover_align_local(S, T, seq1, seq2)` - Traceback for local alignment
- `read_submat_file(filename)` - Load substitution matrix
- `read_fasta(filename)` - Load sequences from FASTA files

## Results Summary

The notebook provides detailed results for each problem including:

- **Problem 1:** Alignment scores, mismatches, and gaps for bacterial recA genes
- **Problem 2:** Local alignment regions and identity percentages for COX1 proteins
- **Problem 3:** Pairwise alignment scores and conserved positions for hemoglobin proteins

## Notes

- The implementations are based on the algorithms taught in Lab 5
- For production use, consider established tools like MUSCLE, ClustalW, or MAFFT
- The MSA approach uses pairwise alignments; dedicated MSA algorithms provide better results
- All sequence data was obtained from NCBI GenBank and UniProt databases

## References

- Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. *Journal of Molecular Biology*, 48(3), 443-453.
- Smith, T. F., & Waterman, M. S. (1981). Identification of common molecular subsequences. *Journal of Molecular Biology*, 147(1), 195-197.
- Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from protein blocks. *Proceedings of the National Academy of Sciences*, 89(22), 10915-10919.

## Author

Assignment completed for Computational Biology (SCIE6062001)
Odd Semester 2025–2026 AY

## License

This is an academic assignment. Please refer to your institution's academic integrity policies.

