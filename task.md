# Computational Biology (SCIE6062001)
**Odd Semester 2025–2026 AY**

## Assignment 3: Pairwise and Multiple Sequence Alignment

---

### 🧬 Problem 1 – Global Alignment of Homologous Genes (Needleman–Wunsch)

**Objective:**  
Apply a global alignment algorithm (Needleman–Wunsch) to measure overall similarity between homologous genes from different bacterial species.

**Dataset:**
- *E. coli* recA gene – GenBank accession **NC_000913.3**, coordinates **2,634,390–2,635,806**
- *Salmonella enterica* recA gene – GenBank accession **NC_003198.1**, coordinates **2,711,800–2,713,200**

**Task:**
1. Retrieve both sequences in FASTA format from NCBI.  
2. Perform a global alignment using the Needleman–Wunsch algorithm.  
3. Report:
   - Alignment score  
   - Number of mismatches  
   - Number of gaps  

---

### 🔍 Problem 2 – Local Alignment of Similar Protein Regions (Smith–Waterman)

**Objective:**  
Apply a local alignment algorithm (Smith–Waterman) to identify regions of strong local similarity between two related proteins.

**Dataset:**
- Human cytochrome c oxidase subunit I (COX1) – UniProt **P00395**
- Mouse cytochrome c oxidase subunit I (COX1) – UniProt **P00405**

**Task:**
1. Retrieve both protein sequences in FASTA format.  
2. Perform a local alignment to detect the most similar subsequences.  
3. Report:
   - Local alignment region and score  
   - Gaps in the aligned segment  

---

### 🧩 Problem 3 – Multiple Sequence Alignment of Hemoglobin Proteins (MSA)

**Objective:**  
Perform a multiple sequence alignment (MSA) on homologous hemoglobin beta subunit proteins to study evolutionary conservation.

**Dataset:**  
Collect protein sequences from UniProt:
- Human HBB – UniProt **P68871**  
- Mouse HBB – UniProt **P02088**  
- Chicken HBB – UniProt **P01994**  
- Whale HBB – UniProt **P02109**  
- Frog HBB – UniProt **P02125**

**Task:**
1. Download sequences in FASTA format.  
2. Align all sequences using an MSA tool (e.g., **Biopython**, **ClustalW**, **MUSCLE**).  
3. Report:
   - Alignments  
   - Locations (indices) of identical amino acids  

---

### 📄 Notes
- Make sure all FASTA files are properly labeled and included in your submission.  
- Clearly indicate which tool or Python library was used for each alignment.  
- Include screenshots or visualizations of alignments if available.  
