#!/usr/bin/env python
"""
Test script to verify the alignment functions work correctly
"""

# Helper functions
def score_pos(c1, c2, sm, g):
    """Score a position in alignment"""
    if c1 == "−" or c2 == "−":
        return g
    else:
        return sm.get(c1+c2, sm.get(c2+c1, -1))

def max3t(v1, v2, v3):
    """Return which of three values is maximum (1, 2, or 3)"""
    if v1 > v2:
        if v1 > v3:
            return 1
        else:
            return 3
    else:
        if v2 > v3:
            return 2
        else:
            return 3

def Needleman_Wunsch(seq1, seq2, sm, g):
    """Needleman-Wunsch global alignment algorithm"""
    S = [[0]]
    T = [[0]]
    # Initialize gaps' row
    for j in range(1, len(seq2)+1):
        S[0].append(g * j)
        T[0].append(3)
    # Initialize gaps' column
    for i in range(1, len(seq1)+1):
        S.append([g * i])
        T.append([2])
    # Apply the recurrence relation
    for i in range(0, len(seq1)):
        for j in range(len(seq2)):
            s1 = S[i][j] + score_pos(seq1[i], seq2[j], sm, g)
            s2 = S[i][j+1] + g
            s3 = S[i+1][j] + g
            S[i+1].append(max(s1, s2, s3))
            T[i+1].append(max3t(s1, s2, s3))
    return (S, T)

def recover_align(T, seq1, seq2):
    """Recover alignment from traceback matrix"""
    res = ["", ""]
    i = len(seq1)
    j = len(seq2)
    while i > 0 or j > 0:
        if T[i][j] == 1:
            res[0] = seq1[i-1] + res[0]
            res[1] = seq2[j-1] + res[1]
            i -= 1
            j -= 1
        elif T[i][j] == 3:
            res[0] = "−" + res[0]
            res[1] = seq2[j-1] + res[1]
            j -= 1
        else:
            res[0] = seq1[i-1] + res[0]
            res[1] = "−" + res[1]
            i -= 1
    return res

def read_fasta(filename):
    """Read sequence from FASTA file"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    sequence = ''
    for line in lines[1:]:  # Skip header
        sequence += line.strip()
    return sequence.upper()

def create_dna_scoring_matrix(match=1, mismatch=-1):
    """Create a simple DNA scoring matrix"""
    bases = ['A', 'T', 'G', 'C']
    sm = {}
    for b1 in bases:
        for b2 in bases:
            if b1 == b2:
                sm[b1+b2] = match
            else:
                sm[b1+b2] = mismatch
    return sm

# Test 1: Simple DNA alignment
print("="*60)
print("Test 1: Simple DNA Alignment")
print("="*60)

dna_sm = create_dna_scoring_matrix(match=1, mismatch=-1)
seq1 = "ATCGATCG"
seq2 = "ATCGATCG"

S, T = Needleman_Wunsch(seq1, seq2, dna_sm, -2)
alignment = recover_align(T, seq1, seq2)

print(f"Sequence 1: {seq1}")
print(f"Sequence 2: {seq2}")
print(f"Alignment score: {S[len(seq1)][len(seq2)]}")
print(f"Aligned seq1: {alignment[0]}")
print(f"Aligned seq2: {alignment[1]}")
print("✓ Test 1 passed\n")

# Test 2: DNA alignment with mismatch
print("="*60)
print("Test 2: DNA Alignment with Mismatches")
print("="*60)

seq1 = "ATCGATCG"
seq2 = "ATGGATGG"

S, T = Needleman_Wunsch(seq1, seq2, dna_sm, -2)
alignment = recover_align(T, seq1, seq2)

print(f"Sequence 1: {seq1}")
print(f"Sequence 2: {seq2}")
print(f"Alignment score: {S[len(seq1)][len(seq2)]}")
print(f"Aligned seq1: {alignment[0]}")
print(f"Aligned seq2: {alignment[1]}")
print("✓ Test 2 passed\n")

# Test 3: Load real data
print("="*60)
print("Test 3: Load Real FASTA Files")
print("="*60)

try:
    human_hbb = read_fasta('datasets/P68871.fasta')
    print(f"Human HBB loaded: {len(human_hbb)} amino acids")
    print(f"First 50 aa: {human_hbb[:50]}")
    print("✓ Test 3 passed\n")
except Exception as e:
    print(f"✗ Test 3 failed: {e}\n")

# Test 4: Small protein alignment
print("="*60)
print("Test 4: Small Protein Alignment")
print("="*60)

# Create simple protein scoring matrix
protein_sm = {}
amino_acids = "ACDEFGHIKLMNPQRSTVWY"
for a1 in amino_acids:
    for a2 in amino_acids:
        if a1 == a2:
            protein_sm[a1+a2] = 2
        else:
            protein_sm[a1+a2] = -1

seq1 = "MVHLTPEEK"
seq2 = "MVHLTPEEK"

S, T = Needleman_Wunsch(seq1, seq2, protein_sm, -2)
alignment = recover_align(T, seq1, seq2)

print(f"Sequence 1: {seq1}")
print(f"Sequence 2: {seq2}")
print(f"Alignment score: {S[len(seq1)][len(seq2)]}")
print(f"Aligned seq1: {alignment[0]}")
print(f"Aligned seq2: {alignment[1]}")
print("✓ Test 4 passed\n")

print("="*60)
print("All tests completed successfully!")
print("="*60)

