
NUCLEOTIDES = ['-', 'A', 'T', 'C', 'G']

CODON_TO_RESIDUE = {'---': '-',
                    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                    'TTA': 'L', 'TCA': 'S', 'TAA': '.', 'TGA': '.',
                    'TTG': 'L', 'TCG': 'S', 'TAG': '.', 'TGG': 'W',
                    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
                    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G',
                    '***': '*',
                    }

CODONS = list(CODON_TO_RESIDUE)

RESIDUES = ['-',
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
            '.',
            ]
