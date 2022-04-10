import  pandas as pd
import numpy as np
from Bio import SeqIO
import subprocess

def process_data_histone(df, res_folder):
    df['chr'] = df['chr'].str.replace('chr', '')
    df['start'] = df['pos'].astype(int) - 999
    df['start'] =df['start'].clip(lower=0)
    df['end'] = df['pos'].astype(int) + 1000
    df['bed'] = df['chr'].astype(str) + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
    df['bed'].to_csv('{}/bed_histone'.format(res_folder), sep='\t', index=False, header=False)
    df[['chr', 'pos', 'strand', 'ref', 'alt']].to_csv('{}/snps_histone'.format(res_folder), sep='\t', index=False,
                                                   header=False)

    exit_code = subprocess.Popen(
        "samtools faidx genome.fa -r {}/bed_histone -o {}/seq_histone".format(
            res_folder, res_folder), shell=True, stdout=subprocess.PIPE).stdout.read()
    print(exit_code)
    code_ =['A','C','G','T','R','Y','S','W','K','M']
    code = {
                    #              A, C, G, T
                    'A': np.array([1, 0, 0, 0]),
                    'C': np.array([0, 1, 0, 0]),
                    'G': np.array([0, 0, 1, 0]),
                    'T': np.array([0, 0, 0, 1]),
                    'R': np.array([0.5 , 0    , 0.5  , 0   ]),
                    'Y': np.array([0   , 0.5  , 0    , 0.5 ]),
                    'S': np.array([0   , 0.5  , 0.5  , 0   ]),
                    'W': np.array([0.5 , 0    , 0    , 0.5 ]),
                    'K': np.array([0   , 0    , 0.5  , 0.5 ]),
                    'M': np.array([0.5 , 0.5  , 0    , 0   ]),
    }

    fasta_sequences = SeqIO.parse(open('{}/seq_histone'.format(res_folder)), 'fasta')
    bed = pd.read_csv('{}/snps_histone'.format(res_folder), header=None, sep='\t')
    bed.columns=['chr', 'pos','strand','ref','alt']
    refs = bed['ref'].tolist()
    alts = bed['alt'].tolist()
    chrs = bed['chr'].tolist()
    strands = bed['strand'].tolist()
    poss=bed['pos'].tolist()
    ref_matrix = []
    alt_matrix = []
    chr_tmp =[]
    strand_tmp =[]
    pos_tmp =[]
    ref_tmp =[]
    alt_tmp =[]
    for i, (fasta, chr, strand, pos, ref, alt) in enumerate(zip(fasta_sequences,chrs, strands, poss, refs,alts)):
        name, sequence = fasta.id, str(fasta.seq).upper()
        coding_ref = []
        coding_alt = []

        tmp = list(sequence)
        if len(tmp)== 2000:

            if ref.upper() not in code_:
                continue

            if alt.upper() not in code_:
                continue


            for idx, j in enumerate(tmp):

                if idx==999:
                    coding_alt.append(code[alt.upper()])
                    coding_ref.append(code[ref.upper()])
                else:
                    coding_ref.append(code[j])
                    coding_alt.append(code[j])

            ref_matrix.append(coding_ref)
            alt_matrix.append(coding_alt)
            chr_tmp.append(chr)
            pos_tmp.append(pos)
            strand_tmp.append(strand)
            ref_tmp.append(ref)
            alt_tmp.append(alt)

    bed_histone = pd.DataFrame({'chr':chr_tmp, 'pos':pos_tmp, 'strand':strand_tmp, 'ref':ref_tmp, 'alt':alt_tmp })

    return bed_histone, np.array(ref_matrix),np.array(alt_matrix)


def process_data_methy(df, res_folder):
    df['chr'] = df['chr'].str.replace('chr', '')
    df['start'] = df['pos'].astype(int) - 499
    df['start'] =df['start'].clip(lower=0)
    df['end'] = df['pos'].astype(int) + 500
    df['bed'] = df['chr'].astype(str) + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
    df['bed'].to_csv('{}/bed_methy'.format(res_folder), sep='\t', index=False, header=False)
    df[['chr', 'pos', 'strand', 'ref', 'alt']].to_csv('{}/snps_methy'.format(res_folder), sep='\t', index=False,
                                                      header=False)
    exit_code = subprocess.Popen(
        "samtools faidx genome.fa -r {}/bed_methy -o {}/seq_methy".format(
            res_folder, res_folder), shell=True, stdout=subprocess.PIPE).stdout.read()

    print(exit_code)

    code_ =['A','C','G','T','R','Y','S','W','K','M']
    code = {
                    #              A, C, G, T
                    'A': np.array([1, 0, 0, 0]),
                    'C': np.array([0, 1, 0, 0]),
                    'G': np.array([0, 0, 1, 0]),
                    'T': np.array([0, 0, 0, 1]),
                    'R': np.array([0.5 , 0    , 0.5  , 0   ]),
                    'Y': np.array([0   , 0.5  , 0    , 0.5 ]),
                    'S': np.array([0   , 0.5  , 0.5  , 0   ]),
                    'W': np.array([0.5 , 0    , 0    , 0.5 ]),
                    'K': np.array([0   , 0    , 0.5  , 0.5 ]),
                    'M': np.array([0.5 , 0.5  , 0    , 0   ]),
    }

    fasta_sequences = SeqIO.parse(open('{}/seq_methy'.format(res_folder)), 'fasta')
    bed = pd.read_csv('{}/snps_methy'.format(res_folder), header=None, sep='\t')
    bed.columns=['chr', 'pos','strand','ref','alt']
    refs = bed['ref'].tolist()
    alts = bed['alt'].tolist()

    chrs = bed['chr'].tolist()
    strands = bed['strand'].tolist()
    poss=bed['pos'].tolist()
    chr_tmp =[]
    strand_tmp =[]
    pos_tmp =[]
    ref_tmp =[]
    alt_tmp =[]

    ref_matrix = []
    alt_matrix = []
    for i, (fasta, chr, strand, pos, ref, alt) in enumerate(zip(fasta_sequences,chrs, strands, poss, refs,alts)):

        if strand=='-':
            # tmp = list(sequence.reverse_complement())
            name, sequence = fasta.id, str(fasta.reverse_complement().seq).upper()
        else:

            name, sequence = fasta.id, str(fasta.seq).upper()
        coding_ref = []
        coding_alt = []

        tmp = list(sequence)
        if len(tmp) == 1000:

            if ref.upper() not in code_:
                continue

            if alt.upper() not in code_:
                continue

            for idx, j in enumerate(tmp):

                if idx==499:
                    coding_alt.append(code[alt])
                    coding_ref.append(code[ref])
                else:
                    coding_ref.append(code[j])
                    coding_alt.append(code[j])

            ref_matrix.append(coding_ref)
            alt_matrix.append(coding_alt)
            chr_tmp.append(chr)
            pos_tmp.append(pos)
            strand_tmp.append(strand)
            ref_tmp.append(ref)
            alt_tmp.append(alt)

    bed_methy = pd.DataFrame({'chr':chr_tmp, 'pos':pos_tmp, 'strand':strand_tmp, 'ref':ref_tmp, 'alt':alt_tmp })

    return bed_methy,np.array(ref_matrix),np.array(alt_matrix)
