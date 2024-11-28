from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt 

def read_fasta_file(file):
    stringio = StringIO(file.read().decode("utf-8"))
    fasta_file = ""
    for seq_record in SeqIO.parse(stringio, "fasta"): #gene.fna is the file name and fasta is type of file
        fasta_file = seq_record.seq
    return fasta_file

def count_nucleotides(filename):
    #create dictionary to store nucleotides
    gene = filename
    nucleotides = {"A": 0, 'C': 0, 'G': 0, 'T':0 }
    #"ACTGTGTCGATCGATCAGCTGGATCAAT"
    #count the nucleotides, then add to the dictionary
    nucleotides["A"] = gene.count('A') #count every string "A" in the gene and store it in the dictionary
    nucleotides["C"] = gene.count('C')
    nucleotides["G"] = gene.count('G')
    nucleotides["T"] = gene.count('T')

    # return dictionary nucleotides
    return nucleotides


def count_codons(filename):
    gene = filename
    gene = gene.transcribe()
    codons = {} #Empty dictionary
    # "ACTGTGTCGATCGATCAGCTGGATCAAT"

    # Loop through the data gene
    for i in range(0, len(gene), 1):
        # determine seq of 3 nucleotides = codon
        codon = gene[i:i+3]
        #Check if codon already in codons dictionary, if codon in codons dictionary, we add 1 to the count
        if codon in codons:
            codons[codon] += 1
        
        # if codon not in codons dictionary, we add it and start the count from 1
        else:
            codons[codon] = 1

    #sorted codon total
    codons = dict(sorted(codons.items(), key=lambda item: item[1] )) 

    # Delete seq that is not a codon (seq that is not 3 nucleotides)
    codons = dict((k, v) for k, v in codons.items() if len(k) == 3)

    cleaned_data = {str(seq): value for seq, value in codons.items()}
    codons = cleaned_data
    #codons as result
    return codons

AMINO_ACIDS = [
            'Alanine', 'Arginine', 'Asparagine', 'Aspartic acid', 'Cysteine',
            'Glutamine', 'Glutamic acid', 'Glycine', 'Histidine', 'Isoleucine',
            'Leucine', 'Lysine', 'Methionine', 'Phenylalanine', 'Proline',
            'Serine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valine'
        ]

AMINO_ACIDS_CODE = {# Phenylalanine
    'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine',
    # Leucine
    'UUA': 'Leucine', 'UUG': 'Leucine', 'CUU': 'Leucine', 'CUC': 'Leucine', 'CUA': 'Leucine', 'CUG': 'Leucine',
    # Isoleucine
    'AUU': 'Isoleucine', 'AUC': 'Isoleucine', 'AUA': 'Isoleucine',
    # Methionine (Start Codon)
    'AUG': 'Methionine(Start)',
    # Valine
    'GUU': 'Valine', 'GUC': 'Valine', 'GUA': 'Valine', 'GUG': 'Valine',
    # Serine
    'UCU': 'Serine', 'UCC': 'Serine', 'UCA': 'Serine', 'UCG': 'Serine', 'AGU': 'Serine', 'AGC': 'Serine',
    # Proline
    'CCU': 'Proline', 'CCC': 'Proline', 'CCA': 'Proline', 'CCG': 'Proline',
    # Threonine
    'ACU': 'Threonine', 'ACC': 'Threonine', 'ACA': 'Threonine', 'ACG': 'Threonine',
    # Alanine
    'GCU': 'Alanine', 'GCC': 'Alanine', 'GCA': 'Alanine', 'GCG': 'Alanine',
    # Tyrosine
    'UAU': 'Tyrosine', 'UAC': 'Tyrosine',
    # Stop Codons
    'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',
    # Histidine
    'CAU': 'Histidine', 'CAC': 'Histidine',
    # Glutamine
    'CAA': 'Glutamine', 'CAG': 'Glutamine',
    # Asparagine
    'AAU': 'Asparagine', 'AAC': 'Asparagine',
    # Lysine
    'AAA': 'Lysine', 'AAG': 'Lysine',
    # Aspartic Acid
    'GAU': 'Aspartic Acid', 'GAC': 'Aspartic Acid',
    # Glutamic Acid
    'GAA': 'Glutamic Acid', 'GAG': 'Glutamic Acid',
    # Cysteine
    'UGU': 'Cysteine', 'UGC': 'Cysteine',
    # Tryptophan
    'UGG': 'Tryptophan',
    # Arginine
    'CGU': 'Arginine', 'CGC': 'Arginine', 'CGA': 'Arginine', 'CGG': 'Arginine', 'AGA': 'Arginine', 'AGG': 'Arginine',
    # Glycine
    'GGU': 'Glycine', 'GGC': 'Glycine', 'GGA': 'Glycine', 'GGG': 'Glycine'}

def bar_chart_viz(dataframe, title, color):
    fig, ax = plt.subplots(figsize=(20, 10)) #initialize plot size

    colors = plt.cm.viridis(np.linspace(0, 1, len(dataframe['Codon']))) #set color of bars
    codons_plot = plt.bar(dataframe["Codon"], dataframe["Count"], color=colors) #make barplot

    for bar in codons_plot:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2.0, yval + 25, "--- "+str(yval), ha="center", fontsize=10, rotation=90)

    plt.margins(x=0.01) #set distance betweene bars and axis
    plt.title(title)
    plt.xticks(rotation=90)

    ax.set_ylim(0, max(dataframe['Count']) * 1.1)

    plt.show() #show the plot





    def dump_code():
        fig, ax = plt.subplots(figsize=(20, 10))
        ax1 = ax.bar(df["Nucleotide"], df['Count'], color=['#2589BD', '#FE4A49', '#B9E3C6', '#F8C537'])
        for bar in ax1:
            yval = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2.0, yval + (yval*0.02), str(yval), ha='center', fontsize=12)
        # Add grid horizontal line
        ax.yaxis.grid(True, linestyle='-.', linewidth=0.5, color='black')

        plt.margins(x=0.01)
        ax.set_xticks(range(len(df["Nucleotide"])))
        ax.set_xticklabels(['Adenine', 'Cytosine', 'Guanine', 'Thymine'])
        plt.title('Number of Nucleotides in Gene', fontsize=16)
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        # Remove axis lines
        #ax.xaxis.set_visible(False)
        #ax.yaxis.set_visible(False)

        # Adjust y-axis limits to provide more space for labels
        
        ax.set_ylim(0, max(df['Count']) * 1.1)
        st.pyplot(fig)