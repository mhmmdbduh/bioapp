import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from utils import read_fasta_file, count_nucleotides, count_codons, AMINO_ACIDS_CODE


st.title("DNA Sequence Tools and Analysis")


filename = st.file_uploader("Upload your sequence file", type=["txt", "fasta", "fa", "fastq", "fq","fna", "faa"])

st.subheader("Sequence Tools")
st.write("Essential tool for sequence manipulationn. Transcription, Translation, finding Reverse compliment of your sequence.")
transcribe, translate, read_fasta = st.columns(3)



if transcribe.button("Transcribe"):    
    if filename is not None:
        fasta_file = read_fasta_file(filename)
        st.write("Transcribe your sequence")        
        st.write(fasta_file.transcribe())
    else:
        st.write("Please upload a fasta file")

if translate.button("Translate"):
    if filename is not None:
        fasta_file = read_fasta_file(filename)
        st.write("Translate your sequence")        
        st.write(fasta_file.translate())
    else:
        st.write("Please upload a fasta file")
    
if read_fasta.button("Read File"):
    if filename is not None:
        fasta_file = read_fasta_file(filename)
        st.write(fasta_file)     
    else:
        st.write("Please upload a fasta file") 


st.subheader("Sequence Analysis")

count_nucleotide, count_codon, count_aminoacids = st.columns(3)

if count_nucleotide.button("Count Nucleotide"):
    if filename is not None:
        #fasta_file = read_fasta_file(filename)
        
        count = count_nucleotides(filename)
        st.write(f"Your sequence has {count['A']} Adenine (A), {count['C']} Cytosine (C), {count['G']} Guanine (G), {count['T']} Thymine (T)")
        df = pd.DataFrame([count])
        st.dataframe(df)

    else:
        st.write("Please upload a fasta file") 
if count_codon.button("Count Codons"):
    if filename is not None:
        count = count_codons(filename)
        df = pd.DataFrame(list(count.items()), columns=['Codon', 'Count'])
        count_codon.write("Count Codons") 
        count_codon.dataframe(df)
if count_aminoacids.button("Count Amino Acids"):
    if filename is not None:
        count = count_codons(filename)
        
        amino_group = list(count.keys())
        df = pd.DataFrame(list(count.items()), columns=['codon', 'Count'])
        amino_group = [AMINO_ACIDS_CODE[amino] for amino in amino_group]
        df["Amino Acid"] = amino_group
        df = df.groupby("Amino Acid").sum()
        df = df.drop(columns=["codon"])
        count_aminoacids.write("Count Amino Acids")
         
        count_aminoacids.dataframe(df)

st.subheader("Sequence Visualization")
st.write("Visualize your sequence data using different visualization techniques")


        
    