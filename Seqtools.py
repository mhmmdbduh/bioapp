import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import numpy as np
from utils import read_fasta_file, count_nucleotides, count_codons, AMINO_ACIDS_CODE


st.title("DNA Sequence Tools and Analysis")


filename = st.file_uploader("Upload your sequence file", type=["txt", "fasta", "fa", "fastq", "fq","fna", "faa"])

st.subheader("Sequence Tools")
st.write("Essential tool for sequence manipulationn. Transcription, Translation, finding Reverse compliment of your sequence.")
transcribe, translate, read_fasta = st.columns(3)

if filename is not None:
    fasta_file = read_fasta_file(filename)

if transcribe.button("Transcribe"):    
    if filename is not None:
        
        st.write("Transcribe your sequence")        
        st.write(fasta_file.transcribe())
    else:
        st.write("Please upload a fasta file")

if translate.button("Translate"):
    if filename is not None:
        
        st.write("Translate your sequence")        
        st.write(fasta_file.translate())
    else:
        st.write("Please upload a fasta file")
    
if read_fasta.button("Read File"):
    if filename is not None:
        
        st.write(fasta_file)     
    else:
        st.write("Please upload a fasta file") 


st.subheader("Sequence Analysis")

count_nucleotide, count_codon, count_aminoacids = st.columns(3)

if count_nucleotide.button("Count Nucleotide"):
    if filename is not None:
        #fasta_file = read_fasta_file(filename)
        
        count = count_nucleotides(fasta_file)
        st.write(f"Your sequence has {count['A']} Adenine (A), {count['C']} Cytosine (C), {count['G']} Guanine (G), {count['T']} Thymine (T)")
        
        df = pd.DataFrame(list(count.items()),columns=['Nucleotide', 'Count'])
        st.dataframe(df)
        fig = px.bar(df, x='Nucleotide', y='Count', color='Nucleotide', color_continuous_scale='viridis')
        fig.update_layout(
            title=dict(text='Count of Nucleotides in Gene', font=dict(size=30), automargin=True, yref='paper', xref='paper', x=0.4),
)
        fig.update_layout(
            xaxis=dict(
            tickmode='array',
            tickvals=[i for i in range(len(df['Nucleotide']))],
            ticktext=df['Nucleotide'],
            tickangle=0,
            tickfont=dict(size=10)  # Change font size here
            ),
            # bargap=0.2,
            height=720  # Custom height
        )
        fig.update_traces(texttemplate='%{y}', textposition='outside', textangle=0, textfont_size=10)
        st.plotly_chart(fig, use_container_width=False)
    else:
        st.write("Please upload a fasta file")



if count_codon.button("Count Codons"):
    if filename is not None:
        count = count_codons(fasta_file)
        df = pd.DataFrame(list(count.items()), columns=['Codon', 'Count'])
        count_codon.write("Count Codons") 
        count_codon.dataframe(df)
        
        fig = px.bar(df, x='Codon', y='Count', color='Count', color_continuous_scale='viridis')
        fig.update_layout(
            title=dict(text='Count of Codons in Gene', font=dict(size=30), automargin=True, yref='paper', xref='paper', x=0.4),
)
        fig.update_layout(
            xaxis=dict(
            tickmode='array',
            tickvals=[i for i in range(len(df['Codon']))],
            ticktext=df['Codon'],
            tickangle=90,
            tickfont=dict(size=10)  # Change font size here
            ),
            # bargap=0.2,
            height=720  # Custom height
        )
        #fig.update_traces(texttemplate='--- %{y}', textposition='outside', textangle=-90, textfont_size=100)
        st.plotly_chart(fig, use_container_width=False)

        colors = plt.cm.viridis(np.linspace(0, 1, len(df['Codon']))) #set color of bars
        fig, ax = plt.subplots(figsize=(20, 10))
        ax1 = ax.bar(df["Codon"], df['Count'], color = colors)
        for bar in ax1:
            yval = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2.0, yval + (yval * 0.1), "--- "+str(yval), ha="center", fontsize=10, rotation=90)
        # Add grid horizontal line
        ax.yaxis.grid(True, linestyle='-.', linewidth=0.5, color='black')

        plt.margins(x=0.01)

        plt.title('Number of Nucleotides in Gene', fontsize=16)
        plt.xticks(rotation=90)
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
    else:
        st.write("Please upload a fasta file")




if count_aminoacids.button("Count Amino Acids"):
    if filename is not None:
        count = count_codons(fasta_file)
        
        amino_group = list(count.keys())
        df = pd.DataFrame(list(count.items()), columns=['codon', 'Count'])
        amino_group = [AMINO_ACIDS_CODE[amino] for amino in amino_group]
        df["Amino Acid"] = amino_group
        df = df.groupby("Amino Acid").sum().reset_index()
        df = df.drop(columns=["codon"])
        df = df.sort_values(by='Count', ascending=True)
        count_aminoacids.write("Count Amino Acids")
         
        count_aminoacids.dataframe(df)
        fig = px.bar(df, x='Amino Acid', y='Count', color='Amino Acid')
        fig.update_layout(
            title=dict(text='Count of Codons in Gene', font=dict(size=30), automargin=True, yref='paper', xref='paper', x=0.4),
)
        fig.update_layout(
            xaxis=dict(
            tickmode='array',
            tickvals=[i for i in range(len(df['Amino Acid']))],
            ticktext=df['Amino Acid'],
            tickangle=90,
            tickfont=dict(size=10)  # Change font size here
            ),
            # bargap=0.2,
            height=720,  # Custom height
            margin=dict(l=0, r=0, t=100, b=0)
        )
        fig.update_traces(texttemplate='--- %{y}', textposition='outside', textangle=-90, textfont_size=1000)
        st.plotly_chart(fig,theme="streamlit", use_container_width=False)
    else:
        st.write("Please upload a fasta file")

if st.button("visualize Sequences"):
    if filename is not None:
        
        # Create a grid of boxes representing the nucleotides
        total_nucleotides = len(fasta_file)
        grid_size = (int(np.ceil(np.sqrt(total_nucleotides))), int(np.ceil(np.sqrt(total_nucleotides))))
        # grid_size = (5, 5)  # Define the size of the grid
        fig, ax = plt.subplots(figsize=(20, 20))

        # Create a color map for nucleotides
        color_map = {'A': 'red', 'C': 'blue', 'G': 'green', 'T': 'orange'}

        # Plot each nucleotide in the grid
        for i in range(grid_size[0]):
            for j in range(grid_size[1]):
                if i * grid_size[1] + j < len(fasta_file):
                    nucleotide = fasta_file[i * grid_size[1] + j]
                    color = color_map.get(nucleotide, 'white')
                    rect = plt.Rectangle((j, grid_size[0] - i - 1), 1, 1, facecolor=color)
                    ax.add_patch(rect)
                    # ax.text(j + 0.5, grid_size[0] - i - 1 + 0.5, nucleotide, ha='center', va='center', fontsize=12)

        #Set the limits and labels
        ax.set_xlim(0, grid_size[1])
        ax.set_ylim(0, grid_size[0])
        ax.set_xticks(np.arange(grid_size[1]) + 0.5)
        ax.set_yticks(np.arange(grid_size[0]) + 0.5)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        # Remove the outline
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        # remove background
        ax.set_facecolor('black')
        # make outside plot also black
        fig.patch.set_facecolor('black')
        # Remove axis lines
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        # Create legend
        legend_handles = [plt.Rectangle((0, 0), 1, 1, color=color_map[nuc]) for nuc in color_map]
        legend_labels = list(color_map.keys())
        ax.legend(legend_handles, legend_labels,  loc='upper center', bbox_to_anchor=(0.1, 1.0), ncol=4)
        plt.title('Visualization of Nucleotides in Chromosome 17', fontdict={'color': 'white', 'fontsize': 20})
        st.pyplot(fig)




        
    