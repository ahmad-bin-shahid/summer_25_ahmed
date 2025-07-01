# importing libraries
import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import plotly.express as px
from Bio.Restriction import EcoRI
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
        
def valid(a):
    try:
        validation = all(base in 'ATGC' for base in a.upper())
        if validation:
            return True
        else:
            print("In-valid sequence")
    except FileNotFoundError:
        print("the file is not found")
def gc_content(a):
    try:
        gc_count=gc_fraction(a.upper())
        return gc_count
    except FileExistsError:
        print("invalid sequences")
        
def codon_freq(a):
    try:
        freq = len(a)//3
        return freq     
    except FileNotFoundError:
        print("file not open")    
def transcription(a):
    try:
        rna = a.transcribe()
        return rna[:50]
    except KeyError:
        print("invalid key")
        
def translation(a):
    try:
        protein = a.translate()
        return protein[:50]
    except KeyError:
        print("invalid key")
        
def Restriction(a):
    try:
        cut_site = EcoRI.search(a)
        return cut_site[:20]
    except KeyError:
        print("Invalid key")

def save_to_csv(results):
    try:
        if results:
            df = pd.DataFrame(results)
            df.to_csv("sequence_analysis.csv")
            print("Results saved to sequence_analysis.csv")
        else:
            print("No results to save.")
    except NameError:
        print("there is same name file already exist")


def plot_gc_content_bar(df):
    try:
        plt.figure(figsize=(12, 6))
        plt.bar(df['Sequence_ID'], df['GC_Content'], color='skyblue')
        plt.xlabel('Sequence ID')
        plt.ylabel('GC Content (%)')
        plt.title('GC Content per Sequence')
        plt.xticks(rotation=90, fontsize=8)
        plt.tight_layout()
        plt.savefig('gc_content_bar.png')
        plt.close()
    except FileNotFoundError:
        print("Data frame doesn't exist")

def plot_length_histogram(df):
    try:
        plt.figure(figsize=(10, 6))
        plt.hist(df['Length'], bins=20, color='lightcoral', edgecolor='black')
        plt.xlabel('Sequence Length')
        plt.ylabel('Frequency')
        plt.title('Distribution of Sequence Lengths')
        plt.tight_layout()
        plt.savefig('length_histogram.png')
        plt.close()
    except FileExistsError:
        print("Data frame doesn't exist")

def plot_GC_EcoRI_scatter(df):
    try:
        fig = px.scatter(df, x="GC_Content", y="EcoRI_Sites", hover_name = "Sequence_ID",title = "GC_Content VS EcoRI Sites")
        fig.write_html("restriction_scatter.html")
    except FileNotFoundError:
        print("Data Frame doesn't exist")
def analysis_and_visualization(a):
    try:
        results = []
        fasta_seq = SeqIO.parse(a,'fasta')
        for i in fasta_seq:
            i_id=i.id
            i_seq=i.seq
            if valid(i_seq) == True:
                length = len(i_seq)
                GC_Content = gc_content(i_seq)
                EcoRI_Sites = Restriction(i_seq)
                rna = transcription(i_seq)
                Protein = translation(i_seq)
                codon_usage_frequencies = codon_freq(i_seq)
                result = {
                    "Sequence_ID"   :   i_id,
                    "Length"        :   length,
                    "GC_Content"    :   GC_Content,
                    "EcoRI_Sites"   :   EcoRI_Sites,
                    "RNA (first 50 bases)"  : rna,
                    "Protein"       :   Protein,
                    "Codon usage frequencies"   : codon_usage_frequencies
                    }
            results.append(result)
        save_to_csv(results)
        df = pd.DataFrame(results)
        plot_gc_content_bar(df)
        plot_GC_EcoRI_scatter(df)
        plot_length_histogram(df)
    except FileNotFoundError:
        print("your file is not found")
        
if __name__ == '_main_':
    if len(sys.argv) != 2:
        sys.exit("Invalid arguments")
    seq = sys.argv[1]
    analysis_and_visualization()
