# Project Summaries
# 1. Basic FASTA Handling and Conversion
Scripts: FASTA_to_GB_biopy.py, fasta_file_handling.py, fasta_script.py
Goal: Automate reading, writing, and converting FASTA files into GenBank format using Biopython.
Skills Used:
-> FASTA parsing with SeqIO
-> File format conversion
-> Basic sequence manipulation
# 2. Biopython Practice Suite
Folder: /Biopython/
Goal: Learn sequence manipulation and pairwise alignment operations.
Highlights:
basic_biopython.py: Demonstrates sequence object operations (translation, reverse complement, slicing).
sequence_alignment.py: Performs global/local alignments using Biopython’s pairwise2 module.
# 3. Comparative Analysis of SNP Distribution and Gene Expression
Folder: /Comparative_Analysis_of_SNP_Distribution_and_Gene_Expression_Across_Diverse_Organism/
Goal: Explore the relationship between SNP density and gene expression levels across different organisms.
Methods:
Data integration from SNP and expression datasets
Visualization with Matplotlib and Seaborn
Output plots: snp_vs_expression.png, top_expression_bar.png, etc.
# 4. Bacterial Genome Analysis
Folder: /bacterial_genome_analysis/
Goal: Perform exploratory analysis of bacterial genomes.
Analyses Performed:
-> GC content and sequence length distribution
-> Restriction enzyme site mapping
-> HTML visualization of restriction data
Outputs: gc_content_bar.png, restriction_scatter.html
# 5. Drug Likeliness Prediction
Folder: /drug_likeliness_prediction/
Goal: Predict drug-likeness of compounds using RDKit descriptors.
Dataset: Drug_like.sdf vs Drug_decoys.sdf
Notebook: drug_likeliness_pred.ipynb
Techniques:
-> Molecular feature extraction
-> Classification using machine learning models
-> Evaluation metrics and feature importance visualization
# 6. Promoter Prediction Using Machine Learning
Folder: /promoter_prediction/
Goal: Predict promoter activity from DNA sequences.
Modules:
extract_features.py: Extracts GC content, k-mers, motifs (TATA, CAAT), CpG ratio, etc.
model.py: Builds regression model (e.g., Random Forest or Linear Regression).
main.py: Integrates feature extraction and prediction pipeline.
Outputs: features.csv, model metrics, and prediction reports.
# 7. SNV Classification
Folder: /snv_classification/
Goal: Classify Single Nucleotide Variants (SNVs) as deleterious or neutral.
Pipeline:
vcf_filter.py: Cleans and filters variant data
extract_features.py: Computes sequence-based and structural features
training.py & model.py: Train ML models (Logistic Regression, Random Forest)
evaluate.py & visualization.py: Evaluate using ROC and confusion matrix plots
Outputs: mapped_variants.csv, ROC and confusion matrix visualizations
