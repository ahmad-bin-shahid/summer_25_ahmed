import pandas as pd  # For data manipulation and analysis
import matplotlib.pyplot as plt  # For creating plots and visualizations
from scipy.stats import pearsonr  # For calculating Pearson correlation coefficient
from Bio import SeqIO  # For biological sequence input/output operations
from collections import defaultdict  # For creating dictionaries with default values

def parse_gff_file(gff_path, max_exons=5000):
    """Parse GFF file to extract exon information"""
    exons = []  # Initialize empty list to store exon data
    try:
        with open(gff_path, 'r') as file:  # Open the GFF file for reading
            for line in file:  # Iterate through each line in the file
                if line.startswith("#"):  # Skip comment lines that start with #
                    continue
                parts = line.strip().split('\t')  # Split line by tab characters
                if len(parts) != 9:  # Check if line has exactly 9 columns (GFF format)
                    continue
                # Extract the 9 standard GFF columns
                chrom, source, feature_type, start, end, score, strand, phase, attributes = parts
                if feature_type != "exon":  # Only process exon features
                    continue

                chrom = chrom.strip()  # Remove whitespace from chromosome name
                if chrom == "AE014298.5":  # Convert specific chromosome identifier
                    chrom = "chr4"
                elif chrom.startswith("AE"):  # Keep AE prefixed chromosomes as is
                    chrom = chrom
                elif not chrom.startswith("chr"):  # Add chr prefix if missing
                    chrom = "chr" + chrom

                gene_id = None  # Initialize gene ID variable
                for field in attributes.split(';'):  # Split attributes by semicolon
                    if field.startswith("Dbxref="):  # Look for database cross-reference
                        values = field.split('=')[1].split(',')  # Extract values after =
                        for v in values:  # Check each value
                            if v.startswith("FLYBASE:FBgn"):  # Look for FlyBase gene ID
                                gene_id = v.split(':')[1]  # Extract gene ID after colon
                                break

                if gene_id:  # Only add exon if we found a gene ID
                    exons.append({  # Add exon data to list
                        "GeneID": gene_id,
                        "Chrom": chrom,
                        "Start": int(start),  # Convert start position to integer
                        "End": int(end)  # Convert end position to integer
                    })

                if len(exons) >= max_exons:  # Stop if we've reached the maximum exons
                    break

        return exons  # Return the list of exon data
    except FileNotFoundError:  # Handle case where file doesn't exist
        print(f"Error: File not found - {gff_path}")
        return []

def parse_vcf_data(vcf_file_path):
    """Parse VCF file to extract SNP positions"""
    snps = []  # Initialize empty list to store SNP data
    try:
        with open(vcf_file_path, 'r') as file:  # Open VCF file for reading
            for line in file:  # Iterate through each line
                if line.startswith("#") or not line.strip():  # Skip header lines and empty lines
                    continue
                parts = line.strip().split('\t')  # Split line by tab characters
                if len(parts) < 2:  # Check if line has at least 2 columns
                    continue
                chrom = parts[0].strip()  # Extract chromosome from first column
                pos = parts[1]  # Extract position from second column
                if pos.replace('.', '', 1).isdigit():  # Check if position is numeric
                    pos = int(float(pos))  # Convert position to integer
                    snps.append((chrom, pos))  # Add chromosome and position as tuple
        return snps  # Return list of SNP tuples
    except FileNotFoundError:  # Handle case where file doesn't exist
        print(f"Error: File not found - {vcf_file_path}")
        return []

def group_snps_by_chrom(snps):
    """Group SNPs by chromosome for efficient lookup"""
    grouped = defaultdict(list)  # Create dictionary with list as default value
    for chrom, pos in snps:  # Iterate through each SNP tuple
        grouped[chrom].append(pos)  # Add position to the chromosome's list
    return grouped  # Return dictionary with chromosome as key, positions as values

def map_snps_to_exons(exons, snp_dict):
    """Map SNPs to exons and count SNPs per gene"""
    gene_snp_counts = {}  # Dictionary to store gene-level SNP counts
    for exon in exons:  # Iterate through each exon
        gene = exon["GeneID"]  # Extract gene ID
        chrom = exon["Chrom"]  # Extract chromosome
        start = exon["Start"]  # Extract start position
        end = exon["End"]  # Extract end position
        exon_length = end - start  # Calculate exon length
        count = 0  # Initialize SNP count for this exon

        if chrom in snp_dict:  # Check if chromosome has SNPs
            for snp_pos in snp_dict[chrom]:  # Iterate through SNPs on this chromosome
                if start <= snp_pos <= end:  # Check if SNP is within exon boundaries
                    count += 1  # Increment SNP count

        if gene not in gene_snp_counts:  # If gene not yet in dictionary
            gene_snp_counts[gene] = {"SNPs": 0, "Length": 0}  # Initialize gene entry
        gene_snp_counts[gene]["SNPs"] += count  # Add SNP count to gene total
        gene_snp_counts[gene]["Length"] += exon_length  # Add exon length to gene total
    return gene_snp_counts  # Return dictionary with gene-level data

def calculate_density(gene_snp_counts):
    """Calculate SNP density per kilobase for each gene"""
    data = []  # Initialize list to store calculated data
    for gene, values in gene_snp_counts.items():  # Iterate through each gene
        snp_count = values["SNPs"]  # Extract SNP count
        length = values["Length"]  # Extract total exon length
        length_kb = length / 1000 if length > 0 else 0  # Convert length to kilobases
        density = snp_count / length_kb if length_kb > 0 else 0  # Calculate density
        data.append({  # Add gene data to list
            "GeneID": gene,
            "SNP_Count": snp_count,
            "Exon_Length": length,
            "SNP_Density": density
        })
    return pd.DataFrame(data)  # Convert list to pandas DataFrame

def merge_with_expression(snp_df, expression_file):
    """Merge SNP data with gene expression data"""
    try:
        expr_df = pd.read_csv(expression_file, sep='\t')  # Read expression file
        expr_df.rename(columns={expr_df.columns[0]: "GeneID"}, inplace=True)  # Rename first column
        expr_df = expr_df[['GeneID', 'F_Brain']]  # Select only gene ID and brain expression
        expr_df.rename(columns={'F_Brain': 'Expression'}, inplace=True)  # Rename expression column
        merged = pd.merge(snp_df, expr_df, on="GeneID", how="inner")  # Merge datasets on gene ID
        return merged  # Return merged DataFrame
    except FileNotFoundError:  # Handle case where expression file doesn't exist
        print(f"Error: Expression file not found - {expression_file}")
        return pd.DataFrame()  # Return empty DataFrame

def perform_pearson_correlation(df):
    """Calculate Pearson correlation between SNP density and expression"""
    try:
        x = df['SNP_Density']  # Extract SNP density values
        y = df['Expression']  # Extract expression values
        r, p = pearsonr(x, y)  # Calculate Pearson correlation coefficient and p-value
        print("\nPearson Correlation Results:")
        print(f"Correlation coefficient (r): {r:.4f}")
        print(f"P-value: {p:.4e}")
        return r, p  # Return correlation coefficient and p-value
    except Exception as e:  # Handle any calculation errors
        print("Error calculating Pearson correlation:", e)
        return None, None

def plot_scatter(df, output_file="snp_vs_expression.png"):
    """Create scatter plot of SNP density vs expression, removing outliers"""
    try:
        # Calculate 95th percentile cutoffs to remove extreme outliers
        expr_limit = df['Expression'].quantile(0.95)
        snp_limit = df['SNP_Density'].quantile(0.95)
        
        # Filter data to remove top 5% of values in both dimensions
        filtered_df = df[(df['Expression'] <= expr_limit) & (df['SNP_Density'] <= snp_limit)]
        
        plt.figure(figsize=(8,6))  # Create figure with specified size
        plt.scatter(filtered_df['SNP_Density'], filtered_df['Expression'], alpha=0.5)  # Create scatter plot
        plt.title("SNP Density vs Gene Expression")  # Add title
        plt.xlabel("SNP Density (SNPs per kb of exons)")  # Add x-axis label
        plt.ylabel("Gene Expression")  # Add y-axis label
        plt.grid(True)  # Add grid for better readability
        
        r, _ = pearsonr(filtered_df['SNP_Density'], filtered_df['Expression'])  # Calculate correlation
        plt.annotate(f'r = {r:.2f}', xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12)  # Add correlation annotation
        
        plt.tight_layout()  # Adjust layout to prevent label cutoff
        plt.savefig(output_file)  # Save plot to file
        print(f"Scatter plot saved as {output_file}")
        print(f"Removed {len(df) - len(filtered_df)} outlier genes")
        
    except Exception as e:  # Handle any plotting errors
        print("Error creating scatter plot:", e)

def plot_histogram(df, column='SNP_Density', output_file="snp_density_hist.png"):
    """Create histogram of specified column, removing outliers"""
    try:
        # Calculate 95th percentile cutoff to remove extreme outliers
        limit = df[column].quantile(0.95)
        filtered_data = df[df[column] <= limit][column]  # Filter data below cutoff
        
        plt.figure(figsize=(8,6))  # Create figure with specified size
        plt.hist(filtered_data, bins=30, alpha=0.7, color='skyblue', edgecolor='black')  # Create histogram
        plt.title(f"Histogram of {column}")  # Add title
        plt.xlabel(column)  # Add x-axis label
        plt.ylabel("Frequency")  # Add y-axis label
        plt.grid(True)  # Add grid for better readability
        plt.tight_layout()  # Adjust layout to prevent label cutoff
        plt.savefig(output_file)  # Save plot to file
        print(f"Histogram saved as {output_file}")
        print(f"Removed {len(df) - len(filtered_data)} outlier values")
        
    except Exception as e:  # Handle any plotting errors
        print(f"Error creating histogram for {column}:", e)

def plot_top_expression_bar(df, top_n=20, output_file="top_expression_bar.png"):
    """Create bar plot of top expressed genes, removing outliers first"""
    try:
        # Calculate 95th percentile cutoff to remove extreme outliers
        expr_limit = df['Expression'].quantile(0.95)
        filtered_df = df[df['Expression'] <= expr_limit]  # Filter data below cutoff
        
        # Sort by expression and take top N genes from filtered data
        top_df = filtered_df.sort_values(by='Expression', ascending=False).head(top_n)
        
        plt.figure(figsize=(10,6))  # Create figure with specified size
        plt.bar(top_df['GeneID'], top_df['Expression'], color='green')  # Create bar plot
        plt.xticks(rotation=90)  # Rotate x-axis labels for better readability
        plt.title(f"Top {top_n} Expressed Genes")  # Add title
        plt.xlabel("GeneID")  # Add x-axis label
        plt.ylabel("Expression Level")  # Add y-axis label
        plt.tight_layout()  # Adjust layout to prevent label cutoff
        plt.savefig(output_file)  # Save plot to file
        print(f"Bar plot saved as {output_file}")
        print(f"Removed {len(df) - len(filtered_df)} outlier genes before selecting top genes")
        
    except Exception as e:  # Handle any plotting errors
        print("Error creating bar plot:", e)

def main():
    """Main function to orchestrate the entire analysis pipeline"""
    # Define input file paths
    gff_file = "genomic.gff"
    vcf_file = "output_snps.vcf"
    expression_file = "expression_matrix.tsv"

    # Parse exon data from GFF file
    exons = parse_gff_file(gff_file, max_exons=5000)
    print(f"Parsed {len(exons)} exons from GFF")
    
    # Parse SNP data from VCF file
    snps = parse_vcf_data(vcf_file)
    print(f"Parsed {len(snps)} SNPs from VCF file")
    
    # Group SNPs by chromosome for efficient lookup
    grouped_snps = group_snps_by_chrom(snps)
    
    # Map SNPs to exons and calculate per-gene counts
    gene_snp_counts = map_snps_to_exons(exons, grouped_snps)
    
    # Calculate SNP density per kilobase for each gene
    snp_df = calculate_density(gene_snp_counts)
    
    # Merge SNP data with expression data
    merged_df = merge_with_expression(snp_df, expression_file)

    if not merged_df.empty:  # Check if merged data contains any records
        output_csv = "snp_expression_output.csv"
        try:
            merged_df.to_csv(output_csv, index=False)  # Save merged data to CSV
            print(f"Output saved to {output_csv}")
        except Exception as e:  # Handle any file writing errors
            print("Error saving output file:", e)

        # Perform statistical analysis and create visualizations
        perform_pearson_correlation(merged_df)
        plot_scatter(merged_df)
        plot_histogram(merged_df)
        plot_top_expression_bar(merged_df)
    else:  # Handle case where no data was successfully merged
        print("Merged data is empty. Check gene ID alignment or input content.")

if __name__ == "__main__":  # Check if script is run directly (not imported)
    main()  # Execute main function