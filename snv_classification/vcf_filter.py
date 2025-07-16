import pandas as pd
from sklearn.utils import resample

def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    return info_dict

vcf_path = "clinvar.vcf"  # or "clinvar.vcf.gz"
variants = []

# Step 1: Parse VCF and extract SNVs with Likely_benign or Likely_pathogenic
with open(vcf_path, 'r') as f:
    for line in f:
        if line.startswith("#"):
            continue

        cols = line.strip().split('\t')
        if len(cols) < 8:
            continue

        chrom, pos, var_id, ref, alt, qual, filt, info = cols[:8]
        info_dict = parse_info(info)

        clnsig = info_dict.get('CLNSIG', '')
        clnvc = info_dict.get('CLNVC', '').lower()

        if clnsig not in ['Likely_benign', 'Likely_pathogenic']:
            continue
        if clnvc != 'single_nucleotide_variant':
            continue

        label = 1 if clnsig == 'Likely_pathogenic' else 0
        gene = info_dict.get('GENEINFO', '').split(':')[0]
        variant_type = info_dict.get('MC', 'NA').split('|')[0]

        variants.append({
            'Chrom': chrom,
            'Pos': int(pos),
            'Ref': ref,
            'Alt': alt,
            'Gene': gene,
            'Variant_Type': variant_type,
            'Label': label
        })

# Step 2: Convert to DataFrame
df = pd.DataFrame(variants)

# Step 3: Balance the classes to 10,000 each
df_benign = df[df['Label'] == 0]
df_pathogenic = df[df['Label'] == 1]

# Make sure there are at least 10,000 of each class
min_count = min(len(df_benign), len(df_pathogenic), 10000)

df_benign_sampled = resample(df_benign, replace=False, n_samples=min_count, random_state=42)
df_pathogenic_sampled = resample(df_pathogenic, replace=False, n_samples=min_count, random_state=42)

df_balanced = pd.concat([df_benign_sampled, df_pathogenic_sampled]).sample(frac=1, random_state=42).reset_index(drop=True)

# Step 4: Save to CSV
df_balanced.to_csv("snv.csv", index=False)

# Step 5: Summary
print("Saved balanced dataset as 'snv_balanced_10000_each.csv'")
print(df_balanced['Label'].value_counts())
print(df_balanced.head())
