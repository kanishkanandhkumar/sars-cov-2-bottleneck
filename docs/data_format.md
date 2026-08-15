# Data Format Documentation

## Input Data Format

The analysis expects a CSV file with the following columns:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| pair_id | string | Transmission pair identifier | "Pair_01" |
| variant_id | string | Variant identifier | "chr1_12345_A_G" |
| donor_freq | float | Donor minor allele frequency (0-1) | 0.15 |
| recipient_depth | int | Recipient sequencing depth at locus | 1500 |
| recipient_count | int | Recipient minor allele count | 23 |
| transmitted | int | 1 if transmitted, 0 if absent | 1 |

## Example Data

```csv
pair_id,variant_id,donor_freq,recipient_depth,recipient_count,transmitted
Pair_01,chr1_12345_A_G,0.15,1500,23,1
Pair_01,chr1_67890_C_T,0.08,1450,5,1
Pair_02,chr2_23456_G_A,0.25,1800,60,1

##Quality Control Metrics

#Required QC metrics for each pair:

Donor coverage ≥ 1,000×

Recipient coverage ≥ 1,000×

Mapping quality ≥ 30

Base quality ≥ 30

Strand bias p-value ≥ 0.05
