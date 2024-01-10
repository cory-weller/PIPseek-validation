# PIPseek-validation

## Retrieve files from Google Drive
```bash
module load rclone
rclone copy --drive-shared-with-me "dtigoogledrive:/DR_cellranger test/" .
```

## Check gRNA distribution
Looking at `DR-9163_25_S25_L001_R2_001.fastq.gz`. See [count_gRNAs.R](scripts/count_gRNAs.R) for detail. Executed as:
```bash
Rscript scripts/count_gRNAs.R
```

Briefly, the script 
1. reads in the sequence-containing rows from the `fastq` file
2. excludes reads with `N`s (there aren't many anyway)
3. finds reads with expected pattern of `'GTTG[ACTG]{20}GTTT'` where the gRNA is precisely 20 nt sandwiched between `GTTG` and `GTTT`
4. matches exctracted gRNA sequence with those expected from `data/adt-tags-183sgRNAsequences-Tian2019.csv`
5. plots distribution of guides with at least 1 occurence

![](plots/gRNA_distribution_linear.png)

![](plots/gRNA_distribution_log10.png)

See table: [gRNA_distribution.tsv](output/gRNA_distribution.tsv)
