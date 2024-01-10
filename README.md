# PIPseek-validation

## Retrieve files from Google Drive
```bash
module load rclone
rclone copy --drive-shared-with-me "dtigoogledrive:/DR_cellranger test/" .
```

## Count gRNA distribution
see [count_gRNAs.R](scripts/count_gRNAs.R) for detail. Executed as:
```bash
Rscript scripts/count_gRNAs.R
```

![](plots/gRNA_distribution_linear.png)

![](plots/gRNA_distribution_log10.png)