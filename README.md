# PIPseek-validation

## Retrieve files from Google Drive
```bash
module load rclone
rclone copy --drive-shared-with-me "dtigoogledrive:/DR_cellranger test/" .
```

## Run basic cellranger
Rscript scripts/count_gRNAs.R

```R

```

![](plots/gRNA_distribution_linear.png)

![](plots/gRNA_distribution_log10.png)