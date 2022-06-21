# leech_DE
Differential expression analysis of leech developmental stages

## Overview
1. Trimmomatic was run in paired end mode, with seed mismatches set to 2, palindrome clipping threshold set to 30 matches, simple clip threshold to 10 matches, and minimum adapter length to 1bp with the “keepBothPairs” option enabled. Leading and Trailing qualities were set to a value of 3. Sliding window size was set to 4bp, with a required quality of 15, and minimum read length was thresholded at 30bp.


2. STAR was run with default parameters. The output was sorted and indexed with samtools prior to further analysis.



3. Feature counts were run in non-stranded fashion to collect counts from the provided reference annotation (passed to us in Dan’s emails) to summarize features across gene features based on the corresponding gene_name or gene_id annotation in the reference annotation files (depending on analyzed species). The largest overlap option was turned on to ensure that reads were assigned to genes that they overlapped with best in the annotation.


4. While the DESeq2 code is not something I can share, please refer to the DESeq2 vignette here, which can provide you with a suitable idea of the manner of the code that was run.