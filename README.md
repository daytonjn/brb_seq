# Bulk RNA Barcoding RNA-seq (BRB-seq)

**Description:** Early development scripts to process multiplexed bulk 3'-end RNA-seq data with early inline barcoding (Poly(T) oligo w/ barcode incorporated during first strand synthesis), determine cycling transcripts/genes, and identify DEGs between treatments. This library preparation/sequencing strategy was previously benchmarked against TruSeq by Alpern et al. ([2019](https://doi.org/10.1186/s13059-019-1671-x)) and leveraged to quantify tissue-specific circadian gene expression across *Drosophila* Genetic Reference Panel lines (Litovchenko et al. [2021](https://doi.org/10.1126/sciadv.abc3781)). 

**Purpose:** Past research associated the ~2-3 week difference in the post-diapause development (PDD) time of the European corn borer (*Ostrinia nubilalis*) caterpillar with variation in the core circadian clock gene *period* (*per*) and a ~1 hr difference in the free-running period length of circadian activity. To quantify how/when differences in *per* genotype and the clock interact with the environment to contribute to these phenotypes, we extracted RNA from *O. nubilalis* brains & performed BRB-seq to identify DEGs & DR pathways. Samples from each *per* background were collected at the photosensitive stage of the 5<sup>th</sup> larval instar every 3 hr (12L:12D) for a 24 hr period and compared with samples collected at a single timepoint (1 hr before lights-on/"dawn") that were reared in other photoperiods (15L:09D; 12L:8D:2L:2D). There were 3-4 biological replicates/treatment, each containing an average of ~4 brains. \
\
Sequencing a large number of samples with biological replicates is required for sufficient power to quantify the effect(s) of genotype (G), environment (E), time (T), GxE, GxT, & ExT. Although traditional TruSeq approaches capture both mature mRNAs (polyadenylated transcripts) and the regulatory/non-coding transcripts too, TruSeq library preparation & sequencing is comparatively low-throughput and prohibitively expensive (compared to a 3'-end gene-counting approach); hence, this BRB-seq approach was used.

**BRB-seq library preparation workflow:**

<p align="center">
  <img src="https://github.com/user-attachments/assets/151ea57a-1641-4a43-a839-4184895ae8b1" alt="brb-seq workflow"/>
</p>

Figure from Alpern et al. (2019).

**Analysis Workflow:**
 1. Filter rRNA genes from GFF (alternative to removing any rRNA reads)
    - 01_remove_rrna_from_gff.sh (awk)
 2. Index *O. nubilalis* reference genome (GCF_963855985.1) using the no-rRNA GFF
    - 02_index_ref_norrna.sh (STAR)
 3. Align .fastq files to genome & quantify per-barcode gene expression
    - 03_map_reads_norrna.sh (STAR/STARsolo)
    - barcodes.txt
 4. Generate sample read count matrix in R & output files
    - 04_generate_read_counts.R
 5. Perform base RNAseq analysis for DEGs
     -analysis_brbseq.R (limma-voom & edgeR, limorhyde)

