# pci_snATAC

This repository contains scripts used in Han et al. (2025) for analyzing single-nucleus ATAC-seq data, identifying regulatory elements, and integrating with genetic and functional genomic datasets.

All scripts were developed and tested under the following environments: R v4.2.2 and Python v3.9.7.

All dependencies can be installed using `BiocManager::install()` in R and pip in python.

Example: install ArchR and other dependencies
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ArchR", "Seurat", "limma", "enrichR", ...))  # Add more as needed
```
Each section below corresponds to major analysis steps in the manuscript. Scripts are organized by number for clarity.

# 0. Data Processing and Cell Type Annotation
Scripts for preprocessing scATAC-seq data, clustering, and annotating cell types and subtypes.

00_archrproject/

`01.qc.R`

`02.makeProject.R`

`03.umap_callpeak.R`

`04.subset_to_define_subtype.R`

# 1. Identification of cis-regulatory elements
Define and validate cis-regulatory elements and their target genes.

10_regulatory_element/

`11.subC_cA.R` : Identify peak-peak coaccessibility (CAs)

`12.subC_peak2gene.R` : Identify peak-to-gene relationships

`13.get.CAs_subG.R` : Extract CAs from computation result

`14.get.P2G_subG.R` : Extract P2G links from computation result

`15.classify.PEG.R` : Identify and summarize the enhancers

`16.ValidateCRE.R`: Validate CREs using external datasets


# 2. Identification of differentially accessible regions (DARs)
Determine disease vs. control DARs in cell types and subtypes.

20_dar_dgs/

`21.DAR_byRegion_celltype_discovery.R`

`22.DAR_byRegion_celltype_downsample.R`

`23.subC_markerPeaks.R`

`24.DifferentialGeneScores.R`


# 3. Identification of dynamic CREs
Identify CREs dynamically changed across four conditions and overlap with GWAS SNPs.

30_dynamicpeaks/

`31.dynamic_peaks_disease.R`

`32.dynamic_peaks_w_PSPGWASsnps_3types.R`


# 4. Heritability partition by LDSC
Scripts for LD score regression and heritability enrichment.

40_ldsc/

`41.extract_clusterSpecificPeak_for_LDSC.R`

`42.batch.make_annot.sh`

`43.batch.compute_LDscore.sh`

`44.batch.partition_h2.sh`

`45.analyze_h2Out_stauc.R`


# 5. Enrichment of MPRA and sn-eQTLs
Integration with functional genomic data to examine the distribution of MPRA variants and assess their enrichment in dynamic peaks.

50_MPRA_eQTL/

`02.Plink_LDblock.sh`

`03.Plink_calLD_r2.sh`

`04.CREs_TFBS_FIMO.R`

`51.MPRA_raw_fishertest.R`

`52_MPRA_raw_in_genes.R`

`53.MPRA_effect_distri.R`

`54.MPRA_fsVar_peak_distri.R`

`55.CREs_w_frVar_Module-subC.R`

`56.MPRA_frVar_motif_discover.R`

`57.QTL_cal_fishertest_eGene_match.R`


# 6. Regulatory TFs and target genes
Footprinting analysis and TF-target inference using TOBIAS.

60_TFBS_TOBIAS/

`00.get_subc_barcodes.R`

`01.Bam_to_Sam_86samples.sh`

`02.Filter_Bam_perSubC_perSample.sh`

`03.Merge_Bam.sh`

`04.TOBIAS_ATACorrect.sh`

`05.TOBIAS_footprint.sh`

`06.TOBIAS_BINDetect.sh`

`meme.R`


# 7. Epigenomic stability analysis with single-nucleus methylation data
Integration with single-nucleus methylation (snmC-seq) data to assess the epigenomic stability of cell types and regulatory regions.

70_sn3mc_allcools/

`71.loop_mcds.py`

`72.extract_mcds.py`

`73.allcools_region.sh`

`74.allcools_generate_dataset.sh`

`75.sbatchPermute.R`

`76.Permutation_byGroup.R`


# 8. Cell-cell interactions analysis
Analyze ligand–receptor interactions using CellPhoneDB.

80_cellphoneDB/

`80.DownloadDB.py`

`81.DEG_seurat.R`

`82.Run_statistical_analysis.py`

`83.Run_deg_analysis.py`

`84.sumNb_degCCI_subC.R`

# 9. Topic modeling of variant-linked CREs in microglia
Topic modeling of microglial CREs containing either frVars or sn-eQTLs was performed to determine whether variant-linked regulatory elements converge into coherent and biologically meaningful regulatory programs associated with disease-related microglial states.

90_topic_modeling_variant_linked_CREs_microglia/

`91.extract_PeakMat.R`

`92.extract_CREset_w_SNPs.R`

`93.run_cisTopic.R`

`94.topics_plot_TF_enrichR.R`

`95.sum_topicst_snps-ADgwas.R`

`95.sum_topicst_snps-FTDgwas.R`

`96.topics_refine_by_frVar.R`
