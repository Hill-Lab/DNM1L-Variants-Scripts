# DNM1L-Variants-Scripts
R-scripts &amp; ImageJ macros used in the submitted manuscript "Novel DNM1L variants impair mitochondrial dynamics through divergent mechanisms"

# **R-scripts for biophysical studies**:
Contributed by K. Nolden

1. 20200417_Drp1_variant_additional_boxplots
* Calculating and plotting Kcat,K0.5, Kcat/K0.5, and vmax for the following: 
* Drp1 WT, L230dup, G363D, G401S, and R710G GTPase results
* Experiment initially performed with recombinant L230dup
* However, later MS + SDS-PAGE confirmed protein was truncated and not usable

2. 20210608_Drp1_variants_analysis.R
* Script for analyzing and visualizing SEC-MALS data
* Data collected on 20210608
* Drp1 constructs dialyzed overnight into SEC-MALS buffer
* Exact buffer match to running buffer used
* 6xHis tag present

3. 20210804_CD_Drp1_UK.R
* Analyzing circular dichroism data
* experiment: CD of WT Drp1 + UK pathological variants
* data collected on 20210804
* data corrected for empty cuvette + buffer signal
* scaled in Excel prior to loading to provide comparable baselines
* MRE calculated for each sample at each wavelength evaluated

4. UK_Drp1_WT_variants_TSA_merged.R
* Plotting combined Thermal Shift Assay (TSA) data
* Amplification data = Temperature vs Fluorescence
* Sypro-Orange thermofluor TSA on Drp1 WT and UK variants +/- 100 or 500 uM GDP or GTP
* Data collected on 20210514 and 20210602
* 1st derivative of fluorescence with respect to temperature used to determine
* melting temperature (Tm)

# **R-scripts for analysing co-localization data**:

1. Drp1_coloc_analysis.R
* Extract Pearson's R values from ImageJ Output
* for DRP1-TOM20 colocalization experiments using patient-derived
* fibroblasts immunostained for DRP1, TOM20, and DAPI
* performed by Ollie of UK team
* Two data sets per sample: C1, C2, P1, P2, P3, and P4
* C1: pediatric control
* C2: adult control
* P1: Drp1 G401S variant
* P2: Drp1 G363D variant
* P3: Drp1 L230dup variant
* P4: Drp1 R710G variant

2. Drp1_coloc_analysis_PMP.R
* Extract Pearson's R values from ImageJ Output
* for DRP1-PMP70 colocalization experiments using patient-derived
* fibroblasts immunostained for DRP1, TOM20, and DAPI
* performed by Ollie of UK team
* Two replicates per sample: C1, C2, P1, P2, P3, and P4
* Collected in same experiment/data set
* C1: pediatric control
* C2: adult control
* P1: Drp1 G401S variant
* P2: Drp1 G363D variant
* P3: Drp1 L230dup variant
* P4: Drp1 R710G variant
