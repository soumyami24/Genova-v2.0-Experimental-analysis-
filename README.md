# Genova: Experimental Validation, Statistical Benchmarking & Feature Ablation on South Asian Cohorts

> **Note on Intellectual Property:** The core multi-ancestry architecture, weight optimization algorithms, and parameter fusion mechanics underlying **Genova** are currently **patent-pending**. Consequently, raw model weights and proprietary aggregation scripts are omitted from this public repository. This repository provides the complete experimental validation matrices, mathematical derivations of the evaluation metrics, and feature-ablation spreadsheets demonstrating the model's empirical behavior on South Asian populations.
> **Recently we change the name of Genova-v2.0 to SAKSHAM-PRS so in some files instead of Genova,SAKSHAM-PRS is mentioned.**
---

##  Executive Summary & Research Motivation

Standard Polygenic Risk Score (PRS) models trained predominantly on European (EUR) discovery cohorts exhibit severe **ancestry based biasness** when deployed on South Asian (SAS) genomes due to differences in allele frequencies and linkage disequilibrium (LD) structure. **Genova** was engineered to address this cross-ancestry generalization gap for ischemic stroke risk,diabetes and other polygenic diseases stratification while rigorously controlling for sub-continental population substructure.

This repository documents the **in-silico validation, statistical mechanics, and targeted feature ablation** of the Genova inference engine across **n = 182** South Asian individuals from the 1000 Genomes Phase 3 dataset, benchmarked against a pan-South Asian reference panel (**N = 489**) and an external `PRS-CS` (European reference) baseline.Here we did this for ischemic stroke mainly.

### Main Results

* **64.98% Reduction in Sub-Population Drift:** Across distinct South Asian subpopulations-Punjabi from Lahore (`PJL`, n = 96) and Bengali from Bangladesh (`BEB`, n = 86)-Genova compresses inter-cohort mean shift from **0.5368 SD** (observed in `PRScs-European`) down to **0.1880 SD**.
* **56.26% Variance Regularization:** In the combined South Asian cohort (n = 182), Genova reduces genome-wide score variance from **1.0818** (SD = 1.0401 in `PRScs-European`) to **0.4732** (SD = 0.6879), alongside a **37.53% compression in interquartile range (IQR)** (0.8609 vs. 1.3782), eliminating extreme tail-risk distortion.
* **Sub-Population Feature Invariance (Delta = 0.000012):** The marginal contribution of 6 literature-validated South Asian stroke SNPs remains virtually invariant between `PJL` (Mean = -0.000657) and `BEB` (Mean = -0.000645).
* **Polygenic Robustness & Ablation Stability (p = 0.4424):** Exact leave-out perturbation of the 6 known stroke loci induces a bounded mean absolute shift of **0.1494%** of the genome-wide standard deviation (max = **0.5448%**), with paired Wilcoxon signed-rank testing (Z = -0.7682, p = 0.4424) confirming zero systematic directional bias or oligogenic collapse.
* **Mathematical & Algorithmic Verification (r = 0.99655):** Variant-level risk contributions between Genova and `PRScs-European` achieve **R² = 99.31%** concordance across all 182 individuals, verifying the exactness of the custom allele-harmonization and dosage-scoring engine, while the genome-wide scores share **R² = 76.56%** variance (r = 0.875)-leaving **23.44%** orthogonal variance driven by Genova's ancestry-aware calibration.

---

## ⚙️Bioinformatics Pipeline designed by me

Evaluating genome-wide polygenic models across multi-cohort whole-genome sequencing (WGS) data required building a fault-tolerant data engineering pipeline in `Python` (`pandas`, `numpy`, `scipy`) and `Bash` (`PLINK 2.0`, `GNU Awk`):

1. **Genotype Matrix Conversion & Coordinate Indexing:**
   * Extracted South Asian cohorts (`PJL` n = 96, `BEB` n = 86, and full reference N = 489) from 1000 Genomes Phase 3 VCFs into binary `PLINK` sets (`.bed`, `.bim`, `.fam`) partitioned across chromosomes 1 to 22.
   * Standardized internal variant identifiers into canonical `CHR:POS:REF:ALT` strings to bypass inconsistent `rsID` builds between `GRCh37`/`GRCh38` weight matrices and target genotypes.

2. **Automated Allele Harmonization & Strand-Flip Resolution:**
   * Engineered a high-speed chromosome-matching engine that joins multi-million-row posterior weight files against `.bim` coordinate tables on `(CHR, POS)`.
   * Automatically classified every locus into three states:
     * `SAME`: Exact orientation match (`A1 == G_A1` and `A2 == G_A2`).
     * `REVERSED`: Ref/Alt allele swap (`A1 == G_A2` and `A2 == G_A1`), preserving true effect-allele mapping without corrupting beta weight signs.
     * `MISMATCH`: Incompatible tri-allelic or strand-ambiguous variants, which were programmatically quarantined and filtered prior to scoring.

3. **Batch Chromosomal Linear Scoring & Reference Panel Standardization:**
   * Executed parallel `PLINK2 --score` runs (`cols=+scoresums`) across 22 autosomal chromosomes per cohort and aggregated chromosomal allelic dosages:
     ```text
     PRS_raw(i) = SUM_{chr = 1..22} [ SUM_{snp in chr} ( Dosage(i, snp) * Beta(snp) ) ]
     ```
   * Rather than normalizing `PJL` and `BEB` independently (which artificially masks population stratification by forcing both cohorts to Mean = 0, SD = 1), raw scores were standardized against population moments (`Mean_ref`, `SD_ref`) computed over the full **N = 489** South Asian super-population:
     ```text
     Z_i = ( PRS_raw(i) - Mean_ref_489 ) / SD_ref_489
     ```

4. **Targeted Sub-Pathway Scoring (`WithVariants` Isolation):**
   * Built a dedicated `subprocess` orchestration module targeting chromosomes `{1, 2, 5, 8, 19}` to isolate the exact marginal dosage contribution of curated stroke polymorphisms on the same reference scale (`SD_ref_489`), enabling exact additive decomposition.

---

##  Snp

To verify whether Genova captures clinically established South Asian ischemic stroke biology while maintaining a well-regularized polygenic architecture, 8 candidate single-nucleotide polymorphisms (SNPs) were curated from South Asian stroke meta-analyses and case-control studies (*Yadav et al., 2013; Saleheen et al., 2005; Munshi et al., 2012*).

Following automated intersection against our harmonized chromosomal weight matrices and `.bim` files, **6 primary functional variants** were retained and 2 were excluded due to absence in reference discovery weights:

| Gene | Variant | rsID | Chr:Pos (GRCh37) | Genotype ID (CHR:POS:REF:ALT) | Case-Control Evidence in SAS | Biological Mechanism | Pipeline Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **PDE4D** | SNP83 | `rs966221` | `5:59502520` | `5:59502520:A:G` | OR = 2.20 (p = 0.01, N = 1338); Pakistani OR = 1.64 (p < 0.005) | Regulates intracellular cAMP degradation & vascular smooth muscle proliferation | **Retained** (`FOUND` in PJL & BEB) |
| **IL10** | -1082 G/A | `rs1800896` | `1:206946897` | `1:206946897:T:C` | OR = 1.44 (p = 0.01, N = 1414 Indian meta-analysis) | Anti-inflammatory cytokine suppressing pro-inflammatory vascular cascades | **Retained** (`FOUND` in PJL & BEB) |
| **MTHFR** | C677T | `rs1801133` | `1:11856378` | `1:11856378:G:A` | OR = 2.50 (p = 0.08); +5.18 umol/L homocysteine (p = 0.001) | Folate/homocysteine metabolism & endothelial oxidative stress | **Retained** (`FOUND` in PJL & BEB) |
| **CYP4F2** | G1347A | `rs2108622` | `19:15990431` | `19:15990431:C:T` | OR = 1.58 (p = 0.004, 507 cases / 487 controls in India) | Arachidonic acid / 20-HETE synthesis & cerebral vascular tone | **Retained** (`FOUND` in PJL & BEB) |
| **CYP11B2** | -344 C/T | `rs1799998` | `8:143999600` | `8:143999600:A:G` | OR = 0.43 (p = 0.0001, 403 cases / 394 controls) | Aldosterone synthase activity, RAAS pathway & hypertension | **Retained** (`FOUND` in PJL & BEB) |
| **IL1A** | -889 C/T | `rs1800587` | `2:113542960` | `2:113542960:G:A` | OR = 0.53 (p = 0.009, 176 cases / 212 controls) | Pro-inflammatory interleukin-1 alpha signaling in vascular endothelium | **Retained** (`FOUND` in PJL & BEB) |
| **PDE4D** | SNP32 | `rs456009` | - | - | OR = 1.57 (p = 0.045, monomorphic in South Indian cohorts) | Phosphodiesterase 4D signaling | *Excluded (Absent in reference weights)* |
| **ABCB1** | C3435T | `rs1045642` | - | - | OR = 2.66 (p < 0.000001, 560 cases / 560 controls) | ATP-binding cassette P-glycoprotein transporter | *Excluded (Absent in reference weights)* |

---

##  Statistical work and Results

The uploaded Excel workbooks (`COMBINED_experiment_results.xlsx` and sub-cohort sheets) contain the complete individual-level results. Below is the detailed mathematical breakdown and verification of every analytical metric.

### 1. Distributional Moments, Variance Compression & Drift Analysis

We evaluated the six core distributional parameters-**Mean**, **Median (Q50)**, **Standard Deviation (SD)**, **Variance (SD²)**, **Interquartile Range (IQR = Q75 - Q25)**, **Minimum**, and **Maximum**-across **PJL** (n = 96), **BEB** (n = 86), and the **Combined South Asian Cohort** (n = 182):

| Cohort | Scoring Variable | Mean | Median (Q50) | Std Dev (SD) | Variance (SD²) | IQR (Q75 - Q25) | Min | Max | Dynamic Range |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **PJL (n = 96)** | **Genova (Full)** | -0.065849 | -0.081760 | 0.676538 | 0.457704 | 0.851331 | -2.431873 | +1.468526 | 3.900399 |
| **PJL (n = 96)** | **PRScs-European (Full)** | +0.044803 | -0.005707 | 1.061450 | 1.126676 | 1.425180 | -2.768420 | +2.098660 | 4.867080 |
| **PJL (n = 96)** | *Genova 6-SNP (Delta Z)* | -0.000657 | +0.000094 | 0.001517 | 0.00000230 | 0.003000 | -0.003580 | +0.001074 | 0.004654 |
| **PJL (n = 96)** | *PRScs-EUR 6-SNP (Delta Z)* | -0.001287 | +0.000000 | 0.002741 | 0.00000751 | 0.005528 | -0.006659 | +0.001916 | 0.008575 |
| **BEB (n = 86)** | **Genova (Full)** | -0.253812 | -0.256959 | 0.690865 | 0.477294 | 0.853256 | -2.896791 | +1.238439 | 4.135230 |
| **BEB (n = 86)** | **PRScs-European (Full)** | -0.491975 | -0.495961 | 0.943760 | 0.890683 | 1.131715 | -2.921190 | +1.326170 | 4.247360 |
| **BEB (n = 86)** | *Genova 6-SNP (Delta Z)* | -0.000645 | +0.000097 | 0.001505 | 0.00000226 | 0.000697 | -0.003748 | +0.000760 | 0.004508 |
| **BEB (n = 86)** | *PRScs-EUR 6-SNP (Delta Z)* | -0.001179 | +0.000051 | 0.002682 | 0.00000719 | 0.001481 | -0.006536 | +0.001598 | 0.008134 |
| **Combined (n = 182)** | **Genova (Full)** | -0.154667 | -0.138529 | 0.687917 | 0.473230 | 0.860920 | -2.896791 | +1.468526 | 4.365317 |
| **Combined (n = 182)** | **PRScs-European (Full)** | -0.208839 | -0.221578 | 1.040115 | 1.081839 | 1.378208 | -2.921190 | +2.098660 | 5.019850 |
| **Combined (n = 182)** | *Genova 6-SNP (Delta Z)* | -0.000651 | +0.000094 | 0.001507 | 0.00000227 | 0.002885 | -0.003748 | +0.001074 | 0.004822 |
| **Combined (n = 182)** | *PRScs-EUR 6-SNP (Delta Z)* | -0.001236 | +0.000007 | 0.002706 | 0.00000732 | 0.005427 | -0.006659 | +0.001916 | 0.008575 |

#### Quantitative Calculations & ML Inference:

1. **Exact Sample-Weighted Mean Verification:**
   Combining the two sub-cohorts (`n_PJL = 96`, `n_BEB = 86`, `N = 182`) yields exact additive consistency:
   ```text
   Mean_Combined(Genova) = [ 96 * (-0.065849) + 86 * (-0.253812) ] / 182
                         = ( -6.321504 - 21.827832 ) / 182
                         = -0.154667

   Mean_Combined(PRScs)  = [ 96 * (0.044803) + 86 * (-0.491975) ] / 182
                         = ( 4.301088 - 42.309850 ) / 182
                         = -0.208839
   ```

2. **64.98% Reduction in Sub-Continental Distributional Shift:**
   When a risk model is deployed clinically across South Asia, geographic drift between Northwestern (Punjabi/Pakistani) and Eastern (Bengali) subpopulations can cause severe false-positive misclassifications if one group's baseline distribution is artificially shifted.
   ```text
   Baseline Drift (PRScs-European) = | +0.044803 - (-0.491975) | = 0.536778 SD units
   Genova Drift                    = | -0.065849 - (-0.253812) | = 0.187963 SD units

   Relative Drift Reduction        = ( 1 - 0.187963 / 0.536778 ) * 100%
                                   = 64.98%
   ```
   **Inference:** Genova removes nearly two-thirds of the inter-population baseline drift between `PJL` and `BEB`.

3. **56.26% Variance Compression & 37.53% IQR Tightening:**
   ```text
   SD Compression (Combined)       = ( 1 - 0.687917 / 1.040115 ) * 100% = 33.86%
   Variance Compression (Combined) = ( 1 - 0.473230 / 1.081839 ) * 100% = 56.26%
   IQR Compression (Combined)      = ( 1 - 0.860920 / 1.378208 ) * 100% = 37.53%
   ```
   **Cross-Group Homoscedasticity:** The standard deviation of Genova between `PJL` (SD = 0.676538) and `BEB` (SD = 0.690865) differs by only **0.014327** (2.07%), whereas `PRScs-European` exhibits a **0.117690** (11.09%) SD disparity between the two populations.

---

### 2. Targeted Feature Ablation & Sensitivity Bounds

A critical requirement in high-dimensional genomic modeling is ensuring that the predictor represents a truly distributed polygenic signal rather than suffering from **oligogenic collapse** (where a handful of candidate SNPs dominate the prediction).

To test this, we isolated the standardized contribution `Delta_Z_variants(i)` of the 6 curated stroke SNPs on the exact same reference standard deviation scale (`SD_ref_489`) used for the full genome-wide score, and computed the counterfactual **Ablated Score** for each individual `i = 1..182`:

```text
Genova_ablated(i)  = Genova_full(i) - Delta_Z_variants(i)
Ablation_Change(i) = Genova_full(i) - Genova_ablated(i) = Delta_Z_variants(i)
Absolute_Change(i) = | Genova_full(i) - Genova_ablated(i) | = | Delta_Z_variants(i) |
```

#### Step-by-Step Perturbation Calculations:

1. **Relative Mean Absolute Perturbation:**
   Across the combined cohort (n = 182), the mean absolute change is `0.0010280`. Expressed as a percentage of the full genome-wide standard deviation (`SD_Genova = 0.687917`):
   ```text
   Relative Mean Perturbation = ( 0.0010280 / 0.687917 ) * 100% = 0.1494% of SD
   ```

2. **Worst-Case Individual Perturbation Bound (Max Sensitivity):**
   The maximum absolute perturbation observed across all 182 individuals is `0.0037480` (occurring in `BEB`):
   ```text
   Max Individual Perturbation = ( 0.0037480 / 0.687917 ) * 100% = 0.5448% of SD
   ```

3. **Standard Deviation Ratio of Target Loci vs. Genome-Wide Background:**
   Comparing the spread of the 6-variant sub-score (`SD_variants = 0.001507`) against the genome-wide score (`SD_Genova = 0.687917`):
   ```text
   Variant-to-Genome SD Ratio = ( 0.001507 / 0.687917 ) * 100% = 0.2191%
   ```

4. **Cross-Cohort Invariance of Target Stroke Loci:**
   Comparing the mean and standard deviation of the 6-SNP contribution between Pakistani (`PJL`) and Bangladeshi (`BEB`) cohorts:
   ```text
   Mean Difference (PJL vs BEB) = -0.000657 - (-0.000645) = -0.000012  (|Delta| = 1.2e-5)
   SD Difference   (PJL vs BEB) =  0.001517 -   0.001505  = +0.000012  (|Delta| = 1.2e-5)
   ```
   This near-zero divergence (`0.000012` in both first and second moments) proves that the 6 stroke loci exert highly uniform population-level genetic effects across distinct South Asian groups.

---

### 3. Non-Parametric Paired Hypothesis Testing (Wilcoxon Signed-Rank Test)

Because `Genova_full` and `Genova_ablated` represent paired observations on the same individuals-and because multi-SNP dosage sums over 6 discrete loci exhibit mild skewness (`Mean = -0.000651` vs. `Median = +0.000094`)-we performed a paired, non-parametric **Wilcoxon Signed-Rank Test** to evaluate the null hypothesis `H0: Median(Delta_Z_variants) = 0`.

For `n_nz` non-zero paired differences (excluding individuals whose net dosage contribution across the 6 SNPs equals exactly zero), let `W` denote the sum of signed ranks. Under `H0`, the asymptotic normal approximation has expected value `E[W]`, standard error `SE(W)`, and Z-statistic `Z_W`:

```text
E[W]  = n_nz * (n_nz + 1) / 4
SE(W) = SQRT[ n_nz * (n_nz + 1) * (2 * n_nz + 1) / 24 ]
Z_W   = ( W - E[W] ) / SE(W)
```

| Cohort | Total Sample (N) | Non-Zero Pairs (n_nz) | Zero-Delta Pairs | Observed Wilcoxon W | Expected E[W] | Std Error SE(W) | Standardized Z_W | Two-Sided p-value | Statistical Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **PJL** | 96 | 93 | 3 | 2014.0 | 2185.5 | 260.99 | -0.6571 | 0.5111 | **No systematic shift** (p > 0.05) |
| **BEB** | 86 | 83 | 3 | 1648.0 | 1743.0 | 220.26 | -0.4313 | 0.6662 | **No systematic shift** (p > 0.05) |
| **Combined SAS** | 182 | 176 | 6 | 7268.0 | 7788.0 | 676.90 | -0.7682 | 0.4424 | **No systematic shift** (p > 0.05) |

#### Why Non-Significance Here Is a Key Positive Validation:

A common misconception is that an ablation experiment should yield `p < 0.05`. In a **genome-wide polygenic model**, a statistically significant Wilcoxon test upon removing just 6 SNPs would indicate a severe pathology: namely, that 6 variants systematically bias the entire population's risk distribution in one direction.

* Instead, our results show that **176 of 182 individuals (96.70%)** carry non-zero risk perturbations from these 6 stroke loci (confirming active biological integration), while positive and negative risk allele combinations balance out across the population (`Z_W = -0.7682`, `p = 0.4424`).
* Consequently, Genova is **imputation-resilient**: if any of these candidate stroke SNPs fail genotyping quality control on a clinical microarray, the patient's overall polygenic risk percentile remains stable within **±0.55%** of a standard deviation.

---

### 4. Global & Local Concordance Analysis (`Genova` vs. `PRScs-European`)

To validate both our software implementation and the statistical behavior of Genova relative to an established external benchmark, we computed Pearson linear correlation coefficients (`r`) and coefficients of determination (`R² = r²`) at both the **genome-wide scale** and the **6-SNP variant scale**:

| Analytical Scale | Cohort | Sample Size (n) | Pearson Correlation (r) | Shared Linear Variance (R²) | Orthogonal / Calibrated Variance (1 - R²) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1. Genome-Wide PRS** (`Genova_full` vs. `PRScs_full`) | **PJL** | 96 | 0.8880 | 78.85% | 21.15% |
| **1. Genome-Wide PRS** (`Genova_full` vs. `PRScs_full`) | **BEB** | 86 | 0.8700 | 75.69% | 24.31% |
| **1. Genome-Wide PRS** (`Genova_full` vs. `PRScs_full`) | **Combined SAS** | 182 | 0.8750 | 76.56% | 23.44% |
| **2. Target 6-SNP Sub-Score** (`Delta_Z_Genova` vs. `Delta_Z_PRScs`) | **PJL** | 96 | 0.9970 | 99.40% | 0.60% |
| **2. Target 6-SNP Sub-Score** (`Delta_Z_Genova` vs. `Delta_Z_PRScs`) | **BEB** | 86 | 0.9970 | 99.40% | 0.60% |
| **2. Target 6-SNP Sub-Score** (`Delta_Z_Genova` vs. `Delta_Z_PRScs`) | **Combined SAS** | 182 | 0.99655 | 99.31% | 0.69% |

#### Engineering & Scientific Interpretation:

1. **99.31% Variant-Level Verification (r = 0.99655):**
   At the local 6-SNP level, Genova and `PRScs-European` achieve near-unity linear correlation across all 182 individuals. Because these two scores were computed through separate weight-processing pipelines, `r = 0.99655` provides an airtight **software verification proof**: our custom coordinate-mapping, strand-orientation, and `PLINK2` dosage-scoring scripts capture the exact biological directionality of the individual risk alleles without any sign flips or numerical corruption. Moreover, Genova scales the variant-level standard deviation down by **44.31%** (`0.001507` vs. `0.002706`), preventing uncalibrated European effect sizes from over-amplifying single-locus contributions.

2. **76.56% Shared Global Signal vs. 23.44% Ancestry Calibration (r = 0.8750):**
   At the genome-wide level, Genova preserves **76.56%** of the foundational polygenic variance present in `PRScs-European`, confirming that it remains anchored to established genetic architecture. Crucially, the remaining **23.44% orthogonal variance** (21.15% in `PJL`, 24.31% in `BEB`) represents the ancestry-aware corrections introduced by Genova's proprietary multi-ancestry architecture-which directly accounts for the **64.98% reduction in inter-cohort drift** and **56.26% compression in score variance** demonstrated in Section 1.

---

## 📂 Repository File Structure & Excel column description

The accompanying spreadsheets (`COMBINED_experiment_results.xlsx`, `PJL_experiment_results.xlsx`, `BEB_experiment_results.xlsx`) contain both the individual-level scoring matrices (n = 182) and the embedded statistical test summaries.

### Column Schema in Individual-Level Sheets:

* `IID`: 1000 Genomes Phase 3 sample identifier (e.g., `HG00096` for PJL, `HG03006` for BEB).
* `SAKSHAM-PRS`: Complete, reference-standardized genome-wide polygenic risk score (`SAKSHAM-PRS_full`).
* `SAKSHAM-PRS variants`: Standardized marginal contribution (`Delta_Z_variants`) of the 6 South Asian stroke-associated SNPs scaled against the full-genome reference standard deviation (`SD_ref_489`).
* `PRScs-European`: External baseline genome-wide PRS computed using standard European `PRS-CS` weights.
* `PRScs-European variants`: Marginal contribution of the same 6 variants under the European baseline framework.
* `Ablated SAKSHAM-PRS`: Counterfactual genome-wide score after removing the 6 target stroke variants (`= B2 - C2`, i.e., `Genova_full - Delta_Z_variants`).
* `Ablation Change`: Signed perturbation induced by variant removal (`= B2 - G2`, verifying exact numerical recovery of `Delta_Z_variants`).
* `Absolute Change`: Magnitude of the variant ablation effect (`= ABS(H2)`, i.e., `|Delta_Z_variants|`), used to bound individual-level sensitivity.

---

## Software & Computational Stack

* **Genomic Processing & Scoring:** `PLINK v2.00a` (`--score`, binary `.bed/.bim/.fam` manipulation), `VCFtools`
* **Data Engineering & Statistical Computing:** `Python 3.10+` (`pandas`, `numpy`, `scipy.stats`), `GNU Awk`, `Bash` automation
* **Reference & Validation Datasets:** 1000 Genomes Project Phase 3 (`GRCh37`/`GRCh38`) - South Asian Super-Population (`SAS` N = 489; validation cohorts `PJL` n = 96 and `BEB` n = 86)
