# 🧬 Breast Cancer Gene Expression ( Microarray data) analysis and Machine Learning modules practice outcomes

**Author**: Gayatri Sunil Samant  
**Contact**: 8431036658  
**Email**: gayatrisamant05@gmail.com  
**Project**: practice project for NyberMan 15-day internship  
**Duration**: 19th July to 3rd August 2025

This project involves the analysis of breast cancer microarray data (GSE21422) to identify **differentially expressed genes (DEGs)** and apply **machine learning models** for predictive insights. The data was analyzed using R/Bioconductor and includes preprocessing, normalization, visualization, DEG identification, and classification modeling.

---

## 📁 Dataset

- **Source**: [GEO - GSE21422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21422)
- **Platform**: Affymetrix Human Genome U133 Plus 2.0 Array
- **Samples**: Breast cancer and normal tissue samples in `.CEL` format

---

## 📌 Objectives

- Import and normalize microarray data
- Perform quality control and exploratory analysis
- Identify differentially expressed genes (DEGs)
- Visualize DEGs using Volcano plot, PCA, heatmaps
- Apply machine learning algorithms to classify cancer vs normal

---

## 🔧 Tools and Libraries

- **R Programming Language**
- **Bioconductor Packages**:
  - `affy`, `limma`, `oligo`, `GEOquery`, `pheatmap`, `ggplot2`
- **Machine Learning**:
  - `caret`, `randomForest`, `e1071`

---

## 🧪 Methodology

1. **Data Import and Preprocessing**
   - Load `.CEL` files using `affy` or `oligo`
   - Background correction and normalization (RMA)

2. **Quality Control**
   - Boxplots
   - PCA (Principal Component Analysis)

3. **Differential Expression Analysis**
   - Design matrix creation
   - Linear modeling using `limma`
   - Volcano plot for DEGs

4. **Visualization**
   - PCA plot
   - Heatmap of top DEGs
   - Volcano plot

5. **Machine Learning Modeling**
   - Feature selection (top DEGs)
   - Model training (Random Forest, SVM)
   - Accuracy evaluation

---

## 📊 Output

- Normalized expression matrix
- List of DEGs with logFC and adjusted p-values
- PCA and heatmap visualizations
- ML classification performance metrics

---

## 📎 Folder Structure

Breast_Cancer_data_analysis_using_ML<br>
├── CEL_files/ # Raw microarray data<br>
├── Scripts/ # R scripts for each analysis step<br>
├── Plots/ # PCA, heatmaps, volcano plots<br>
├── Results/ # DEG tables and ML results<br>
└── README.md

---

## 📚 Internship Project

This analysis was completed as part of the **NyBerMan Bioinformatics Internship (July 2025)** under the theme _"AI-Driven Genomic Data Analysis"_.


---

## 🔗 References

- GEO Accession: [GSE21422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21422)
- Bioconductor: https://www.bioconductor.org/
- LIMMA User Guide: https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

---

## 🧠 Keywords

`Breast Cancer`, `Microarray`, `Machine Learning`, `DEG`, `R`, `GSE21422`, `limma`, `affy`, `SVM`, `R
