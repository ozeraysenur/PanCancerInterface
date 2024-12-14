# 🧬 **Pan-Cancer Analysis User Interface**

The **Pan-Cancer Analysis User Interface** is an advanced interactive tool designed to facilitate the analysis of genomic variations, mutation patterns, and multi-omics data across various cancer types. This platform empowers researchers, bioinformaticians, and clinicians with powerful visualizations and data insights to uncover meaningful patterns in cancer-related data.

---

## 🌟 **Key Features**

### 🔬 **Multi-Omics Analysis**
- **Copy Number Variations (CNV)**: Visualize large-scale changes in genomic structure with customizable heatmaps.
- **Transcriptome Data**: Explore gene expression profiles and identify potential biomarkers.
- **Single Nucleotide Variations (SNV)**: Pinpoint mutation hotspots across different cancer types.

### 📊 **Interactive Visualizations**
- **Dynamic Heatmaps**: Customize color schemes, annotations, and clustering options.
- **Correlation Matrices**: Analyze relationships between genomic features and expression levels.
- **Pie Charts & Bar Graphs**: Gain deeper insights into mutation distributions and gene expression profiles.

### 🏰 **Global Overview**
- Integration of geographic data for global comparisons in mutation prevalence.

### ⚙️ **Advanced Customization**
- Fully customizable visualization parameters to suit diverse analytical needs.
- Real-time data filtering for targeted analysis.

---

## 💻 **Technology Stack**

- **Frontend**: Shiny framework for a responsive and user-friendly interface.
- **Backend**: R for data processing, machine learning, and visualization.
- **Visualization Tools**: ggplot2, Plotly for interactive and publication-quality graphics.
- **Data Handling**: data.table for efficient manipulation of large datasets.

---

## 🚀 **Getting Started**

### **Prerequisites**
1. **Install R** (version 4.0 or newer).
2. Install the required R packages:
   ```R
   install.packages(c("shiny", "ggplot2", "data.table", "plotly", "DT", "reshape2"))
   ```

### **Clone the Repository**
   ```bash
   git clone https://github.com/ozeraysenur/SeniorProject.git
   cd SeniorProject
   ```

### **Run the Application**
   ```R
   shiny::runApp("app")
   ```

### **Access the App**
Navigate to `http://localhost:8100` in your web browser.

---

## 📂 **Repository Structure**

```
SeniorProject/
├── app/                   # Shiny app components
│   ├── ui.R               # User interface
│   ├── server.R           # Backend server logic
│   └── helpers.R          # Helper functions
├── data/                  # Input datasets
├── assets/                # Images and custom CSS
├── scripts/               # Data preprocessing and analysis scripts
└── README.md              # Project documentation
```

---

## 📊 **Dataset Description**

- **Source**: Public cancer genome datasets.
- **Structure**:
  - Mutation profiles per sample.
  - Gene expression levels across cancer types.
  - Clinical annotations for patient stratification.

---

## ✨ **Advanced Features**

### 🔍 **Data Insights**
- Heatmap visualizations with hierarchical clustering.
- Gene and sample filtering for focused analysis.
- Mutation density distributions.

### 📦 **Data Export**
- Export filtered and analyzed data in CSV and JSON formats.

### 🧠 **Machine Learning Integration**
- Incorporation of clustering algorithms to identify sample subtypes.
- Regression analysis for gene expression predictions.

---

## 📜 **How to Contribute**

We welcome contributions!  
1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Submit a pull request with a clear description of your changes.

---

## 💜 **Acknowledgements**

This project was developed as part of a senior thesis project to advance cancer research through data visualization and analysis. A heartfelt thank you to all mentors, collaborators, and the open-source community for their support.

Made with 💜 by **Ayşenur Özer**
