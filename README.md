# prioritizer-pmatch
Prioritize personalized treatment options for cancer patients based on integrated clinical and genomic profiles.

**`prioritizer-pmatch`** is an R-based framework designed to identify high-priority therapeutic options for patients using multi-omic signatures and clinical biomarkers. It supports ranking algorithms, signature-based scoring, and cohort stratification for precision oncology applications.

---

## Features

- Computes gene signature scores using methods like GSVA, ssGSEA, Weighted Mean, and specific methods.
- Loads and processes gene expression matrices from raw CSV.
- Scales signature scores for visualization and comparative analysis.
- Generates boxplots with fixed threshold cutoffs to highlight selected patients.
- Supports integration with curated signature metadata.

---

## Repository Structure

```
RibaA/
├── data/
│   ├── rawdata/                  # Raw expression and signature files
|   ├── procata/                  # Selected IO signatures
│   ├── results/                  # Output directory for scores and plots
├── scripts/
│   ├── runSigAnalysis.R          # Main script for signature scoring
│   ├── runPatientSelection.R     # Select patients
│   ├── runVisualization.R        # Script for creating boxplots
├── README.md
```

---