# Multi-Scale Chromatin Polymer Physics & Scaling Analysis

This repository contains a reproducible computational pipeline to analyze the spatial contact probability scaling laws of the human genome. Using Hi-C chromosomal conformation capture data from the ENCODE project, we investigate the biophysical folding states of chromatin by calculating the contact scaling exponent ($\alpha$) across varying genomic distances to test the **Fractal Globule** vs. **Equilibrium Globule** polymer models.

---

## 🧬 Scientific Background

The spatial architecture of the genome is fundamentally linked to gene regulation. In their seminal work, *Lieberman-Aiden et al. (2009)* demonstrated that chromatin packaging obeys physical scaling laws:

$$\text{Probability of Contact } P(s) \propto s^\alpha$$

Where $s$ represents the genomic distance along the chromosome, and $\alpha$ is the scaling exponent:
*   **$\alpha \approx -1.0$ (Fractal Globule):** An unknotted, highly accessible, self-similar folding state optimized for dynamic transcriptional activation.
*   **$\alpha \approx -1.5$ (Equilibrium Globule):** A maximally tangled, highly dense, collapsed polymer state driven by thermodynamic entropy.

---

## 📊 Summary of Key Findings

This pipeline analyzes the sub-telomeric region of Chromosome 1 (`chr1:0-10,000,000`) in human K562 myeloid leukemia cells at a 500kb resolution.

| Analysis State | Scaling Exponent ($\alpha$) | Biophysical Interpretation |
| :--- | :--- | :--- |
| **Raw Uncorrected Data** | $-1.45$ | Dominated by sequencing and mappability bias. |
| **KR-Balanced (Normalized)** | $-1.56$ | Matches collapsed polymer physics at large scales. |
| **Short-Range (500kb - 2.5Mb)** | **$-1.59$** | Classically inactive sub-telomeric heterochromatin. |
| **Long-Range (3.0Mb - 10.0Mb)** | **$-1.65$** | Random-walk entropy scaling at macroscopic distances. |

### Biophysical Discussion
Our observation of $\alpha \approx -1.6$ at the sub-telomeric tip of Chromosome 1 represents the physical signature of **constitutive heterochromatin**. Lacking active loop-extrusion machinery (e.g., cohesin loops), this gene-poor, transcriptionally silent region physically collapses into an entropy-dominated **equilibrium globule**. Furthermore, our multi-scale analysis demonstrates the transition toward steeper, random-walk polymer states ($\alpha \to -1.65$) at macroscopic distances (>3Mb).

---

## 🛠️ Pipeline Architecture & Dependencies

This pipeline is optimized for cloud Linux platforms (like Google Colab) and relies on robust, industry-standard biophysical file formats and APIs:
*   **`hic2cool`**: Binary format parser to convert multi-resolution `.hic` archives to `.cool` HDF5 containers.
*   **`cooler`**: Python API to fetch and query sparse contact matrices.
*   **`matplotlib` & `numpy`**: Array processing, data filtering, and vector fitting.

### Installation
```bash
pip install cooler hic2cool numpy pandas matplotlib
