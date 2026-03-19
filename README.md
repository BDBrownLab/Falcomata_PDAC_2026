SERPINE1 and SERPINB2 orchestrate spatially dominant immunosuppressive niches in pancreatic cancer
===

> In this repository, we report how to reproduce the figures from our 2026 manuscript: Falcomatà, C., et al. (2026), SERPINE1 and SERPINB2 orchestrate spatially dominant immunosuppressive niches in pancreatic cancer. Additionally, we provide our data processing pipelines to pre-process single-cell, 10X Visium, and 10X Xenium datasets used in this study.

Table of Contents
===

### Abstract

Pancreatic ductal carcinoma (PDAC) is characterized by a highly immunosuppressive, ECM-rich microenvironment, yet tumors display striking heterogeneity. This raises the question of whether immune resistance is a global tumor property or is organized within spatially restricted niches. Using Perturb-map spatial functional genomics, we determine how different genes shape the growth and cellular environments of PDAC clones across space and time. This revealed early gene-driven remodeling of local immune neighborhoods preceding late-stage spatial clonal dominance. We identify SERPINE1 and SERPINB2 as dominant regulators of tumor microenvironment (TME) control and immune evasion. These tumor-derived fibrinolysis regulators promote the stabilization of fibrin-rich ECM niches that spatially retain and program macrophages toward immunosuppressive states while excluding cytotoxic T cells. Consistently, loss of Serpine1 or Serpinb2, or pharmacologic inhibition of PAI-1, improves tumor control and synergizes with PD-1 blockade to prolong survival. Multimodal spatial analysis of patient tumors revealed that immunosuppressive niches form around rare SERPINB2- and SERPINE1-expressing PDAC cell subpopulations, dominated by proliferating and SPP1+MARCO+ macrophages. Together, our findings identify cancer cell-derived fibrinolysis regulators as local spatial organizers of immune suppression, linking tumor-intrinsic heterogeneity to local microenvironmental control and therapeutic resistance in PDAC. 

### Setup and Environments

| Environment Purpose        | Pertaining Figures           |
| ------------- |:-------------:|
|  Perturb-map analysis | 1C, 1D, S4A-C, 5K, 5L |
|  Perturb-map Multi-modal analysis | 3A, 3B, S5D, S5G, S5H |
|  ScRNAseq analysis | 2F, 5A-C, S5G, S10 |
|  Visium analysis  |  4P, 4Q, S9  |
|  BulkRNAseq analysis |    2A, 2B, S5A, S5I   |

### Preprocessing

1. [Visium download and pre-processing pipeline]()

