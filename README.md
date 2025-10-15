# Iron Homeostasis Model — MATLAB Implementation

This repository contains the MATLAB reproduction of the mechanistic iron homeostasis model developed in our lab.  
It reproduces the steady-state simulations under standard breastfeeding conditions described in the manuscript:

> *Modeling the influence of nutrition and supplementation on iron homeostasis during infancy.*

## 🔬 Contents

| Folder | Description |
|---------|--------------|
| `src/` | MATLAB scripts for ODEs, parameter import, and Hb estimation |
| `data/` | Parameter tables (Excel, CSV) |
| `results/` | Simulation outputs (steady-state, Hb traces) |
| `docs/` | Manuscript and supplementary materials |

---

## ⚙️ How to Run

1. Open MATLAB (R2022b+ recommended).  
2. Set working directory to this folder.  
3. Run:

   ```matlab
   run_standard_breastfeeding

