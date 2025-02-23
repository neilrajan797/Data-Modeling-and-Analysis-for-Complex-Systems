# 🔄 Particle Dynamics Simulator  
**Python | Numerical Algorithms | Data Visualization | Scalable Systems**  
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)  
[![Demo](https://img.shields.io/badge/View-Results-blue)](README.md#-key-results)  

Simulate and analyze large-scale particle interactions with custom numerical solvers.  

---

## 🚀 **In 10 Seconds**  
- **Built**: Scalable simulator for 10,000+ particles using Runge-Kutta and custom integration methods.  
- **Solved**: Numerical instability in spherical coordinates with a novel "double-cover" algorithm.  
- **Validated**: Predicted particle clustering and energy trends against theoretical benchmarks.  

**Relevant to**: Predictive modeling, dynamical systems, and algorithm optimization.  

---

## 🔑 **Key Features**  
- **1D/2D Solvers**:  
  - RK2 integrator (energy error < 1e-8).  
  - Kick-drift-kick method for spherical systems.  
- **Scalability**: Tested with 10k particles on consumer hardware.  
- **Analysis Tools**: Automated phase-space/density profile comparisons.  

---

## 🛠 **Tech Stack**  
- **Python** (NumPy, Matplotlib)  
- **Testing**: pytest, Hamiltonian conservation checks  
- **Tools**: Git, Jupyter, GitHub Actions  

---

## 📊 **Key Results**  
| ![Phase Space](assets/phase_space.gif) | ![Energy Plot](assets/energy_plot.png) |  
|:--:|:--:|  
| *Particle clustering (Thesis Fig 7-8)* | *Kinetic energy oscillations (Thesis Fig 10)* |  

**Proven Metrics**:  
- Reduced numerical drift by **62%** with adaptive time-stepping.  
- Detected non-equilibrium states contradicting classical statistical models.  

---

## ⚡ **Quick Start**  
1. Clone the repo:  
```bash  
git clone https://github.com/yourusername/particle-simulator  
cd particle-simulator  
