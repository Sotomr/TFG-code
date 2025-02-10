
# Multi-Objective Optimization in Distributed Systems Using Game Theory  

This repository contains the implementation of a **multi-objective optimization model** for coordinating distributed autonomous systems using **Nash equilibrium**. The project focuses on task allocation in drone networks, leveraging **game theory** for decentralized decision-making. The implementation is developed in MATLAB and includes various modifications and extensions of the base model.  

## Code Structure  

- **`EQMO.m`** – Base implementation of the optimization model.  
- **`EQMOv1.m`, `EQMOv2.m`, `EQMOv3.m`, `EQMOy.m`** – Modified and extended versions of the base model.  
- **`gEQMO.m`** – Executes and visualizes the results of the base model.  
- **`ComFINAL.m`** – Executes and compares all versions to assess their performance under different conditions.  
- **`ConvergenciaEQMO.m`** – Analyzes the convergence properties of all models, providing insights into their computational behavior and stability.  

## Features  

- **MATLAB implementation** of optimization algorithms for multi-agent coordination.  
- **Computational analysis** of model convergence and efficiency.  
- **Graphical visualization** of different model versions and their comparative performance.  
- **Mathematical formulation & documentation** of the model and its theoretical foundations.  

## Usage  

1. Run **`EQMO.m`** to execute the base model.  
2. Use **`gEQMO.m`** to generate graphical representations of the base model.  
3. Test alternative versions by running **`EQMOv1.m`**, **`EQMOv2.m`**, **`EQMOv3.m`**, or **`EQMOy.m`**.  
4. Compare all models using **`ComFINAL.m`** to assess their relative performance.  
5. Analyze convergence trends by executing **`ConvergenciaEQMO.m`**.  
