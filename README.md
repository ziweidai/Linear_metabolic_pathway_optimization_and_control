# Linear_metabolic_pathway_optimization_and_control

This repository contains MATLAB codes for simulating enzyme allocation, flux efficiency, and flux-control trade-offs in linear metabolic pathways under different rate laws.

---

## First_order

Codes simulating flux-control trade-off under first-order kinetics.

- compare_flux_efficiency_control_efficiency.m
  Computes Pareto front balancing flux efficiency and control efficiency for glycolysis and TCA cycle under first-order kinetics. Parameter values are based on literature.

---

## Michaelis_Menten

Codes simulating maximal flux efficiency solution under Michaelis-Menten kinetics.

- RandSample_Linear.m
  Simulates enzyme allocation profiles maximizing pathway flux under Michaelis-Menten kinetics. Parameters are sampled from literature-based values.
