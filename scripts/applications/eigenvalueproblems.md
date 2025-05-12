# List with (Multiparameter) Eigenvalue Applications

The different applications are listed within one of the different application domains.
There are three different types of problems: system of square multiparameter eigenvalue problems (SMEP), single rectangular multiparameter eigenvalue problem (RMEP), and unified multiparameter eigenvalue problem (UMEP).

## System identification

- Least-squares identification of autoregressive moving-average (ARMA) models
    - type: RMEP
    - database: none
    - application: `arma.m` 
    - reference: vermeersch2019arma
- Least-squares realization of linear time-invariant autonomous systems
    - type: RMEP
    - database: none
    - application: `lsr.m`, `lsrwalsh.m`
    - reference: demoor2019least, demoor2020least

## Model order reduction 

- Globally optimal H2-norm model order reduction for single-input single-output systems
    - type: RMEP
    - database: `h2agudelo.csv`, `h2f*.csv`, `h2fom2r*.csv`, `h2fourdisk.csv`, `h2simple.csv`, `h2spanos3r4.csv`
    - application: `h2sisored*.m`
    - reference: agudelo2021globally, lagauw2023globally, hanzon2007optimal, hanzon1998model, alsubaie2019h2

## Systems theory

- Dynamic model updating
    - type: SMEP
    - database: none
    - application: none
    - reference: cottin2001dynamic
- Electromagnetism
    - type: SMEP
    - database: none
    - application: none
    - reference: morrison1976charge
- Analysis of aereolasticity, structural stability, and dynamics
    - type: SMEP
    - database: none
    - application: none
    - reference: pons2017multiparameter, pons2019nonlinear
- Hopf bifurcations
    - type: SMEP
    - database: none
    - application: none
    - reference: meerbergen2010shift
- Stability of delay-differential equations
    - type: SMEP
    - database: none
    - application: none
    - reference: jarlebring2009polynomial

## Cryptography

- Minimum-rank matrix problem
    - type: RMEP
    - database: none
    - application: none
    - reference: faugere2008cryptanalysis

## Differential equations

- Heine--Stieltjes problem
    - type: RMEP
    - database: none
    - application: none
    - reference: shapiro2009eigenvalues, shapiro2010algebro
- Seperable partial differential equations
    - type: SMEP
    - database: none
    - application: none
    - reference: gheorghiu2012spectral, plestenjak2015spectral, volkmer1988multiparameter, atkinson2011multiparameter, willatzen2011separable, willatzen2005numerical
- Wave and vibration simulation through structures
    - type: SMEP
    - database: none
    - application: none
    - reference: gravenkamp2025computation, kiefer2023computing
- Electromagnetic wave propagation
    - type: SMEP
    - database: none
    - application: none
    - reference: valovik2022multiparameter

## Multivariate polynomial systems and optimization problems

- Common roots of bivariate polynomial systems
    - type: SMEP
    - database: none
    - application: none 
    - reference: plestenjak2016roots, boralevi2017uniform, plestenjak2017minimal 
- Multivariate polynomial optimization of real-valued cost functions in complex variables
    - type: RMEP
    - database: none
    - application: none
    - reference: vermeersch2023multivariate, sorber2014numerical

