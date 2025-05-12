# Database of Test Problems

## Project description

## Folder structure

```
|-> functions
|	|-> files to manipulate the database and use the problems
|-> problems
|	|-> eigenvalueproblems
|	|	|-> eigenvalue problems
|	|-> polynomialsystems
|		|-> systems of multivariate polynomial equations
|-> scripts
| 	|-> papers
|   |   |-> scripts to recreate results in papers
|   |-> applications
|       |-> scripts to generate problems (e.g., arma.m or lsrwalsh.m)
|       |-> a list of multiparameter eigenvalue applications
|-> Solutions
|	|-> solutions of problems solved via MacaulayLab
|-> database.yaml
|-> bibids.bib
|-> readme.md
```

- readproblem.m: read problem.csv (there is also an option to load information from the database)
- writeproblem.m: write problem.csv without consulting database.yaml
- createreport.m: export the problem in LaTeX format
- linearize.m: linearize multiparameter eigenvalue problems

- problem.csv: file that contains the entries and support of the problem

- arma.m: script that creates an MEP to identify an autoregressive moving-average model
- lagauw2023h2.m: script that corresponds to Sibren's paper

- dreesen1.sol: solutions obtained via MacaulayLab when solving dreesen1.csv

## Related reference(s)
