# Database of Test Problems

## Folder structure
```
|-> Functions
|	|-> files to manipulate the database and use the problems
|-> Problems
|	|-> EigenvalueProblems
|	|	|-> eigenvalue problems
|	|-> PolynomialSystems
|		|-> systems of multivariate polynomial equations
|-> Scripts
| 	|-> Papers
|   |   |-> scripts to recreate results in papers
|   |-> Applications
|       |-> scripts to generate problems (e.g., arma.m or lsrwalsh.m)
|-> Solutions
|	|-> solutions of problems solved via MacaulayLab
|-> database.yaml
|-> bibids.bib
|-> readme.md
```
## Explanation of the different files
- readproblem.m: read problem.csv (there is also an option to load information from the database)
- writeproblem.m: write problem.csv without consulting database.yaml

- problem.csv: file that contains the entries and support of the problem

- arma.m: script that creates an MEP to identify an autoregressive moving-average model
- lagauw2023h2.m: script that corresponds to Sibren's paper

- dreesen1sol.csv: solutions obtained via MacaulayLab when solving dreesen1.csv
