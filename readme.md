# Database of Test Problems

## File structure
|
|-> Functions
|	|-> files to manipulate the database and use the problems
|-> Problems
|	|-> Applications
|	|	|-> files to solve applications
|	|-> EigenvalueProblems
|	|	|-> eigenvalue problems
|	|-> PolynomialSystems
|		|-> systems of multivariate polynomial equations
|-> Scripts
| 	|-> scripts to recreate results in papers
|-> Solutions
| 	|-> MacaulayLab
|		|-> solutions of problems solved via MacaulayLab
|-> database.csv
|-> readme.md

## Explanation of the different files
- addproblem.m: personal function to add a problem from MATLAB (note that it uses a hardcoded path -> excluded from final database)
- loaddatabase.m: load and display the current database
- loadproblem.m: load a problem from the database + copy all information from database.csv into the problemstruct
- readproblem.m: read problem.csv without consulting database.csv (will be much faster than loadproblem  when the database is large)
- writeproblem.m: write problem.csv without consultin database.csv

- arma11.m: function to create MEP that solves ARMA identification problem

- arma11l4.csv: MEP (constructed via arma11.m) that solves a particular identification problem
- dreesen1.csv: system
- dreesen4.csv: system
- toy1.csv: system

- 23-47.m: script that corresponds to Sibren's paper

- dreesen1.sol: solutions obtained via MacaulayLab when solving dreesen1.csv
