function [system] = chandra6()
    eqs = cell(6,1);
    eqs{1} = [-0.25617 2 0 0 0 0 0; -0.17078 1 1 0 0 0 0; -0.128085 1 0 1 0 0 0; -0.102468 1 0 0 1 0 0; -0.08539 1 0 0 0 1 0; 11.48766 1 0 0 0 0 0; -12 0 0 0 0 0 0];
    eqs{2} = [-0.34156 1 1 0 0 0 0; -0.25617 0 2 0 0 0 0; -0.204936 0 1 1 0 0 0; -0.17078 0 1 0 1 0 0; -0.146382857142857 0 1 0 0 1 0; 11.48766 0 1 0 0 0 0; -12 0 0 0 0 0 0];
    eqs{3} = [-0.384255 1 0 1 0 0 0; -0.307404 0 1 1 0 0 0; -0.25617 0 0 2 0 0 0; -0.219574285714286 0 0 1 1 0 0; -0.1921275 0 0 1 0 1 0; 11.48766 0 0 1 0 0 0; -12 0 0 0 0 0 0];
    eqs{4} = [-0.409872 1 0 0 1 0 0; -0.34156 0 1 0 1 0 0; -0.292765714285714 0 0 1 1 0 0; -0.25617 0 0 0 2 0 0; -0.227706666666667 0 0 0 1 1 0; 11.48766 0 0 0 1 0 0; -12 0 0 0 0 0 0];
    eqs{5} = [-0.42695 1 0 0 0 1 0; -0.365957142857143 0 1 0 0 1 0; -0.3202125 0 0 1 0 1 0; -0.284633333333333 0 0 0 1 1 0; -0.25617 0 0 0 0 2 0; 11.48766 0 0 0 0 1 0; -12 0 0 0 0 0 0];
    eqs{6} = [-0.439148571428571 1 0 0 0 0 1; -0.384255 0 1 0 0 0 1; -0.34156 0 0 1 0 0 1; -0.307404 0 0 0 1 0 1; -0.279458181818182 0 0 0 0 1 1; 11.48766 0 0 0 0 0 1; -12 0 0 0 0 0 0];

    system = systemstruct(eqs);
end