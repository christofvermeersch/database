function [system] = interval()
    eqs = cell(10,1);
    eqs{1} = [1 1 0 0 0 0 0 0 0 0 0; -0.18324757 0 0 0 1 1 0 0 0 0 1; -0.25428722 0 0 0 0 0 0 0 0 0 0];
    eqs{2} = [-0.16275449 1 1 0 0 0 0 1 0 0 0; 1 0 0 1 0 0 0 0 0 0 0; -0.37842197 0 0 0 0 0 0 0 0 0 0];
    eqs{3} = [-0.16955071 1 1 1 0 0 0 0 0 0 0; 1 0 0 0 1 0 0 0 0 0 0; -0.27162577 0 0 0 0 0 0 0 0 0 0];
    eqs{4} = [-0.15585316 1 0 0 0 0 0 1 1 0 0; 1 0 0 0 0 1 0 0 0 0 0; -0.19807914 0 0 0 0 0 0 0 0 0 0];
    eqs{5} = [-0.1995092 0 0 0 1 0 0 1 1 0 0; 1 0 0 0 0 0 1 0 0 0 0; -0.44166728 0 0 0 0 0 0 0 0 0 0];
    eqs{6} = [-0.18922793 0 1 0 0 0 1 0 0 1 0; 1 0 0 0 0 0 0 1 0 0 0; -0.14654113 0 0 0 0 0 0 0 0 0 0];
    eqs{7} = [-0.21180484 0 0 1 0 0 1 0 0 1 0; 1 0 0 0 0 0 0 0 1 0 0; -0.42937161 0 0 0 0 0 0 0 0 0 0];
    eqs{8} = [-0.17081208 1 0 0 0 0 0 1 1 0 0; 1 0 0 0 0 0 0 0 0 1 0; -0.07056438 0 0 0 0 0 0 0 0 0 0];
    eqs{9} = [-0.1961274 0 1 0 0 0 0 1 0 1 0; 1 0 0 0 0 0 0 0 0 0 1; -0.34504906 0 0 0 0 0 0 0 0 0 0];
    eqs{10} = [-0.21466544 1 0 0 0 1 0 0 0 1 0; 1 0 1 0 0 0 0 0 0 0 0; -0.42651102 0 0 0 0 0 0 0 0 0 0];

    system = systemstruct(eqs);
end