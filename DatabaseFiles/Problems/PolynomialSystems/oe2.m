function [system] = oe2(u,y) 
    eqs = cell(7,1); 
    eqs{1} = [y(3) 0 0 0 0 0 0 0; y(2) 0 1 0 0 0 0 0; y(1) 0 0 1 0 0 0 0; ...
        -u(2) 1 0 0 0 0 0 0; -6 0 0 0 1 0 0 0; -6 0 1 0 0 2 0 0; -6 0 0 1 0 0 1 0; ...
        -6 0 2 0 1 0 0 0; -6 0 1 1 0 1 0 0; -6 0 0 2 1 0 0 0]; 
    eqs{2} = [y(4) 0 0 0 0 0 0 0; y(3) 0 1 0 0 0 0 0; y(2) 0 0 1 0 0 0 0; ...
        -u(3) 1 0 0 0 0 0 0; -6 0 0 0 0 1 0 0; -6 0 1 0 0 2 0 0; -6 0 0 1 0 0 0 1; ...
        -6 0 1 0 1 0 0 0; -6 0 2 0 0 1 0 0; -6 0 1 1 0 0 1 0; -6 0 1 1 1 0 0 0; ...
        -6 0 0 2 0 1 0 0];
    eqs{3} = [y(5) 0 0 0 0 0 0 0; y(4) 0 1 0 0 0 0 0; y(3) 0 0 1 0 0 0 0; ...
        -u(4) 1 0 0 0 0 0 0; -6 0 0 0 0 0 1 0; -6 0 1 0 0 0 0 1; -6 0 1 0 0 1 0 0; ...
        -6 0 2 0 0 0 1 0; -6 0 1 1 0 0 0 1; -6 0 0 1 1 0 0 0; -6 0 1 1 0 1 0 0; ...
        -6 0 0 2 0 0 1 0];
    eqs{4} = [y(6) 0 0 0 0 0 0 0; y(5) 0 1 0 0 0 0 0; y(4) 0 0 1 0 0 0 0; ...
        -u(5) 1 0 0 0 0 0 0; -6 0 0 0 0 0 0 1; -6 0 1 0 0 0 1 0; -6 0 2 0 0 0 0 1; ...
        -6 0 0 1 0 1 0 0; -6 0 1 1 0 0 1 0; -6 0 0 2 0 0 0 1];
    eqs{5} = [y(2) 0 0 0 1 0 0 0; -6 0 1 0 2 0 0 0; -6 0 0 1 1 1 0 0; y(3) 0 0 0 0 1 0 0; ...
        -6 0 0 0 1 1 0 0; -6 0 1 0 0 2 0 0; -6 0 0 1 0 1 1 0; y(4) 0 0 0 0 0 1 0; ...
        -6 0 0 0 0 1 1 0; -6 0 1 0 0 0 2 0; -6 0 0 1 0 0 1 1; y(5) 0 0 0 0 0 0 1; ...
        -6 0 0 0 0 0 1 1; -6 0 1 0 0 0 0 2];
    eqs{6} = [y(1) 0 0 0 1 0 0 0; -6 0 0 1 2 0 0 0; y(2) 0 0 0 0 1 0 0; -6 0 1 0 1 1 0 0; -6 0 0 1 0 2 0 0; ...
        y(3) 0 0 0 0 0 1 0; -6 0 0 0 1 0 1 0; -6 0 1 0 0 1 1 0; -6 0 0 1 0 0 2 0; ...
        y(4) 0 0 0 0 0 0 1; -6 0 0 0 0 1 0 1; -6 0 1 0 0 0 1 1; -6 0 0 1 0 0 0 2];
    eqs{7} = [-u(2) 0 0 0 1 0 0 0; -u(3) 0 0 0 0 1 0 0; -u(4) 0 0 0 0 0 1 0; ...
        -u(5) 0 0 0 0 0 0 1];

    system = systemstruct(eqs); 
end