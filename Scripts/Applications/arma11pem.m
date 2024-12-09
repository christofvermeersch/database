function [system] = armapem(y) 
    % TODO: extend code to work with any number of data points.
    eqs = cell(9,1); 
    eqs{1} = [1 zeros(1,6) 1 0 0; 1 0 1 zeros(1,3) 1 zeros(1,3); -y(2) zeros(1,9); -y(1) 1 zeros(1,8)]; 
    eqs{2} = [1 zeros(1,7) 1 0; 1 0 1 zeros(1,4) 1 zeros(1,2); -y(3) zeros(1,9); -y(2) 1 zeros(1,8)];
    eqs{3} = [1 zeros(1,8) 1; 1 0 1 zeros(1,5) 1 zeros(1,1); -y(4) zeros(1,9); -y(3) 1 zeros(1,8)];
    eqs{4} = [y(1) 0 0 1 0 0 zeros(1,4); y(2) 0 0 0 1 0 zeros(1,4); y(3) 0 0 0 0 1 zeros(1,4)];
    eqs{5} = [ones(3,1) zeros(3,2) eye(3) eye(3) zeros(3,1)];
    eqs{6} = [1/4 zeros(1,5) 1 0 0 0; 1 0 1 1 0 0 zeros(1,4)];
    eqs{7} = [1/4 zeros(1,5) 0 1 0 0; 1 0 1 0 1 0 zeros(1,4); 1 0 0 1 0 0 zeros(1,4)];
    eqs{8} = [1/4 zeros(1,5) 0 0 1 0; 1 0 1 0 0 1 zeros(1,4); 1 0 0 0 1 0 zeros(1,4)];
    eqs{9} = [1/4 zeros(1,5) 0 0 0 1; 1 0 0 0 0 1 zeros(1,4)];

    system = systemstruct(eqs); 
end