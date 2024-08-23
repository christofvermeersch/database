function [system] = dreesen7() 
    % DREESEN7 contains a system of multivariate polynomial equations.
    % 
    %   [system] = DREESEN7() returns the system of multivariate polynomial 
    %   equations.
    %
    %   Note: Because of the multiplicity of the affine solutions, the 
    %   optional clustering step is essential while solving this system of 
    %   multivariate polynomial equations.

    % MacaulayLab (2022) - Christof Vermeersch. 

    eqs = cell(2,1); 
    eqs{1} = [1 0 2; -4 0 1; 4 0 0]; 
    eqs{2} = [1 2 0; -2 1 1; 2 1 0; 1 0 2; -2 0 1; 1 0 0]; 

    system = systemstruct(eqs); 
end