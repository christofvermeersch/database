function [mat,supp] = h2sisored3(an,bn,r)
    %H2SISORED3   H2-norm model order reduction problem.
    %   mat = H2SISORED3(an,bn,r) creates the coefficient matrices of the
    %   multiparameter eigenvalue problem that solves the H2-norm model 
    %   order reduction problem of an n-th order model given by its 
    %   coefficients of the transfer function (bn/an) to an r-th order 
    %   model. It uses the so-called "third (Walsh) approach" from [1].
    %
    %   [mat,supp] = H2SISORED3(...) also gives the support of the 
    %   coefficient matrices.
        
    % Copyright (c) 2024 - Christof Vermeersch
    %
    % This function implements the construction of the numerical
    % coefficient matrices as described in [1].
    %
    % [1] Lagauw, S., Agudelo, O. M., and De Moor, B. (2023) Globally
    % Optimal SISO H2-norm Model Reduction Using Walsh's Theorem, IEEE
    % Control System Letters, vol. 7, p. 1970-1675.

    error("Problem is not yet implemented.")
end