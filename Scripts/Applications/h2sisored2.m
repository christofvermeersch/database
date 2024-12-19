function [mat,supp] = h2sisored2(an,bn,r)
    %H2SISORED2   H2-norm model order reduction problem.
    %   mat = H2SISORED2(an,bn,r) creates the coefficient matrices of the
    %   multiparameter eigenvalue problem that solves the H2-norm model 
    %   order reduction problem of an n-th order model to an r-th order 
    %   model. It uses the so-called "second approach" from [1]. The given
    %   n-th order model is given by its transfer function b(1)*s^{n-1} + 
    %   b(1)*s^{n-2} + ... + b(n) / s^n + a(1)*s^{n-1} + ... + a(n).
    %
    %   [mat,supp] = H2SISORED2(...) also gives the support of the 
    %   coefficient matrices.
        
    % Copyright (c) 2024 - Christof Vermeersch
    %
    % This function implements the construction of the numerical
    % coefficient matrices as described in [1].
    %
    % Note: this approach results in the interpolation points of the 
    % reduced order model. You find more information about them in [1].
    %
    % [1] Alsubaie, F. F. F. (2019), H2 Optimal Model Reduction for Linear
    % Dynamic Systems and the Solution of Multiparameter Matrix Pencil 
    % Problems, PhD thesis, Imperial College London, London, UK.

    % Determine the degree:
    n = length(an);

    % Create the auxiliary matrices:
    A = [-an.'; eye(n-1) zeros(n-1,1)];
    B = zeros(n,1); B(1) = 1;
    c = bn.';
    C = zeros(r,n);
    for i = 1:r
        C(i,:) = c*A^(i-1);
    end
    

    % Derive the (numerical) coefficient matrices:
    mat = cell(r+1,1);
    mat{1} = [zeros(n,n) (-A)^r A^r*B; (-A)^r -eye(n) zeros(n,1); ...
        C zeros(r,n+1)];
    for i = 1:r
        temp = [zeros(n,n) (-A)^(r-i) A^(r-i)*B; (-A)^(r-i) ...
            zeros(n,n+1); zeros(r,2*n+1)];
        mat{i+1} = temp;
    end
    supp = [zeros(1,r); eye(r)];
end