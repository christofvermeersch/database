function [mat,supp] = lsr(y,n)
    %LSR   Least-squares realization problem.
    %   mat = LSR(y,n) creates the coefficient matrices of the
    %   multiparameter eigenvalue problem that solves the n-th order least-
    %   squares realization problem for the output data in the vector y.
    %
    %   [mat,supp] = LSR(...) also gives the support of the coefficient
    %   matrices.
        
    % Copyright (c) 2024 - Christof Vermeersch
    %
    % This function implements the construction of the numerical
    % coefficient matrices as described in [1,2].
    %
    % [1] De Moor, B. (2019) Least Squares Realization of LTI Models Is an
    % Eigenvalue Problem, Proc. of the 18th European Control Conference
    % (ECC), Naples, Italy.
    % [2] De Moor, B. (2020) Least Squares Optimal Realization of
    % Autonomous LTI Systems Is an Eigenvalue Problem, Communications in
    % Information Systems, vol. 20, no. 2, p. 163-207.

    % Define the (symbolic) variables:
    N = length(y);
    a = sym("a",[n 1]);
    for i = 1:n
        assume(a(i),"real")
    end

    % Create the auxiliary matrices:
    t = [fliplr(a') 1 zeros(1,N-n-1)];
    Ta = sym('Ta',[N-n,N]);
    for i = 1:N-n
        Ta(i,:) = circshift(t,i-1);
    end
    Da = Ta*Ta';
    
    % Construct matrix polynomial:
    A = sym("A",[(N-n)*(n+1)+n (N-n)*(n+1)+1]);
    A(1:N-n,1:N-n) = Da;
    A(1:N-n,N-n+1:(N-n)*(n+1)) = zeros(N-n,(N-n)*n);
    A(1:N-n,end) = Ta*y;
    for i = 1:n
        A(i*(N-n)+1:(i+1)*(N-n),:) = zeros(N-n,(N-n)*(n+1)+1);
        A(i*(N-n)+1:(i+1)*(N-n),1:N-n) = diff(Da,a(i));
        A(i*(N-n)+1:(i+1)*(N-n),i*(N-n)+1:(i+1)*(N-n)) = Da;
        A(i*(N-n)+1:(i+1)*(N-n),end) = diff(Ta,a(i))*y;
    end 
    for i = 1:n
        A((N-n)*(n+1)+i,:) = zeros(1,(N-n)*(n+1)+1);
        A((N-n)*(n+1)+i,1:N-n) = y.'*diff(Ta,a(i)).';
        A((N-n)*(n+1)+i,i*(N-n)+1:(i+1)*(N-n)) = y.'*Ta.';
    end 

    % Derive the numerical coefficient matrices:
    [monomials,supp] = powers(a);
    mat = cell(length(monomials),1);
    for i = 1:length(monomials)
        mat{i} = zeros(size(A));
    end
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            [coefs,mons] = coeffs(A(i,j));
            for k = 1:length(mons)
                idx = find(mons(k) == monomials);
                temp = mat{idx};
                temp(i,j) = coefs(k);
                mat{idx} = temp;
            end
        end
    end
end

% Define the required subfunction:
function [monomials,supp] = powers(a)
    %POWERS   Monomials and support.
    %   [monomials,supp] = POWERS(a) returns the symbolic monomials and the
    %   (numerical) support for the different powers of the variables in a.
    %   All monomials of total degree up to two are required by the main
    %   function.

    n = length(a);
    I = fliplr(eye(n,n));
    supp = zeros(nchoosek(n+2,2),n);
    supp(2:n+1,:) = I;
    blockLength = n:-1:1;
    blockPos = [0; cumsum(blockLength')];
    for i = n:-1:1
        supp(n+blockPos(i)+2:n+blockPos(i+1)+1,:) = ...
            supp(n-blockLength(i)+2:n+1,:) + I(i,:);
    end
    monomials = prod(a'.^supp,2);
end