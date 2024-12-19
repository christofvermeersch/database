function [mat,supp] = h2sisored1(an,bn,r)
    %H2SISORED1   H2-norm model order reduction problem.
    %   mat = H2SISORED1(an,bn,r) creates the coefficient matrices of the
    %   multiparameter eigenvalue problem that solves the H2-norm model 
    %   order reduction problem of an n-th order model given by its 
    %   coefficients of the transfer function (bn/an) to an r-th order 
    %   model. It uses the so-called "first approach" from [1].
    %
    %   [mat,supp] = H2SISORED1(...) also gives the support of the 
    %   coefficient matrices.
        
    % Copyright (c) 2024 - Christof Vermeersch
    %
    % This function implements the construction of the numerical
    % coefficient matrices as described in [2], which is the homogeneous
    % adaptation of the inhomogeneous approach in [1].
    %
    % [1] Agudelo, O. M., Vermeersch, C., and De Moor, B. (2021)  Globally
    % Optimal H2-norm Model Reduction: A Numerical Linear Algebra Approach,
    % IFAC-PapersOnLine, Part of special issue: 24th International 
    % Symposium on Mathematical Theory of Networks and Systems (MTNS), vol.
    % 54, no. 9, p. 564-571.
    % [2] Vermeersch, C. (2023) The (Block) Macaulay Matrix: Solving
    % Systems of Multivariate Polynomial Equations and Multiparameter
    % Eigenvalue Problems, PhD thesis, KU Leuven, Leuven, Belgium.

    % Define the (symbolic) variables:
    n = length(an);
    ar = sym("ar",[r 1]);
    for i = 1:r
        assume(ar(i),"real")
    end
    br = sym("br",[r 1]);
    for i = 1:r
        assume(br(i),"real")
    end

    % Create the auxiliary matrices:
    An = [-an.'; eye(n-1) zeros(n-1,1)];
    Bn = zeros(n,1); Bn(1) = 1;
    Cn = bn.';
    Ar = [-ar.'; eye(r-1) zeros(r-1,1)];
    Br = zeros(r,1); Br(1) = 1;
    Cr = br.';
    Tr = kron(Ar,eye(r)) + kron(eye(r),Ar);
    Fr = Br*Br.';
    fr = Fr(:);
    Gr = Cr.'*Cr;
    gr = Gr(:);
    Tm = kron(An,eye(r)) + kron(eye(n),Ar);
    Fm = Br*Bn.';
    fm = Fm(:);
    Gm = Cr.'*Cn;
    gm = Gm(:);
    
    % Construct matrix polynomial:
    A = sym("A",[r^3+r^2*(n+1)+r*(n+2) r^3+r^2*(n+1)+r*n+1]);
    A(1:r,1:r^3) = -kron(eye(r),gr.');
    A(1:r,r^3+1:r^3+r^2*n) = -kron(eye(r),gm.');
    A(1:r,r^3+r^2*n+1:end) = zeros(r,r^2+r*n+1);
    for i = 1:r
        A(r+i,1:r^3+r^2*n) = zeros(1,r^3+r^2*n);
        A(r+i,r^3+r^2*n+1:r^3+r^2*n+r^2) = -diff(gr,br(i)).';
        A(r+i,r^3+r^2*n+r^2+1:r^3+r^2*n+r^2+r*n) = -2*diff(gm,br(i)).';
        A(r+i,end) = 0;
    end
    A(2*r+1:(2+r^2)*r,1:r^3) = kron(eye(r),Tr);
    A(2*r+1:(2+r^2)*r,r^3+1:r^3+r^2*n) = zeros(r^3,r^2*n);
    A(2*r+1:(2+r^2)*r,r^3+r^2*n+r^2+1:end) = zeros(r^3,r*n+1);
    for i = 1:r
        A(2*r+(i-1)*r^2+1:2*r+i*r^2,r^3+r^2*n+1:r^3+r^2*n+r^2) = ...
            diff(Tr,ar(i));
    end
    A(2*r+r^3+1:2*r+r^3+r^2*n,1:r^3) = zeros(r^2*n,r^3);
    A(2*r+r^3+1:2*r+r^3+r^2*n,r^3+1:r^3+r^2*n) = kron(eye(r),Tm);
    A(2*r+r^3+1:2*r+r^3+r^2*n,r^3+r^2*n+1:r^3+r^2*n+r^2) = ...
        zeros(r^2*n,r^2);
    A(2*r+r^3+1:2*r+r^3+r^2*n,end) = zeros(r^2*n,1);
    for i = 1:r
        A(2*r+r^3+(i-1)*r*n+1:2*r+r^3+i*r*n,r^3+r^2*n+r^2+1:end-1) = ...
            diff(Tm,ar(i));
    end
    A(2*r+r^3+r^2*n+1:2*r+r^3+r^2*n+r^2,:) = [zeros(r^2,r^3+r^2*n) Tr ...
        zeros(r^2,r*n) fr];
    A(2*r+r^3+r^2*n+r^2+1:2*r+r^3+r^2*n+r^2+r*n,:) = ...
        [zeros(r*n,r^3+r^2*n+r^2) Tm fm];

    % Derive the numerical coefficient matrices:
    [monomials,supp] = powers(ar,br);
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
function [monomials,supp] = powers(ar,br)
    %POWERS   Monomials and support.
    %   [monomials,supp] = POWERS(a,c) returns the symbolic monomials and 
    %   the (numerical) support for the different powers of the variables 
    %   in a and c. All monomials in c of total degree up to two, but only 
    %   the linear monomials in a, are required by the main function.

    n = length(br);
    I = fliplr(eye(n,n));
    supp = zeros(nchoosek(n+2,2),n);
    supp(2:n+1,:) = I;
    blockLength = cumsum(ones(n,1),"reverse");
    blockPos = [0; cumsum(blockLength)];
    for ni = n:-1:1
        supp(n+blockPos(ni)+2:n+blockPos(ni+1)+1,:) = ...
            supp(n-blockLength(ni)+2:n+1,:) + I(ni,:);
    end
    monomials = [prod(br'.^supp,2); ar];
    supp = [zeros(size(supp,1),length(ar)) supp; ...
        eye(length(ar)) zeros(length(ar),length(br))];
end