
% H2 globally optimal model reduction using Walsh's theorem
% Theory available on: https://ftp.esat.kuleuven.be/pub/stadius/slagauw/23-47.pdf
%
% (c) Sibren Lagauw, 2023.


clear all

% Setup:
model = 'FOM2' % Higher-order model
methodology = 'MEP' % Options: system, MEP, secular
m = 2 % Order of approximant

switch model
    % b = [b_(n), ..., b_1]
    % a = [1, a_(n), ..., a_1]
    case 'Mauricio'
        % Paper: Globally Optimal H2-Norm Model Reduction: A Numerical 
        % Linear Algebra Approach, O.M.Agudelo, V. Christof, B. De Moor 
        % (2021) - IFAC-PapersOnLine.
        b = [1, 9, -10];
        a = [1, 12, 49, 78];
        isContinuous = true;
    case 'F1' 
        % Paper: Second order H2 optimal approximation of linear dynamical
        % systems, M.I. Ahmad, M. Frangos, I.M. Jaimoukha (2011) 
        % - Proc. of the 18th IFAC World Congress, Milano, Italy.
        b = [-2.9239, -39.5525,-97.5270, -147.1508]; 
        a = [1, 11.9584, 43.9119, 73.6759, 44.3821];
        isContinuous = true;
    case 'F2' 
        % Paper: Second order H2 optimal approximation of linear dynamical
        % systems, M.I. Ahmad, M. Frangos, I.M. Jaimoukha (2011) 
        % - Proc. of the 18th IFAC World Congress, Milano, Italy.
        b = [-1.2805, -6.2266, -12.8095, -9.3373];
        a = [1, 3.1855, 8.9263, 12.2936, 3.1987];
        isContinuous = true;
    case 'F3'
        % Paper: Second order H2 optimal approximation of linear dynamical
        % systems, M.I. Ahmad, M. Frangos, I.M. Jaimoukha (2011) 
        % - Proc. of the 18th IFAC World Congress, Milano, Italy.
        b = [-1.3369, -4.8341, -47.5819, -42.7285];
        a = [1, 17.0728, 84.9908, 122.4400, 59.9309];
        isContinuous = true;
    case 'F5'
        % Paper: Second order H2 optimal approximation of linear dynamical
        % systems, M.I. Ahmad, M. Frangos, I.M. Jaimoukha (2011) 
        % - Proc. of the 18th IFAC World Congress, Milano, Italy.
        b = [1, 15, 50];
        a = [1, 5, 33, 79, 50];
        isContinuous = true;
    case 'F6'
        % Paper: Second order H2 optimal approximation of linear dynamical
        % systems, M.I. Ahmad, M. Frangos, I.M. Jaimoukha (2011) 
        % - Proc. of the 18th IFAC World Congress, Milano, Italy.
        b = [41, 50, 140];
        a = conv([1, 1, 1],[1, 10, 100]);
        isContinuous = true;
    case 'FOM2'
        % Paper: H2 Model Reduction for Large-Scale Linear Dynamical 
        % Systems, S. Gugercin, A.C. Antoulas, and C. Beattie (2008)
        % - SIAM Journal on Matrix Analysis and Applications.
        b = [2, 11.5, 57.75, 178.625, 345.5, 323.625, 94.5];
        a = [1, 10, 46, 130, 239, 280, 194, 60];
        isContinuous = true;
    case 'Spanos-3'
        % Paper: A New Algorithm for L_2 Optimal Model Reduction, 
        % J.T. Spanos, M.H. Milman, D.L. Mingori (1992) - Automatica.
        b= [0, 0, 0, 0.00001, 0.0110, 1];
        a = [1, 0.2220, 22.1242, 3.5445, 122.4433,...
            11.3231, 11.1100];
        isContinuous = true;
    case 'FourDisk'
        % # Paper: H2 model reduction for SISO systems, B. De Moor, 
        % P. Van Overschee, and G. Schelfout (1993) 
        % - Proc. of the 12th IFAC World Congress, Sydney, Australia.
        b = [0.0448, 0.2368, 0.0013, 0.0211, 0.2250, 0.0219]
        a = [1, -1.2024, 2.3675, -2.0039, 2.2337, -1.0420, 0.8513]
        isContinuous = false;
end
n = length(a)-1;
roots(a)
sys = tf(b,a)

% Construct system of equations:
syms a_sym [1 m]
syms b_sym [1 m]
syms g_sym [1 n-m]
a_sym = [1, flip(a_sym)]
b_sym = flip(b_sym)
g_sym = flip(g_sym)

% Compose coefficient vector of the polynomial:
tmp1 = conv_sym(b,a_sym)
tmp2 = conv_sym(b_sym,a)

if isContinuous
    % Flip over imaginary axis:
    for i = 0:m
        j = m-i;
        if mod(j,2) ~= 0
            a_m(i+1) = -a_sym(i+1);
        else
            a_m(i+1) = a_sym(i+1);
        end
    end
else
    % Flip over unit circle:
    a_m = flip(a_sym);
end
tmp3 = conv_sym(conv_sym(a_m,a_m),g_sym);
p = addPoly(addPoly(tmp1,-tmp2),-tmp3);
[P,Var] = fcn_sym2PhilFormat(p);

% Compose MEVP:
dinit = 2; % Degree of MEP
nbEigs = m; % Number of eigenvalues
[MEP,Var_mep,EigV_mep] = SYS2MEP(P,Var,[b_sym, g_sym]);
mep = mepstruct(MEP,nbEigs,dinit);

% % Compose polynomial matrix of MEVP:
% mons = monomialsmatrix(2,length(Var_mep));
% K = zeros(size(MEP{1}));
% for i = 1:length(MEP)
%     disp('--------')
%     MEP{i};
%     vars = Var_mep.^(mons(i,:));
%     K = K + MEP{i}.*prod(vars);
% end
% vpa(K) % 'K' contains the polynomial matrix defining the MEVP.
% 
% % Solve for stationary points:
% switch methodology
%     case 'system'
%         dend = 15;
%         
%         options = struct(); 
%         options.algorithm = 'null';
%         options.blocked = true;
%         options.recursive = 'sparse';
%         options.verbose = true;
%         options.posdim = true;
%         [X,output] = solvesystem(systemstruct(P),dend,options)
%         X = X(:,1:m); % Select values for a_i only.
%         
%     case 'MEP'
%         
%         dend = 15;
%         options = struct(); 
%         options.algorithm = 'null';
%         options.blocked = true;
%         options.recursive = 'sparse';
%         options.verbose = true;
%         options.posdim = true;
%         
%         [X,output] = solvemep(mep,dend,options)
%         
%     case 'secular'
%         % CAUTION: number of secular equations grows quickly!!
%         if m == 1 
%             p = det(K);
%             X = double(roots(flip(coeffs(p))))
%         else
%             minor = sym('m', [nchoosek(n+m,m-1),1]);
%             combs = nchoosek(1:n+m,m-1);
%             for i = 1:size(combs,1)
%                 disp(['Computing minor: ',mat2str(i)])
%                 tmp = setdiff(1:m+n,combs(i,:));
%                 minor(i) = det(K(tmp,:));
%             end
% 
%             % Solve system of equations:
%             [P_minor,Var_minor] = fcn_sym2PhilFormat(minor)
%             dend = 30;
%             options = struct(); 
%             options.algorithm = 'null';
%             options.blocked = true;
%             options.recursive = 'sparse';
%             options.verbose = true;
%             options.posdim = false;
%             [X,output] = solvesystem(systemstruct(P_minor),dend,options);
%         end
% end
% 
% % Extract real solutions only:
% X = X(sum(abs(imag(X)),2) <= 1e-5,:);
% disp(['Number real solutions: ',mat2str(size(X,1))])
%        
% % Find corresponding b_i and g_i:
% B = zeros(size(X)); R = zeros(size(X,1),n-m);
% for sol = 1:size(X,1)
%     A = vpa(subs(K,flip(a_sym(2:end)),X(sol,1:m)));
%     [~,s,v] = svd(A); % v = [b_m, ..., b_1, g_i, ...]
%     assert(s(n+1,n+1) < 1e-4*s(n,n)) % Assert that A is rank-deficient.
%     v = v(:,end); % nullspace of A
%     v = v./v(1);
%     B(sol,:) = flip(v(2:m+1)); 
%     R(sol,:) = flip(v(m+2:end));
% end
% X = [X,B,R];
% 
% % Inspect solutions:
% nbSols = size(X,1);
% T_all = array2table(X);
% T_all.Properties.VariableNames = cellstr(string(Var));
% 
% % Calculate the poles of the models to inspect stability:
% poles = zeros(nbSols,m);
% stable = zeros(nbSols,1);
% for i = 1:nbSols
%     poles(i,:) = roots([1,T_all{i,string(a_sym(2:end))}]);
%     if isContinuous
%         stable(i) = max(poles(i,:)) <= 0;
%     else
%         stable(i) = max(abs(poles(i,:))) <= 1;
%     end
% end
% tmp_poles = table(poles,stable,'VariableNames',{'poles','stability'});
% disp(['Number of stable solutions: ',mat2str(sum(stable))])
% T = [T_all, tmp_poles];
% 
% % Check whether solutions satisfy the FONC of MOR:
% G = @(s) polyval(b,s)./polyval(a,s);
% clear Gr
% fonc_zeros = zeros(nbSols,m);
% fonc_poles = zeros(nbSols,m);
% for i = 1:nbSols
%     Gr{i} = @(s) polyval(T{i,string(b_sym)},s)./polyval([1,T{i,string(a_sym(2:end))}],s);
%     for j = 1:m
%         if isContinuous
%             % Optimal zeros:
%             fun = @(w) (G(1i.*w)-Gr{i}(1i.*w)).*...
%                 conj(1./(1i.*w-poles(i,j)));
%             fonc_zeros(i,j) = 1/(2*pi).*integral(@(t) fun(t),-inf,inf);
%             % Optimal poles:
%             fun = @(w) (G(1i.*w)-Gr{i}(1i.*w)).*...
%                 conj(1./((1i.*w-poles(i,j)).^2));
%             fonc_poles(i,j) = 1/(2*pi).*integral(@(t) fun(t),-inf,inf);
%         else 
%             % Optimal zeros:
%             fun = @(w) (G(exp(1i*w))-Gr{i}(exp(1i*w))).*...
%                 conj(1./(exp(1i*w)-poles(i,j)));
%             fonc_zeros(i,j) = 1/(2*pi).*integral(@(t) fun(t),-pi,pi);
%             % Optimal poles:
%             fun = @(w) (G(exp(1i*w))-Gr{i}(exp(1i*w))).*...
%                 conj(1./((exp(1i*w)-poles(i,j)).^2));
%             fonc_poles(i,j) = 1/(2*pi).*integral(@(t) fun(t),-pi,pi);
%         end
%     end
% end
% tmp = table(fonc_zeros,fonc_poles,'VariableNames',...
%     {'FONC_resid_zeros','FONC_resid_poles'});
% T = [T,tmp];
% 
% 
% % Compute values of costfunction:
% fVals = zeros(nbSols,1);
% if isContinuous
%     fun = @(w) (G(1i*w).*conj(G(1i*w)));
%     G_norm = sqrt(1/(2*pi)*integral(@(t) fun(t),-inf,inf))
%     
%     for i = 1:nbSols
%         Gr{i} = @(s) polyval(T{i,string(b_sym)},s)./polyval([1,T{i,string(a_sym(2:end))}],s);
%         fun = @(w) (Gr{i}(1i*w)-G(1i*w)).*conj(Gr{i}(1i*w)-G(1i*w)); 
%         fVals(i) = sqrt(1/(2*pi)*integral(@(t) fun(t),-inf,inf));
%     end
% else
%     fun = @(w) (G(exp(1i*w)).*conj(G(exp(1i*w))));
%     G_norm = sqrt(1/(2*pi)*integral(@(t) fun(t),-pi,pi))
%     
%     for i = 1:nbSols
%         Gr{i} = @(s) polyval(T{i,string(b_sym)},s)./polyval([1,T{i,string(a_sym(2:end))}],s);
%         fun = @(w) (Gr{i}(exp(1i*w))-G(exp(1i*w))).*conj(Gr{i}(exp(1i*w))-G(exp(1i*w))); 
%         fVals(i) = sqrt(1/(2*pi)*integral(@(t) fun(t),-pi,pi));
%     end
% end
% 
% tmp = table(fVals,fVals./G_norm, 'VariableNames', {'error','rel_error'});
% tmp2 = table((1:nbSols).','VariableNames',{'Model'});
% T = [tmp2,tmp,T];
% 
% open('T')

%% Auxiliary functions:
function res = addPoly(in1,in2)
    % Make sure that vectors are given as col vectors:
    %
    % (c) Sibren Lagauw, 2022.
    
    [~,p] = size(in1);
    if p ~= 1
        in1 = in1.';
    end
    [~,p] = size(in2);
    if p ~= 1
        in2 = in2.';
    end
    
    n = length(in1);
    m = length(in2);
    if m < n
        in2 = [zeros(n-m,1);in2];
    else
        in1 = [zeros(m-n,1);in1];
    end
    res = in1+in2;
end

function res = conv_sym(in1,in2)
    % Calculates the convolution of two vectors (elements can be complex).
    %
    % (c) Sibren Lagauw, 2022.

    % Make sure that vectors are given as col vectors:
    [~,p] = size(in1);
    if p ~= 1
        in1 = in1.';
    end
    [~,p] = size(in2);
    if p ~= 1
        in2 = in2.';
    end

    n = length(in1);
    m = length(in2);
    l = n+m-1;
    
    % Pad with zeros:
    in1 = [in1;zeros(l-n,1)];
    in2 = [in2;zeros(l-m,1)];
    
    res = cell(l,1);
    for i = 1:l
       tmp = 0;
        for j = 1:i
            tmp = tmp + in1(j)*in2(i+1-j);
        end
        res{i} = tmp;
    end
    res = cell2sym(res);
end

function [Mep,Var_mep,EigV] = SYS2MEP(P,Var,EigV)
% This function converts a system of polynomial equations (square
% system =>  number of equations = number of unknowns) to the MEP format used 
% in Christof's Toolbox, by extracting a given set of variables. 

% CAUTION: the set of EigV must be chosen carefully, the equations can only
% contain terms that are linear expressions in the EigV parameters! 
% E.g., if l1 and l2 only appear in sys with degree one, both
% variables will can be extracted to the eigenvector. However, if sys contains 
% a term l1*l2 than only one of both variables can be extracted!
%
% INPUT:
% P     = Cell containing monomials which compose the multivariate polynomial equations 
% Var   = All variables in P
% EigV  = Variables in the eigenvector of the MEP
%
% OUTPUT:
% Mep =  System of polynomial equations  in Philippe's format
% Var =  Vector containing the variables of the polynomials (Symbolic
% expressions)
% EigV = Vector containing the variables in the eigenvector of the MEP
%
% Sibren Lagauw
% (c) 2021
%

npol =length(P); % Number of eqs.
nvar=length(Var); % Number of vars.

%if (nvar~=npol)
%   Mep={}; Var = []; EigV = []; 
%   warning('The number of variables is not equal to the number of equations. An empty cell array has been returned.');
%   return
%end

Var_mep = setdiff(Var, EigV);

% Find max degree of each variable:
d_max = zeros(nvar,1);
for poly_index = 1:npol
    poly = P{poly_index};
    for mon_index = 1:size(poly,1)
        for var_index = 1:nvar
            d_max(var_index) = max(d_max(var_index), poly(mon_index,1+var_index));
        end
    end
end
Var_mep_indices = zeros(1,length(Var_mep));
for i = 1:length(Var_mep)
    Var_mep_indices(i) = find(Var == Var_mep(i));
end
EigV_indices = zeros(1,length(EigV));
for i = 1:length(EigV)
    EigV_indices(i) = find(Var == EigV(i));
end
    
% Create monomials (both for MEP-matrices and eigenvector):
monomials_MEP = monomialsmatrix(max(d_max),length(Var_mep));
monomials_EigV = monomialsmatrix(1,length(EigV));

% Compose matrices of MEP one by one: 
Mep = cell(1,size(monomials_MEP,1));
for matrix_index = 1:length(monomials_MEP(:,1))
    M = zeros(npol,length(EigV)+1);
 
    for poly_index = 1:npol
        pol = P{poly_index}; 
        
        for mon_index = 1:size(pol,1)
            % Check whether this is the correct matrix for this monomial:
            tmp = pol(mon_index,Var_mep_indices+1);
            if monomials_MEP(matrix_index,:) == tmp
                % Find correct column of matrix: 
                mon_eig = pol(mon_index, EigV_indices+1);
                [tf_2, index_2] = ismember(mon_eig,monomials_EigV,'rows');
                if ~tf_2 % Should be always true and unique
                    warning('Something went wrong')
                    return
                else
                   % Fill entry in matrix
                   M(poly_index,index_2) = double(pol(mon_index,1));
                end
            end
        end
    end
    
    Mep{matrix_index} = M;
end
Mep = Mep.'; % Adjustment to comply with format of new version Maclab-internal

EigV = [1, EigV];

Var_mep = Var_mep;

end

function [Pol,Var]=fcn_sym2PhilFormat(P)
% INPUT:
% P   =  Vector containing symbolic multivariate polynomial equations 
%
% OUTPUT:
% Pol =  System of polynomial equations  in Philippe's format
% Var =  Vector containing the variables of the polynomials (Symbolic
% expressions)
%
% Example:
% 
% syms x1 x2
% f(1) = 2*x1^2 - x2;
% f(2) = 3*x1- 4*x2 + 5;
% [P] = fcn_sym2PhilFormat(f);  
%
% In this case P is a cell array with the following structure:
% 
% First polynomial
% P{1} = [ 2 2 0 ] %Coefficient, power of x1, power of x2 : first monomial
%        [-1 0 1]  %Coefficient, power of x1, power of x2 : Second monomial 
%
% Second polynomial
%P{2} =   [ 3  1  0] %Coefficient, power of x1, power of x2 : first monomial
%         [ -4 0  1] %Coefficient, power of x1, power of x2 : Second monomial 
%         [ 5 0  0 ] %%Coefficient, power of x1, power of x2 : Third monomial 
%
% The order of the variables is x1, x2.
%
%
% Oscar Mauricio Agudelo
% (c) 2018
%

npol = length(P); %Number of polynomial equations
Var = symvar(P);
nvar=length(Var);
Var_str = char(Var);
%Var_str([1:8, end-1:end]) = []; %removing the extra characteres such as "matrix(" after the conversion with "char"

% if (nvar~=npol)
%   Pol={}; 
%   warning('The number of variables is not equal to the number of equations. An empty cell array has been returned.');
%   return
% end    

for k=1:npol
   t=feval(symengine,'poly2list',vpa(P(k)),Var_str);
   n_mon = length(t); %number of monomials
    ma = zeros(n_mon,nvar+1);
   for j=1:n_mon
       monomial = t(j);
       ma(j,1) = monomial(1);
       ma(j,2:end) = monomial(2);
   end 
   Pol{k} = ma;
end    

end

