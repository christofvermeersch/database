 function [system] = arma11eq(y) 
    N = length(y);
    S = [eye(N-2) zeros(N-2,1)] + [zeros(N-2,1) eye(N-2)];
    eqs = cell(N+1,1); 
    eqs{1} = [y(1:N-1) zeros(N-1,2) eye(N-1)]; 
    eqs{2} = [ones(N-1,1) zeros(N-1,1) ones(N-1,1) 2*eye(N-1); ...
        ones(N-2,1) zeros(N-2,2) S]; 
    eqs{3} = [N 0 1 0 1 zeros(1,N-3); N 0 0 1 zeros(1,N-2); N 0 2 1 zeros(1,N-2); ...
        y(2) zeros(1,N+1); y(1) 1 zeros(1,N)];
    for k = 2:N-2 
        eqs{k+2} = [N 0 1 zeros(1,k) 1 zeros(1,N-k-2); N 0 0 zeros(1,k-1) 1 zeros(1,N-k-1); ...
            N 0 2 zeros(1,k-1) 1 zeros(1,N-k-1); N 0 1 zeros(1,k-2) 1 zeros(1,N-k); y(k+1) zeros(1,N+1); y(k) 1 zeros(1,N)];
    end
    eqs{N+1} = [N zeros(1,N) 1; N 0 2 zeros(1,N-2) 1; N 0 1 zeros(1,N-3) 1 0; y(N) zeros(1,N+1); y(N-1) 1 zeros(1,N)];

    system = systemstruct(eqs); 
end