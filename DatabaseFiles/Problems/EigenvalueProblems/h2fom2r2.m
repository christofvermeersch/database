function [mep] = h2fom2r2()
    b = [2, 11.5, 57.75, 178.625, 345.5, 323.625, 94.5];
    a = [1, 10, 46, 130, 239, 280, 194, 60];

    mat = cell(6,1);
    mat{1} = [b 0 0; -a 0; 0 -a; -1*eye(5) zeros(5,4)]';
    mat{2} = [zeros(1,2) b; zeros(2,9); zeros(5,2) -2*eye(5) zeros(5,2)]';
    mat{3} = [zeros(1,1) b 0; zeros(2,9); zeros(5,1) 2*eye(5) zeros(5,3)]';
    mat{4} = [zeros(3,9); zeros(5,4) -1*eye(5)]';
    mat{5} = [zeros(3,9); zeros(5,3) 2*eye(5) zeros(5,1)]';
    mat{6} = [zeros(3,9); zeros(5,2) -1*eye(5) zeros(5,2)]';
    mep = mepstruct(mat,2,2);
end