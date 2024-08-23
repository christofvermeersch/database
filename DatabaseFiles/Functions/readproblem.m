function problem = readproblem(name)
    file = append(name,'.csv');
    info = readmatrix(file, 'Range', '1:1');

    coef = cell(info(1),1);
    supp = cell(info(1),1);
    for eq = 1:info(1)
        coefvector = readmatrix(file, 'Range', sprintf('%d:%d',2*eq,2*eq));
        coeflength = sum(~isnan(coefvector));
        coefelements = coeflength/(info(3)*info(4));
        coefvector = coefvector(1:coeflength);
        coef{eq} = reshape(coefvector, [coefelements info(3) info(4)]);
        suppvector = readmatrix(file, 'Range', sprintf('%d:%d',2*eq+1,2*eq+1));
        supplength = sum(~isnan(suppvector));
        suppelements = supplength/(info(2));
        suppvector = suppvector(1:supplength);
        supp{eq} = reshape(suppvector, [suppelements info(2)]);
    end
    problem = problemstruct(coef,supp);
end