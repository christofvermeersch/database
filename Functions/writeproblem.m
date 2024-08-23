function [] = writeproblem(problem,name)
    file = append(name,'.csv');
    [~,k,l] = size(problem.coef{1});
    info = [problem.s problem.n k l];

    writematrix(info, file);

    for eq = 1:problem.s
        coeftensor = problem.coef{eq};
        coefvector = coeftensor(:);
        writematrix(coefvector', file, 'WriteMode', 'append');
        suppmatrix = problem.supp{eq};
        suppvector = suppmatrix(:);
        writematrix(suppvector', file, 'WriteMode', 'append');
    end
end