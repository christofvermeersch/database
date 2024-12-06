function writeproblem(problem,name)
    %WRITEPROBLEM   Creation of CSV file with problem definition.
    %   WRITEPROBLEM(problem,name) writes a CSV file (with user-specified
    %   name) that contains the problem definition in the database format.
    %   The problem is given in the problemstruc format of MacaulayLab.
    %
    %   See also READPROBLEM.
    
    % Copyright (c) 2024 - Christof Vermeersch
    
    % Create a new CSV file:
    file = append(name,".csv");

    % Write the header:
    [~,k,l] = size(problem.coef{1});
    info = [problem.s problem.n k l];
    writematrix(info, file);

    for eq = 1:problem.s
        % Write the coefficients:
        coeftensor = problem.coef{eq};
        coefvector = coeftensor(:);
        writematrix(coefvector.',file,"writemode","append");

        % Write the support:
        suppmatrix = problem.supp{eq};
        suppvector = suppmatrix(:);
        writematrix(suppvector.',file,"writemode","append");
    end
end