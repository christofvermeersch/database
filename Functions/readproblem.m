function problem = readproblem(name,flag)
    %READPROBLEM   Initialization of a problem from a CSV file.
    %   problem = READPROBLEM(name) reads the problem from the CSV file 
    %   "name.csv" and initializes the problem in MATLAB. The problem is 
    %   given in the problemstruc format of MacaulayLab.
    %
    %   problem = READPROBLEM(...,flag) also looks the problem up in the
    %   database and adds any additional information to the problem struct.
    %
    %   See also WRITEPROBLEM.
    
    % Copyright (c) 2024 - Christof Vermeersch

    % Access the problem:
    try
        file = append(name,".csv");
        info = readmatrix(file,"range","1:1");
    catch ME
        if strcmp(ME.identifier,"MATLAB:textio:textio:FileNotFound")
            error("Problem is not part of the database.");
        end
        rethrow(ME)
    end

    % Create the problemstruct:
    coef = cell(info(1),1);
    supp = cell(info(1),1);
    for eq = 1:info(1)
        % Read two lines:
        evenline = readmatrix(file,"range",sprintf("%d:%d",2*eq,2*eq));
        oddline = readmatrix(file,"range",sprintf("%d:%d",2*eq+1,2*eq+1));

        % Process the coefficients of one (matrix) polynomial:
        coeflength = sum(~isnan(evenline));
        coefelements = coeflength/(info(3)*info(4));
        evenline = evenline(1:coeflength);
        coef{eq} = reshape(evenline,[coefelements info(3) info(4)]);

        % Process the support of one (matrix) polynomial:
        supplength = sum(~isnan(oddline));
        suppelements = supplength/(info(2));
        oddline = oddline(1:supplength);
        supp{eq} = reshape(oddline, [suppelements info(2)]);
    end
    problem = problemstruct(coef,supp);

    % Add additional information from the database:
    if nargin > 1 && flag  
        info = checkproblem(name);
        problem = informproblem(problem,info);
    end
end

% Define all the required subfunctions:
function info = checkproblem(name)
    %CHECKPROBLEM   Search problem in the database.
    %   info = CHECKPROBLEM(name) checks whether the problem has an entry
    %   in the database and, if so, retrieves the information from the
    %   database.

    options = detectImportOptions("database.csv");
    options.SelectedVariableNames = 1;
    entries = readmatrix("database.csv",options);
    if any(strcmp(entries,name))
        idx = find(strcmp(entries,name)) + 1;
        options = detectImportOptions("database.csv");
        options.DataLines = [idx idx];
        info = readtable("database.csv",options);
    else
        error("Problem has no entry in the database.")
    end
end

function problem = informproblem(problem,info)
    %INFORMPROBLEM   Supplement problem with additional information.
    %   problem = INFORMPROBLEM(problem,info) adds the information from the
    %   database to the problem.

    % Load the different properties of the database:
    properties = info.Properties.VariableNames;

    % Change the problem to the specific sub-class:
    type = info.(properties{2});
    switch type{1}
        case "polynomial system"
            eqs = cell(problem.s,1);
            for i = 1:problem.s
                eqs{i} = [problem.coef{i} problem.supp{i}];
            end
            problem = systemstruct(eqs);
        case "eigenvalue problem"
            tensor = problem.coef{1};
            mat = cell(size(tensor,1),1);
            for i = 1:size(tensor,1)
                mat{i} = squeeze(tensor(i,:,:));
            end
            problem = mepstruct(mat,problem.supp{1});
    end
    
    % Add the polynomial basis to the problem:
    basis = info.(properties{3});
    switch basis{1}
        case "monomial"
            problem.basis = @monomial;
        case "chebyshev"
            problem.basis = @chebyshev;
        case "legendre"
            problem.basis = @legendre;
        otherwise
            problem.basis = "unknown";
    end

    % Add the other properties to the problem:
    for k = 4:numel(properties)
        if any(problem.(properties{k})) % check whether information agrees.
            if problem.(properties{k}) ~= info.(properties{k})
                error("Problem does not agree with database information.");
            end
        else
            problem.(properties{k}) = info.(properties{k});
        end
    end
end