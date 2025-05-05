function createreport(name)
    %CREATEREPORT Create a report of a problem.

    % Copyright (c) 2025 - Christof Vermeersch.

    % Access the problem:
    try
        problem = readproblem(name);
    catch ME
        if strcmp(ME.identifier, 'MATLAB:readproblem:FileNotFound')
            error('Problem file not found.');
        end
        rethrow(ME);
    end

    % Open a new file for writing the report:
    fid = fopen(strcat(name,".tex"),"w");
    printpreamble(name,fid)

    if problem.k > 1 && problem.l > 1 && problem.s == 1
        str = replace(problem.structure(),"l","\\lambda_");
        printstructure(str,problem.coefmat(),problem.supp{1},fid);
    end

    printsolutions(name,problem.n,fid);

    fprintf(fid,"\\end{document}\n");
end

% Define all the required subfunctions:
function printpreamble(name,fid)
    fprintf(fid,"\\documentclass{article}\n");
    fprintf(fid,"\\usepackage{amsmath}\n");
    fprintf(fid,"\\usepackage{float}\n");
    fprintf(fid,"\\title{%s}\n",name);
    fprintf(fid,"\\author{Christof Vermeersch}\n",name);
    fprintf(fid,"\\begin{document}\n");
    fprintf(fid,"\\maketitle\n");
end

function printstructure(structure,mat,supp,fid)
    fprintf(fid,"\\section*{Structure of the problem}\n");
    fprintf(fid,"\\begin{equation}\n");
    structure = splitlines(structure);
    fprintf(fid,structure{2});
    fprintf(fid,"\\end{equation}\n");
    fprintf(fid,structure{4});

    for i = 1:length(mat)
        printmatrices(mat{i},supp(i,:),fid);
    end
end

function printmatrices(A,supp,fid)
    rows = size(A,1);

    fprintf(fid,"\\begin{equation}\n");
    fprintf(fid,"A%d%d = ",supp(1),supp(2));
    fprintf(fid,"\\begin{bmatrix}\n");

    for i = 1:rows
        rowEntries = arrayfun(@(x) num2str(x),A(i, :),"UniformOutput",false);
        fprintf(fid,strjoin(rowEntries," & "));
        fprintf(fid," \\\\ \n");
    end

    fprintf(fid,"\\end{bmatrix}\n");
    fprintf(fid,"\\end{equation}\n");
end

function printsolutions(name,n,fid)
    fprintf(fid,"\\section*{Structure of the problem}\n");
    fprintf(fid,"\\begin{table}[H]\n");
    fprintf(fid,"\\centering\n");
    pos = join(repmat("c",1,n),"");
    fprintf(fid,"\\begin{tabular}{%s}\n",pos);
    fprintf(fid,"\\hline\n");
    for j = 1:n 
        fprintf(fid,"$\\mathbf{\\lambda_%d}$",j);

        if j < n
            fprintf(fid," & ");
        else
            fprintf(fid," \\\\ \n");
            fprintf(fid,"\\hline\n");
        end
    end

    str = fileread(strcat(name,".sol"));
    sol = splitlines(str);
    for i = 1:length(sol)
        str = split(sol{i},",");

        if length(str) == n
            for j = 1:length(str)
                fprintf(fid,str{j});

                if j < n
                    fprintf(fid," & ");
                else
                    fprintf(fid," \\\\ \n");
                end
            end
        end
    end
    fprintf(fid,"\\hline\n");
    fprintf(fid,"\\end{tabular}\n");
    fprintf(fid,"\\end{table}\n");
end
