function [] = addproblem(problem,name)
    warning('off','backtrace')
    warning('This method is only for personal use!')
    warning('on','backtrace')

    cd ~/Documents/Workspace/Database
    try 
       loadproblem(name);
    catch ME
        % An error indicates that the problem is not yet in the database.
        if strcmp(ME.message,'This problem is not part of the database!')
            appenddatabase(problem,name)

            if strcmp(problem.type,'ps')
                cd ~/Documents/Workspace/Database/Problems/PolynomialSystems/
            elseif strcmp(problem.type,'ep')
                cd ~/Documents/Workspace/Database/Problems/EigenvalueProblems/
            else
                error('Unknown problem type!')
            end
            writeproblem(problem,name)
        else
            rethrow(ME)
        end
    end
end

function [] = appenddatabase(problem,name)
    [~,columns] = loaddatabase;
    line = cell(1,length(columns)+1);
    line{1} = name;
    for k = 1:length(columns)
        line{k+1} = problem.(columns{k});
    end
    writecell(line, 'database.csv', 'WriteMode', 'append');
end