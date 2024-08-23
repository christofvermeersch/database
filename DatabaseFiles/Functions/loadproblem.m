function [problem,info] = loadproblem(name)
    [table] = loaddatabase();
    mask = strcmp(table.id,name);
    if any(mask)
        info = table(mask,:);
        info = table2struct(info);
        problem = getproblem(name,info);
    else
        error('This problem is not part of the database!')
    end
end

function [problem] = getproblem(name,info)
    problem = readproblem(name);
    fn = fieldnames(info);
    for k = 1:numel(fn)
        problem.(fn{k}) = info.(fn{k});
    end
end