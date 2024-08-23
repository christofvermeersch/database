function [table,names] = loaddatabase()
    % TODO: only select some columns
    options = detectImportOptions('database.csv', ReadRowNames=true);
    options = setvartype(options,{'type','basis','appid','appid','bibid','comments'},'string');
    options = setvartype(options,{'multi','posdim','solved'},'logical');
    options = setvaropts(options,{'multi','posdim','solved'}, 'TrueSymbols',{'true','1'},'FalseSymbols',{'false', '0'});
    table = readtable('database.csv', options);
    names = table.Properties.VariableNames;
end