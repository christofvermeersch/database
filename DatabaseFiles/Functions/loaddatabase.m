function [table,names] = loaddatabase()
    options = detectImportOptions('database.csv', ReadRowNames=true);
    options = setvartype(options,{'type','basis','appid','appid','bibid','comments'},'string');
    options = setvartype(options,{'multi','posdim','solved'},'logical');
    table = readtable('database.csv', options);
    names = table.Properties.VariableNames;
end