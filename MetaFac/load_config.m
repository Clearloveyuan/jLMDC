function opts = load_config(filename)
% LOAD_CONFIG - read configuration file

[keynames,values]=textread(filename,'%s=%s', 'commentstyle', 'c++');
v = str2double(values); idx = ~isnan(v);
values(idx) = num2cell(v(idx));
opts = cell2struct(values, keynames);