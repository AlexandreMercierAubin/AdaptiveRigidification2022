function [td] = loadTDcellarray(filename)
    tmp = load(filename);
    td = tmp.td;
end

