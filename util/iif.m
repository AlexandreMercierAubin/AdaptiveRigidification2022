function out = iif(cond, a, b)
    %IIF implements a ternary operator
    if cond
        out = a;
    else
        out = b;
end