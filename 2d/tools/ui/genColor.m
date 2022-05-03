function color = genColor(id)
    if id == 2
        color = [1 0 0];
        return;
    end
    rng(id);
    color = rand(1, 3) / 3 + 2 / 3;
end

