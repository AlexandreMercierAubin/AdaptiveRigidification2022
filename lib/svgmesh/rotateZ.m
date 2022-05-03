function [prx, pry] = rotateZ(px, py, th)
    for ii = 1:length(px)
        H = [cos(th),   -sin(th),   0; ...
            sin(th),    cos(th),   0; ...
            0,                0,    1; ...
            ];
        T = H * [px(ii); py(ii); 1];
        prx(ii) = T(1);
        pry(ii) = T(2);
    end
end