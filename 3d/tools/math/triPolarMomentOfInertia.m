function I = triPolarMomentOfInertia(p1, p2, p3)
    b = norm(p2 - p1);
    a = dot(p3 - p1, p2 - p1) / b;
    h = norm(cross([p2 - p1; 0], [p3 - p1; 0])) / b;
    I = (h * b * b * b + h * a * b * b + h * a * a * b + h * h * h * b) / 12;
end

