% ÒıÈëBine»ìãçÓ³Éäº¯Êı
function x = chebyshev_map(x0, a)
    x = zeros(1,x0);
    x(1) = rand;
    for i = 2:x0
        x(i) = a * sin(pi * x(i-1));
    end
end