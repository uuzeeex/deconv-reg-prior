function P = psfMoffat(dim, s, beta)
m = dim(1);
n = dim(2);
x = -fix(n / 2) : ceil(n / 2) - 1;
y = -fix(m / 2) : ceil(m / 2) - 1;
[X, Y] = meshgrid(x, y);
P = exp(1 + (X .^ 2) / (2 * s ^ 2) + (Y .^ 2) / (2 * s ^ 2)) .^ (-beta);
P = P / sum(P(:));
end