addpath(genpath(pwd));
PSF_init = P;
[n, m, d] = size(B);

P = zeros(9, 9, d);

for i = 1 : d
  P(:, :, i) = psfMoffat([9, 9], 2, 1);
end

X_deblur = zeros(n, m, d);
X_deblur_prev = zeros(n, m, d);
P_prev = zeros(9, 9, d);
X_delta = zeros(19, d);
P_delta = zeros(19, d);
for i = 1 : 20
  if i == 3
    PSF_iter3 = P;
  end
  if i == 7
    PSF_iter7 = P;
  end
  if i == 10
    PSF_iter10 = P;
  end
  P_prev = P;
  for k = 1 : d
    X_deblur(:, :, k) = deblur_tv_fista(B(:, :, k), P(:, :, k), [5, 5], 0.001, -Inf, Inf);
    P(:, :, k) = kernel_estimation(X_deblur(:, :, k), B(:, :, k), 9, 5, 'l1ls', true);
  end
  if i < 10
    for j = 1 : d
      P_delta(i, j) = norm(P(:, :, j) - P_prev(:, :, j), 'fro') ./ norm(P_prev(:, :, j), 'fro');
    end
  end
  if i > 1
    for j = 1 : d
      X_delta(i - 1, j) = norm(X_deblur(:, :, j) - X_deblur_prev(:, :, j), 'fro') ./ norm(X_deblur_prev(:, :, j), 'fro');
    end
  end
  X_deblur_prev = X_deblur;
  if i == 3
    X_iter3 = X_deblur;
  end
  if i == 7
    X_iter7 = X_deblur;
  end
end
%imwrite(B, 'jzy_blur.png')
%imwrite(X_iter3, 'jzy_iter3.png')
%imwrite(X_iter7, 'jzy_iter7.png')
%imwrite(X_deblur, 'jzy_deblurred.png')

P4print = P ./ max(max(P));
