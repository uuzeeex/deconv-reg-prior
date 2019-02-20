function [X_out, fun_all] = deblur_tv_fista(Bobs, P, center, lambda, l, u, pars)

flag = exist('pars', 'var');
if (flag && isfield(pars, 'MAXITER'))
  MAXITER = pars.MAXITER;
else
  MAXITER = 100;
end
if(flag && isfield(pars, 'fig'))
  fig = pars.fig;
else
  fig = 1;
end
if (flag && isfield(pars, 'BC'))
  BC = pars.BC;
else
  BC = 'reflexive';
end
if (flag && isfield(pars, 'tv'))
  tv = pars.tv;
else
  tv = 'iso';
end
if (flag && isfield(pars, 'mon'))
  mon = pars.mon;
else
  mon = 0;
end
if (flag && isfield(pars, 'denoiseiter'))
  denoiseiter = pars.denoiseiter;
else
  denoiseiter = 10;
end

if (nargout == 2)
  fun_all = [];
end

[m, n] = size(Bobs);
Pbig = padPSF(P, [m, n]);

switch BC
  case 'reflexive'
    trans = @(X)dct2(X);
    itrans = @(X)idct2(X);
    e1 = zeros(m, n);
    e1(1, 1) = 1;
    Sbig = dct2(dctshift(Pbig, center)) ./ dct2(e1);
  case 'periodic'
    trans = @(X)(1/sqrt(m * n) * fft2(X));
    itrans = @(X)(sqrt(m * n) * ifft2(X));
    Sbig = fft2(circshift(Pbig, 1 - center));
  otherwise
    error('Invalid boundary conditions should be reflexive or periodic');
end

Btrans = trans(Bobs);

L = 2 * max(max(abs(Sbig) .^ 2));

clear parsin
parsin.MAXITER = denoiseiter;
parsin.epsilon = 1e-5;
parsin.print = 0;
parsin.tv = tv;

X_iter = Bobs;
Y = X_iter;
t_new = 1;

fprintf('***********************************\n');
fprintf('*   Solving with FISTA      **\n');
fprintf('***********************************\n');
fprintf('#iter  fun-val         tv          denoise-iter      relative-dif\n===============================================\n');
for i = 1 : MAXITER
  X_old = X_iter;
  t_old = t_new;
  D = Sbig .* trans(Y) - Btrans;
  Y = Y - 2 / L * itrans(conj(Sbig) .* D);
  Y = real(Y);
  if (i == 1)
    [Z_iter, iter, ~, P] = denoise_bound_init(Y, 2 * lambda / L, l, u, [], parsin);
  else
    [Z_iter, iter, ~, P] = denoise_bound_init(Y, 2 * lambda / L, l, u, P, parsin);
	end
  t = tlv(Z_iter, tv);
  fun_val = norm(Sbig .* trans(Z_iter) - Btrans, 'fro') ^ 2 + 2 * lambda * t;
  if (mon == 0)
    X_iter = Z_iter;
	elseif (i>1)
    fun_val_old = fun_all(end);
    if (fun_val > fun_val_old)
      X_iter = X_old;
      fun_val = fun_val_old;
    else
      X_iter = Z_iter;
    end
	end
  if (nargout == 2)
    fun_all = [fun_all; fun_val];
	end
  t_new = (1 + sqrt(1 + 4 * t_old ^ 2)) / 2;
  Y = X_iter + t_old / t_new * (Z_iter - X_iter) + (t_old - 1) / t_new * (X_iter - X_old);
  fprintf('%3d    %15.5f %15.5f           %3d                  %15.5f\n', i, fun_val, t, iter, norm(X_iter - X_old, 'fro') / norm(X_old, 'fro'));
  %if (fig)
  %  figure(314)
  %  imshow(X_iter, [])
  %end
end
X_out = X_iter;
