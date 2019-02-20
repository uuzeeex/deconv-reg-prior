addpath(genpath(pwd));

rng(2333);

std = 1e-2;

%X = im2double(imread('Cameraman256.png'));
%X = im2double(imread('Lena512.png'));
%X = im2double(imread('image_House256rgb.png'));
X = im2double(imread('image_Peppers512rgb.png'));

[n, m, d] = size(X);

K = zeros(9, 9, d);

for i = 1 : d
  K(:, :, i) = psfGauss([9, 9]);
end

X_lr = fliplr(X);
X_ud = flipud(X);
X_x = fliplr(X_ud);

X_ext = [X_x  X_ud X_x
         X_lr X    X_lr
         X_x  X_ud X_x];

B_ext = zeros(size(X_ext));

for i = 1 : d
  B_ext(:, :, i) = conv2(X_ext(:, :, i), K(:, :, i), 'same') + std * randn(size(X_ext(:, :, i)));
end

B = B_ext(n + 1 : 2 * n, n + 1 : 2 * n, :);