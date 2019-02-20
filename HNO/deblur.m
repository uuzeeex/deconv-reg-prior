P = psfMoffat([9, 9], 2, 1);
imagesc(P), colormap gray
X_deblur = deblur_tv_fista(im2double(pumpkinsblurred2), K, [33, 33], 0.001, -Inf, Inf);