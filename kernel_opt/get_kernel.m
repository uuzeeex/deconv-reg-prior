addpath(genpath(pwd));
k = kernel_estimation(im2double(pumpkins), im2double(pumpkinsblurred2), 9, 5, 'l1ls', true);