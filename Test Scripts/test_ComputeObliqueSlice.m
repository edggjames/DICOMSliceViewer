clc
close all 
clearvars

% MPHYGB24 - MATLAB coursework assignment 2017/18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Script for Task 6
% Tests function ComputeObliqueSlice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NB Script requires 'Bioinformatics Toolbox' for function suptitle.m,
% first used in line 124.

% Load all slices of 3D DICOM volume using function from Task 1:
z_limits(1) = 1;
z_limits(2) = 181;

try 

Image = LoadDICOMVolume(z_limits);

% To assign 3D image intensities to a matrix variable
vol = Image.ImageData;
% To extract voxel dimensions in mm (used to scale axes in plots)
vox_dim = Image.VoxelDimensions; % [dy dx dz]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Testing XY plane extraction, Nearest Neighbour Interpolation, High
% Resolution Sampling, and Single Rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set orientation of plane & interpolation technique
orientation = 'X-Y';
method = 'nearest';

% Set Z slice position (specified in mm):
position = 90;

% Set both slice pixel dimensions to be the smallest in-plane pixel
% dimension of the original image, i.e. testing high resolution sampling
res_1 = 0.9765625;
res_2 = 0.9765625;

% define rotation angles for first oblique slice
alpha = 45; % around x-axis
beta = 0; % around y-axis
gamma = 0; % around z-axis

% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);

% compute orthogonal slice (to plot for comparison)
orthogonal_slice = ComputeOrthogonalSliceBlur(Image, orientation, position, ...
    res_1, res_2, method);

% plot the orthogonal and oblique slices in one subplot
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
imagesc(orthogonal_slice)
colormap gray
xlabel('X')
ylabel('Y')
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title({'Normal View of Orthogonal Slice in XY Plane';''})

% plot the first oblique slice
subplot(2,2,2)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(2) rotated_vox_dim(1) rotated_vox_dim(3)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);

% calculate and plot the second oblique slice
% define rotation angles for second oblique slice
alpha = 0; beta = 45; gamma = 0; 
% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);
% plot the second oblique slice
subplot(2,2,3)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(2) rotated_vox_dim(1) rotated_vox_dim(3)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);

% calculate and plot the third oblique slice
% define rotation angles for second oblique slice
alpha = 0; beta = 0; gamma = 45; 
% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);
% plot the first oblique slice
subplot(2,2,4)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(2) rotated_vox_dim(1) rotated_vox_dim(3)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);
suptitle({'Testing Oblique Slice Computation Using XY Plane as Reference - High Resolution Sampling';
    'Nearest-Neighbour Interpolation'; ['(Z Slice Position = ', num2str(position), ' mm)']})
% save image as .jpg file with appropriate name
saveas(gcf,'XY_high_res','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Testing XY plane extraction, Nearest Neighbour Interpolation, Low
% Resolution Sampling, and Multiple Rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set both slice pixel dimensions to be four times the largest in-plane pixel
% dimension of the original image, i.e. testing low resolution sampling
res_1 = 4 * 0.9765625;
res_2 = 4 * 0.9765625; 

% define rotation angles for first oblique slice
alpha = 30; % around x-axis
beta = 30; % around y-axis
gamma = 0; % around z-axis

% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);

% compute orthogonal slice (to plot for comparison)
orthogonal_slice = ComputeOrthogonalSliceBlur(Image, orientation, position, ...
    res_1, res_2, method);

% plot the orthogonal and oblique slices in one subplot
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
imagesc(orthogonal_slice)
colormap gray
xlabel('X')
ylabel('Y')
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title({'Normal View of Orthogonal Slice in XY Plane';''})

% plot the first oblique slice
subplot(2,2,2)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(2) rotated_vox_dim(1) rotated_vox_dim(3)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);

% calculate and plot the second oblique slice
% define rotation angles for second oblique slice
alpha = 0; beta = 30; gamma = 30; 
% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);
% plot the second oblique slice
subplot(2,2,3)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(2) rotated_vox_dim(1) rotated_vox_dim(3)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);

% calculate and plot the third oblique slice
% define rotation angles for second oblique slice
alpha = 30; beta = 30; gamma = 30; 
% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);
% plot the third oblique slice
subplot(2,2,4)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(2) rotated_vox_dim(1) rotated_vox_dim(3)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);
suptitle({'Testing Oblique Slice Computation Using XY Plane as Reference - Low Resolution Sampling';
    'Nearest-Neighbour Interpolation'; ['(Z Slice Position = ', num2str(position), ' mm)']})
% save image as .jpg file with appropriate name
saveas(gcf,'XY_low_res','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Testing YZ plane extraction, Linear Interpolation, High Resolution
% Sampling, and Single Rotations. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set both slice pixel dimensions to be equal to the smallest in-plane pixel
% dimension of the original image
res_1 = 0.9765625;
res_2 = 0.9765625;

% Set orientation of plane & interpolation technique
orientation = 'Y-Z';
method = 'linear';

% Set X slice position (specified in mm):
position = 250;

% define rotation angles for first oblique slice
alpha = 45; % around x-axis
beta = 0; % around y-axis
gamma = 0; % around z-axis

% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);

% compute orthogonal slice (to plot for comparison)
orthogonal_slice = ComputeOrthogonalSliceBlur(Image, orientation, position, ...
    res_1, res_2, method);

% plot the orthogonal and oblique slices in one subplot
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
imagesc(orthogonal_slice)
colormap gray
xlabel('Z')
ylabel('Y')
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title({'Normal View of Orthogonal Slice in YZ Plane';''})

% plot the first oblique slice
subplot(2,2,2)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(1) rotated_vox_dim(3) rotated_vox_dim(2)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);

% calculate and plot the second oblique slice
% define rotation angles for second oblique slice
alpha = 0; beta = 45; gamma = 0; 
% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);
% plot the second oblique slice
subplot(2,2,3)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(1) rotated_vox_dim(3) rotated_vox_dim(2)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);

% calculate and plot the third oblique slice
% define rotation angles for second oblique slice
alpha = 0; beta = 0; gamma = 45; 
% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);
% plot the third oblique slice
subplot(2,2,4)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(1) rotated_vox_dim(3) rotated_vox_dim(2)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);
suptitle({'Testing Oblique Slice Computation Using YZ Plane as Reference - High Resolution Sampling';
    'Linear Interpolation'; ['(X Slice Position = ', num2str(position), ' mm)']})
% save image as .jpg file with appropriate name
saveas(gcf,'YZ_high_res','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Testing YZ plane extraction, Linear Interpolation, Low Resolution
% Sampling, and Multiple Rotations. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set both slice pixel dimensions to be equal to four times the largest 
% in-plane pixel dimension of the original image
res_1 = 4;
res_2 = 4;

% define rotation angles for first oblique slice
alpha = 0; % around x-axis
beta = 30; % around y-axis
gamma = 30; % around z-axis

% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);

% compute orthogonal slice (to plot for comparison)
orthogonal_slice = ComputeOrthogonalSliceBlur(Image, orientation, position, ...
    res_1, res_2, method);

% plot the orthogonal and oblique slices in one subplot
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
imagesc(orthogonal_slice)
colormap gray
xlabel('Z')
ylabel('Y')
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title({'Normal View of Orthogonal Slice in YZ Plane';''})

% plot the first oblique slice
subplot(2,2,2)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(1) rotated_vox_dim(3) rotated_vox_dim(2)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);

% calculate and plot the second oblique slice
% define rotation angles for second oblique slice
alpha = 30; beta = 0; gamma = 30; 
% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);
% plot the second oblique slice
subplot(2,2,3)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(1) rotated_vox_dim(3) rotated_vox_dim(2)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);

% calculate and plot the third oblique slice
% define rotation angles for second oblique slice
alpha = 30; beta = 30; gamma = 30; 
% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);
% plot the third oblique slice
subplot(2,2,4)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(1) rotated_vox_dim(3) rotated_vox_dim(2)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);
suptitle({'Testing Oblique Slice Computation Using YZ Plane as Reference - Low Resolution Sampling';
    'Linear Interpolation'; ['(X Slice Position = ', num2str(position), ' mm)']})
% save image as .jpg file with appropriate name
saveas(gcf,'YZ_low_res','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5) Testing XZ plane extraction, Cubic Spline Interpolation, High Resolution
% Sampling and Single Rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set orientation of plane & interpolation technique
orientation = 'X-Z';
method = 'spline';

% Set both slice pixel dimensions to be equal to the smallest in-plane pixel
% dimension of the original image
res_1 = 0.9765625;
res_2 = 0.9765625;

% define rotation angles for first oblique slice
alpha = 45; % around x-axis
beta = 0; % around y-axis
gamma = 0; % around z-axis

% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);

% compute orthogonal slice (to plot for comparison)
orthogonal_slice = ComputeOrthogonalSliceBlur(Image, orientation, position, ...
    res_1, res_2, method);

% plot the orthogonal and oblique slices in one subplot
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
imagesc(orthogonal_slice)
colormap gray
xlabel('Z')
ylabel('X')
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title({'Normal View of Orthogonal Slice in XZ Plane';''})

% plot the first oblique slice
subplot(2,2,2)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(2) rotated_vox_dim(3) rotated_vox_dim(1)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);

% calculate and plot the second oblique slice
% define rotation angles for second oblique slice
alpha = 0; beta = 45; gamma = 0; 
% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);
% plot the second oblique slice
subplot(2,2,3)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(2) rotated_vox_dim(3) rotated_vox_dim(1)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);

% calculate and plot the third oblique slice
% define rotation angles for second oblique slice
alpha = 0; beta = 0; gamma = 45; 
% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);
% plot the third oblique slice
subplot(2,2,4)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(2) rotated_vox_dim(3) rotated_vox_dim(1)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);
suptitle({'Testing Oblique Slice Computation Using XZ Plane as Reference - High Resolution Sampling';
    'Cubic Spline Interpolation'; ['(Y Slice Position = ', num2str(position), ' mm)']})
% save image as .jpg file with appropriate name
saveas(gcf,'XZ_high_res','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6) Testing XZ plane extraction, Cubic Spline Interpolation, Low Resolution
% Sampling and Multiple Rotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set both slice pixel dimensions to be equal to four times the largest 
% in-plane pixel dimension of the original image
res_1 = 4;
res_2 = 4;

% define rotation angles for first oblique slice
alpha = 30; % around x-axis
beta = 0; % around y-axis
gamma = 30; % around z-axis

% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);

% compute orthogonal slice (to plot for comparison)
orthogonal_slice = ComputeOrthogonalSliceBlur(Image, orientation, position, ...
    res_1, res_2, method);

% plot the orthogonal and oblique slices in one subplot
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
imagesc(orthogonal_slice)
colormap gray
xlabel('Z')
ylabel('X')
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title({'Normal View of Orthogonal Slice in XZ Plane';''})

% plot the first oblique slice
subplot(2,2,2)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(2) rotated_vox_dim(3) rotated_vox_dim(1)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);

% calculate and plot the second oblique slice
% define rotation angles for second oblique slice
alpha = 0; beta = 30; gamma = 30; 
% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);
% plot the second oblique slice
subplot(2,2,3)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(2) rotated_vox_dim(3) rotated_vox_dim(1)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);

% calculate and plot the third oblique slice
% define rotation angles for second oblique slice
alpha = 30; beta = 30; gamma = 30; 
% compute oblique slice and rotated voxel dimensions
[oblique_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    position, res_1, res_2, method, alpha, beta, gamma,'slice');
% strip low value border from oblique slice
oblique_slice = StripBorderObliqueSlice(oblique_slice, stencil_2D);
% plot the third oblique slice
subplot(2,2,4)
imagesc(oblique_slice)
daspect(1./[rotated_vox_dim(2) rotated_vox_dim(3) rotated_vox_dim(1)])
title({'Normal View of Oblique Slice';
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
% remove numbers from axes
set(gca,'XTick',[])
set(gca,'YTick',[])
colormap(gray);
suptitle({'Testing Oblique Slice Computation Using XZ Plane as Reference - High Resolution Sampling';
    'Cubic Spline Interpolation'; ['(Y Slice Position = ', num2str(position), ' mm)']})
% save image as .jpg file with appropriate name
saveas(gcf,'XZ_low_res','jpg');

catch
% Error message from LoadDICOMVOlume is shown if volume is not loaded
% correctly     
end 