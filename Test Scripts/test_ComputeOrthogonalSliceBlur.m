clc
close all 
clearvars

% MPHYGB24 - MATLAB coursework assignment 2017/18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Script for Task 5
% Tests function ComputeOrthogonalSliceBlur when the resolution of the
% displayed slice is relatively low (repeat of task 2c). Compares results
% directly with ComputeOrthogonalSlice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NB Script requires 'Bioinformatics Toolbox' for function suptitle.m,
% first used in line 81.

% Load all slices of 3D DICOM volume using function from Task 1:
z_limits(1) = 1;
z_limits(2) = 181;

try 
    
Image = LoadDICOMVolume(z_limits);

% To extract voxel dimensions in mm (used to scale axes in plots)
vox_dim = Image.VoxelDimensions; % [dy dx dz]

% To extract image dimensions in voxels
image_dim = size(Image.ImageData); % [rows cols slices]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing XY plane extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set orientation of plane
orientation = 'X-Y';

% Set Z slice position (specified in mm):
position = 90;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Nearest-Neighbour Interpolation
method = 'nearest';

% Both slice pixel dimensions are set to four times the largest in-plane
% pixel dimension of the original image
res_1 = 4*0.9765625;
res_2 = 4*0.9765625;

% compute orthogonal slice without blurring
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);
% compute orthogonal slice with blurring
slice_b = ComputeOrthogonalSliceBlur(Image, orientation, position, res_1, res_2, method);

% to display slice_a
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
slice_a = rot90(slice_a,2);
imagesc(slice_a);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title('Without Pre-Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

% to display slice_b
subplot(1,2,2)
slice_b = rot90(slice_b,2);
imagesc(slice_b);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title('With Pre-Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

% put a title above all subplots
suptitle({'XY Plane - Nearest-Neighbour Interpolation'; ...
    ['Resolution (X,Y) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']; ...
    ['Z Slice Position = ', num2str(position), ' mm']})
saveas(gcf,'XY_NN_blur','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Linear Interpolation
method = 'linear';

% compute orthogonal slice without blurring
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);
% compute orthogonal slice with blurring
slice_b = ComputeOrthogonalSliceBlur(Image, orientation, position, res_1, res_2, method);

% to display slice_a
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
slice_a = rot90(slice_a,2);
imagesc(slice_a);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title('Without Pre-Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

% to display slice_b
subplot(1,2,2)
slice_b = rot90(slice_b,2);
imagesc(slice_b);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title('With Pre-Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

% put a title above all subplots
suptitle({'XY Plane - Linear Interpolation'; ...
    ['Resolution (X,Y) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']; ...
    ['Z Slice Position = ', num2str(position), ' mm']})
saveas(gcf,'XY_lin_blur','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Cubic spline Interpolation
method = 'spline';

% compute orthogonal slice without blurring
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);
% compute orthogonal slice with blurring
slice_b = ComputeOrthogonalSliceBlur(Image, orientation, position, res_1, res_2, method);

% to display slice_a
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
slice_a = rot90(slice_a,2);
imagesc(slice_a);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title('Without Pre-Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

% to display slice_b
subplot(1,2,2)
slice_b = rot90(slice_b,2);
imagesc(slice_b);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title('With Pre-Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

% put a title above all subplots
suptitle({'XY Plane - Cubic Spline Interpolation'; ...
    ['Resolution (X,Y) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']; ...
    ['Z Slice Position = ', num2str(position), ' mm']})
saveas(gcf,'XY_cubic_blur','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing YZ plane extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set orientation of plane
orientation = 'Y-Z';

% Set X slice position (specified in mm):
position = 250;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Nearest-Neighbour Interpolation
method = 'nearest';

% Both slice pixel dimensions are set to four times the largest in-plane
% pixel dimension of the original image
res_1 = 4;
res_2 = 4;

% compute orthogonal slice without pre-blurring
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);
% compute orthogonal slice with pre-blurring
slice_b = ComputeOrthogonalSliceBlur(Image, orientation, position, res_1, res_2, method);

% to display slice_a in correct orientation
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
imagesc(fliplr(slice_a'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])
title('Without Pre-Blurring')

% to display slice_b in correct orientation
subplot(1,2,2)
imagesc(fliplr(slice_b'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])
title('With Pre-Blurring')

% put a title above all subplots
suptitle({'YZ Plane - Nearest-Neighbour Interpolation'; ...
    ['Resolution (Y,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']; ...
    ['X Slice Position = ', num2str(position), ' mm']})
saveas(gcf,'YZ_NN_blur','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Linear Interpolation
method = 'linear';

% compute orthogonal slice without pre-blurring
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);
% compute orthogonal slice with pre-blurring
slice_b = ComputeOrthogonalSliceBlur(Image, orientation, position, res_1, res_2, method);

% to display slice_a in correct orientation
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
imagesc(fliplr(slice_a'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])
title('Without Pre-Blurring')

% to display slice_b in correct orientation
subplot(1,2,2)
imagesc(fliplr(slice_b'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])
title('With Pre-Blurring')

% put a title above all subplots
suptitle({'YZ Plane - Linear Interpolation'; ...
    ['Resolution (Y,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']; ...
    ['X Slice Position = ', num2str(position), ' mm']})
saveas(gcf,'YZ_lin_blur','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Cubic spline Interpolation
method = 'spline';

% compute orthogonal slice without pre-blurring
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);
% compute orthogonal slice with pre-blurring
slice_b = ComputeOrthogonalSliceBlur(Image, orientation, position, res_1, res_2, method);

% to display slice_a in correct orientation
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
imagesc(fliplr(slice_a'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])
title('Without Pre-Blurring')

% to display slice_b in correct orientation
subplot(1,2,2)
imagesc(fliplr(slice_b'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])
title('With Pre-Blurring')

% put a title above all subplots
suptitle({'YZ Plane - Cubic Spline Interpolation'; ...
    ['Resolution (Y,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']; ...
    ['X Slice Position = ', num2str(position), ' mm']})
saveas(gcf,'YZ_cubic_blur','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing XZ plane extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set orientation of plane
orientation = 'X-Z';

% Set Y slice position (specified in mm):
position = 250;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Nearest-Neighbour Interpolation
method = 'nearest';

% Both slice pixel dimensions are set to four times the largest in-plane
% pixel dimension of the original image
res_1 = 4;
res_2 = 4; 

% compute orthogonal slice without pre-blurring
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);
% compute orthogonal slice with pre-blurring
slice_b = ComputeOrthogonalSliceBlur(Image, orientation, position, res_1, res_2, method);

% to display slice_a in correct orientation
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
imagesc(fliplr(slice_a'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title('Without Pre-Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

% to display slice_b in correct orientation
subplot(1,2,2)
imagesc(fliplr(slice_b'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title('With Pre-Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

% put a title above all subplots
suptitle({'XZ Plane - Nearest Neighbour Interpolation'; ...
    ['Resolution (X,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']; ...
    ['Y Slice Position = ', num2str(position), ' mm']})
saveas(gcf,'XZ_NN_blur','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Linear Interpolation
method = 'linear';

% compute orthogonal slice without pre-blurring
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);
% compute orthogonal slice with pre-blurring
slice_b = ComputeOrthogonalSliceBlur(Image, orientation, position, res_1, res_2, method);

% to display slice_a in correct orientation
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
imagesc(fliplr(slice_a'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title('Without Pre-Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

% to display slice_b in correct orientation
subplot(1,2,2)
imagesc(fliplr(slice_b'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title('With Pre-Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

% put a title above all subplots
suptitle({'XZ Plane - Linear Interpolation'; ...
    ['Resolution (X,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']; ...
    ['Y Slice Position = ', num2str(position), ' mm']})
saveas(gcf,'XZ_lin_blur','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Cubic spline Interpolation
method = 'spline';

% compute orthogonal slice without pre-blurring
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);
% compute orthogonal slice with pre-blurring
slice_b = ComputeOrthogonalSliceBlur(Image, orientation, position, res_1, res_2, method);

% to display slice_a in correct orientation
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
imagesc(fliplr(slice_a'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title('Without Pre-Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

% to display slice_b in correct orientation
subplot(1,2,2)
imagesc(fliplr(slice_b'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title('With Pre-Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

% put a title above all subplots
suptitle({'XZ Plane - Cubic Spline Interpolation'; ...
    ['Resolution (X,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']; ...
    ['Y Slice Position = ', num2str(position), ' mm']})
saveas(gcf,'XZ_cubic_blur','jpg');

catch
% Error message from LoadDICOMVOlume is shown if volume is not loaded
% correctly 
end