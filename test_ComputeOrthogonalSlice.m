clc
close all 
clearvars

% MPHYGB24 - MATLAB coursework assignment 2017/18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Script for Task 2
% Tests function ComputeOrthogonalSlice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NB Script requires 'Bioinformatics Toolbox' for function suptitle.m,
% first used in line 111.

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

% Set orientation of view plane
orientation = 'X-Y';

% Set Z slice position (specified in mm):
position = 180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Nearest-Neighbour Interpolation
method = 'nearest';

% Case a)
% Both slice pixel dimensions are set to the smallest in-plane pixel
% dimension of the original image
res_1 = 0.9765625;
res_2 = 0.9765625; 

% compute orthogonal slice
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% display slice with scaled colours
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
% rotate 180 degrees to re-orientate the image correctly
slice_a = rot90(slice_a,2);
imagesc(slice_a);
% change colormap to greyscale
colormap(gray);
% label axes
xlabel('X / voxels')
ylabel('Y / voxels')
% ensure that slices are displayed with the correct scaling, according to
% voxel dimensions
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title({'a)';'';['Resolution (X,Y) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
% set axes limits
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

% Case b)
% Both slice pixel dimensions are set to the largest in-plane pixel
% dimension of the original image (this is the same as case a) here)
res_1 = 0.9765625;
res_2 = 0.9765625;

% compute orthogonal slice
slice_b = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice as above
subplot(1,3,2)
slice_b = rot90(slice_b,2);
imagesc(slice_b);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title({'b)';'';['Resolution (X,Y) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

% Case c)
% Both slice pixel dimensions are set to four times the largest in-plane
% pixel dimension of the original image
res_1 = 4*0.9765625;
res_2 = 4*0.9765625;

% compute orthogonal slice
slice_c = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice as above
subplot(1,3,3)
slice_c = rot90(slice_c,2);
imagesc(slice_c);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title({'c)';'';['Resolution (X,Y) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
% put a title above all subplots
suptitle({'XY Plane - Nearest-Neighbour Interpolation'; ['Z Slice Position = ', num2str(position), ' mm']})
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])
% save image as .jpg file with appropriate name
saveas(gcf,'XY_nearest_neighbour','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Linear Interpolation
method = 'linear';

% Case a)
% Both slice pixel dimensions are set to the smallest in-plane pixel
% dimension of the original image
res_1 = 0.9765625;
res_2 = 0.9765625; 

% compute orthogonal slice
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% display slice as above
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
slice_a = rot90(slice_a,2);
imagesc(slice_a);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title({'a)';'';['Resolution (X,Y) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

% Case b)
% Both slice pixel dimensions are set to the largest in-plane pixel
% dimension of the original image (this is the same as case a) here)
res_1 = 0.9765625;
res_2 = 0.9765625;

% compute orthogonal slice
slice_b = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice as above
subplot(1,3,2)
slice_b = rot90(slice_b,2);
imagesc(slice_b);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title({'b)';'';['Resolution (X,Y) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

% Case c)
% Both slice pixel dimensions are set to four times the largest in-plane
% pixel dimension of the original image
res_1 = 4*0.9765625;
res_2 = 4*0.9765625;

% compute orthogonal slice
slice_c = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice as above
subplot(1,3,3)
slice_c = rot90(slice_c,2);
imagesc(slice_c);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title({'c)';'';['Resolution (X,Y) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
suptitle({'XY Plane - Linear Interpolation'; ['Z Slice Position = ', num2str(position), ' mm']})
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])
saveas(gcf,'XY_linear','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Cubic spline Interpolation
method = 'spline';

% Case a)
% Both slice pixel dimensions are set to the smallest in-plane pixel
% dimension of the original image
res_1 = 0.9765625;
res_2 = 0.9765625; 

% compute orthogonal slice
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% display slice as above
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
slice_a = rot90(slice_a,2);
imagesc(slice_a);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title({'a)';'';['Resolution (X,Y) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

% Case b)
% Both slice pixel dimensions are set to the largest in-plane pixel
% dimension of the original image (this is the same as case a) here)
res_1 = 0.9765625;
res_2 = 0.9765625;

% compute orthogonal slice
slice_b = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice as above
subplot(1,3,2)
slice_b = rot90(slice_b,2);
imagesc(slice_b);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title({'b)';'';['Resolution (X,Y) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

% Case c)
% Both slice pixel dimensions are set to four times the largest in-plane
% pixel dimension of the original image
res_1 = 4*0.9765625;
res_2 = 4*0.9765625;

% compute orthogonal slice
slice_c = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice as above
subplot(1,3,3)
slice_c = rot90(slice_c,2);
imagesc(slice_c);
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title({'c)';'';['Resolution (X,Y) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
suptitle({'XY Plane - Cubic Spline Interpolation'; ['Z Slice Position = ', num2str(position), ' mm']})
saveas(gcf,'XY_cubic_spline','jpg');
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing YZ plane extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set orientation of view plane
orientation = 'Y-Z';

% Set X slice position (specified in mm):
position = 250;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Nearest-Neighbour Interpolation
method = 'nearest';

% Case a)
% Both slice pixel dimensions are set to the smallest in-plane pixel
% dimension of the original image
res_1 = 0.9765625;
res_2 = 0.9765625; 

% compute orthogonal slice
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
imagesc(fliplr(slice_a'));
% NB use transpose and fliplr above for correct image orientation
set(gca,'Xdir','reverse')
% Reverse x-axis values so that 0 is at front of patient
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title({'a)';'';['Resolution (Y,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])

% Case b)
% Both slice pixel dimensions are set to the largest in-plane pixel
% dimension of the original image
res_1 = 1;
res_2 = 1;

% compute orthogonal slice
slice_b = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice
subplot(1,3,2)
imagesc(fliplr(slice_b'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title({'b)';'';['Resolution (Y,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])

% Case c)
% Both slice pixel dimensions are set to four times the largest in-plane
% pixel dimension of the original image
res_1 = 4;
res_2 = 4;

% compute orthogonal slice
slice_c = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
subplot(1,3,3)
imagesc(fliplr(slice_c'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title({'c)';'';['Resolution (Y,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
suptitle({'YZ Plane - Nearest-Neighbour Interpolation'; ['X Slice Position = ', num2str(position), ' mm']})
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])
saveas(gcf,'YZ_nearest_neighbour','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Linear Interpolation
method = 'linear';

% Case a)
% Both slice pixel dimensions are set to the smallest in-plane pixel
% dimension of the original image
res_1 = 0.9765625;
res_2 = 0.9765625; 

% compute orthogonal slice
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
imagesc(fliplr(slice_a'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title({'a)';'';['Resolution (Y,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])

% Case b)
% Both slice pixel dimensions are set to the largest in-plane pixel
% dimension of the original image
res_1 = 1;
res_2 = 1;

% compute orthogonal slice
slice_b = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice
subplot(1,3,2)
imagesc(fliplr(slice_b'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title({'b)';'';['Resolution (Y,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])

% Case c)
% Both slice pixel dimensions are set to four times the largest in-plane
% pixel dimension of the original image
res_1 = 4;
res_2 = 4;

% compute orthogonal slice
slice_c = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
subplot(1,3,3)
imagesc(fliplr(slice_c'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title({'c)';'';['Resolution (Y,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
suptitle({'YZ Plane - Linear Interpolation'; ['X Slice Position = ', num2str(position), ' mm']})
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])
saveas(gcf,'YZ_linear','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Cubic spline Interpolation
method = 'spline';

% Case a)
% Both slice pixel dimensions are set to the smallest in-plane pixel
% dimension of the original image
res_1 = 0.9765625;
res_2 = 0.9765625; 

% compute orthogonal slice
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
imagesc(fliplr(slice_a'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title({'a)';'';['Resolution (Y,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])

% Case b)
% Both slice pixel dimensions are set to the largest in-plane pixel
% dimension of the original image
res_1 = 1;
res_2 = 1;

% compute orthogonal slice
slice_b = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice
subplot(1,3,2)
imagesc(fliplr(slice_b'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title({'b)';'';['Resolution (Y,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])

% Case c)
% Both slice pixel dimensions are set to four times the largest in-plane
% pixel dimension of the original image
res_1 = 4;
res_2 = 4;

% compute orthogonal slice
slice_c = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
subplot(1,3,3)
imagesc(fliplr(slice_c'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title({'c)';'';['Resolution (Y,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
suptitle({'YZ Plane - Cubic Spline Interpolation'; ['X Slice Position = ', num2str(position), ' mm']})
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])
saveas(gcf,'YZ_cubic_spline','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing XZ plane extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set orientation of view plane
orientation = 'X-Z';

% Set Y slice position (specified in mm):
position = 250;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Nearest-Neighbour Interpolation
method = 'nearest';

% Case a)
% Both slice pixel dimensions are set to the smallest in-plane pixel
% dimension of the original image
res_1 = 0.9765625;
res_2 = 0.9765625; 

% compute orthogonal slice
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
imagesc(fliplr(slice_a'));
% NB use transpose and fliplr above for correct image orientation
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title({'a)';'';['Resolution (X,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

% Case b)
% Both slice pixel dimensions are set to the largest in-plane pixel
% dimension of the original image
res_1 = 1;
res_2 = 1; 

% compute orthogonal slice
slice_b = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
subplot(1,3,2)
imagesc(fliplr(slice_b'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title({'b)';'';['Resolution (X,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

% Case c)
% Both slice pixel dimensions are set to four times the largest in-plane
% pixel dimension of the original image
res_1 = 4;
res_2 = 4; 

% compute orthogonal slice
slice_c = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
subplot(1,3,3)
imagesc(fliplr(slice_c'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title({'c)';'';['Resolution (X,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
suptitle({'XZ Plane - Nearest-Neighbour Interpolation'; ['Y Slice Position = ', num2str(position), ' mm']})
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])
saveas(gcf,'XZ_nearest_neighbour','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Linear Interpolation
method = 'linear';

% Case a)
% Both slice pixel dimensions are set to the smallest in-plane pixel
% dimension of the original image
res_1 = 0.9765625;
res_2 = 0.9765625; 

% compute orthogonal slice
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
imagesc(fliplr(slice_a'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title({'a)';'';['Resolution (X,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

% Case b)
% Both slice pixel dimensions are set to the largest in-plane pixel
% dimension of the original image
res_1 = 1;
res_2 = 1; 

% compute orthogonal slice
slice_b = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
subplot(1,3,2)
imagesc(fliplr(slice_b'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title({'b)';'';['Resolution (X,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

% Case c)
% Both slice pixel dimensions are set to four times the largest in-plane
% pixel dimension of the original image
res_1 = 4;
res_2 = 4; 

% compute orthogonal slice
slice_c = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
subplot(1,3,3)
imagesc(fliplr(slice_c'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title({'c)';'';['Resolution (X,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
suptitle({'XZ Plane - Linear Interpolation'; ['Y Slice Position = ', num2str(position), ' mm']})
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])
saveas(gcf,'XZ_linear','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Cubic spline Interpolation
method = 'spline';

% Case a)
% Both slice pixel dimensions are set to the smallest in-plane pixel
% dimension of the original image
res_1 = 0.9765625;
res_2 = 0.9765625; 

% compute orthogonal slice
slice_a = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
imagesc(fliplr(slice_a'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title({'a)';'';['Resolution (X,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

% Case b)
% Both slice pixel dimensions are set to the largest in-plane pixel
% dimension of the original image
res_1 = 1;
res_2 = 1; 

% compute orthogonal slice
slice_b = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
subplot(1,3,2)
imagesc(fliplr(slice_b'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title({'b)';'';['Resolution (X,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

% Case c)
% Both slice pixel dimensions are set to four times the largest in-plane
% pixel dimension of the original image
res_1 = 4;
res_2 = 4; 

% compute orthogonal slice
slice_c = ComputeOrthogonalSlice(Image, orientation, position, res_1, res_2, method);

% to display slice in correct orientation
subplot(1,3,3)
imagesc(fliplr(slice_c'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title({'c)';'';['Resolution (X,Z) = (', num2str(res_1), ' mm, ', num2str(res_2),' mm)']})
suptitle({'XZ Plane - Cubic Spline Interpolation'; ['Y Slice Position = ', num2str(position), ' mm']})
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])
saveas(gcf,'XZ_cubic_spline','jpg');

% Error message testing

% In XY plane
ComputeOrthogonalSlice(Image, 'X-Y', 180.1, res_1, res_2, method);
% In YZ plane
ComputeOrthogonalSlice(Image, 'Y-Z', -0.1, res_1, res_2, method);
% In XZ plane
ComputeOrthogonalSlice(Image, 'X-Z', 499.1, res_1, res_2, method);

catch
% Error message from LoadDICOMVOlume is shown if volume is not loaded
% correctly 
end

