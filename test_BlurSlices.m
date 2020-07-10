clc
close all 
clearvars

% MPHYGB24 - MATLAB coursework assignment 2017/18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Script for Task 4
% Tests function BlurSlices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NB Script requires 'Bioinformatics Toolbox' for function suptitle.m,
% first used in line 103.

% Load all slices of 3D DICOM volume using function from Task 1:
z_limits(1) = 1;
z_limits(2) = 181;

try

Image = LoadDICOMVolume(z_limits);

% To assign 3D image intensities to a matrix variable
vol = Image.ImageData;

% To extract voxel dimensions in mm (used to scale axes in plots)
vox_dim = Image.VoxelDimensions; % [dy dx dz]

% To extract image dimensions in voxels
image_dim = size(Image.ImageData); % [rows cols slices]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing XY plane blurring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) To test 2D symmetrical image blurring and ability of function to handle more
% than one slice

% Take 2 slices far apart on XY plane
slab = vol(:,:,[40 160]);
% assign sigmax and sigmay values
sigmax = 5; 
sigmay = 5; 

% blur the two slices in 2D in one call to BlurSlices
slab_blurred = BlurSlices(slab,'X-Y',sigmax,sigmay,'','both');

% Plot slice 40 before and after blurring
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
% rotate 180 degrees to re-orientate the image correctly
slice_a = squeeze(rot90(slab(:,:,1),2));
% display slice with scaled colours
imagesc(slice_a)
% change colormap to greyscale
colormap(gray);
% label axes
xlabel('X / voxels')
ylabel('Y / voxels')
% ensure that slices are displayed with the correct scaling, according to
% voxel dimensions
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title('Slice 40 Before Blurring')
% set axes limits
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

subplot(2,2,2)
slice_b = squeeze(rot90(slab_blurred(:,:,1),2));
imagesc(slice_b)
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title(['Slice 40 After 2D Blurring - \sigma_x = ',num2str(sigmax) ...
    ,' - \sigma_y = ',num2str(sigmay)])
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

% Plot slice 160 before and after blurring
subplot(2,2,3)
slice_c = squeeze(rot90(slab(:,:,2),2));
imagesc(slice_c)
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title('Slice 160 Before Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

subplot(2,2,4)
slice_b = squeeze(rot90(slab_blurred(:,:,2),2));
imagesc(slice_b)
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title(['Slice 160 After 2D Blurring - \sigma_x = ',num2str(sigmax) ...
    ,' - \sigma_y = ',num2str(sigmay)])
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])
suptitle('Testing Symmetrical Gaussian Blurring on Multiple Slices in X-Y Plane')
saveas(gcf,'XY_multiple_blur','jpg');

% 2) Testing ability of BlurSlices to handle just one slice, and checking
% sigmax and sigmay correspond to the correct axes, and that 1D blurring
% works

% form a single slice slab
slab = vol(:,:,90);
% assign sigmay value
sigmay = 20; 
% blur the single slice slab in the Y-direction only
slab_blurred_a = BlurSlices(slab,'X-Y','',sigmay,'','Y');
% assign sigmax value
sigmax = 20; 
% blur the single slice slab in the X-direction only
slab_blurred_b = BlurSlices(slab,'X-Y',sigmax,'','','X');
% blur the single slice slab in both X and Y directions
slab_blurred_c = BlurSlices(slab,'X-Y',sigmax,sigmay,'','both');

% Plot slice 90 before and after three different blurring permutations 
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
slice_a = squeeze(rot90(slab,2));
imagesc(slice_a)
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title('Slice 90 Before Blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

subplot(2,2,2)
slice_b = squeeze(rot90(slab_blurred_a,2));
imagesc(slice_b)
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title('Slice 90 After 1D Blurring - \sigma_y = 20')
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

subplot(2,2,3)
slice_c = squeeze(rot90(slab_blurred_b,2));
imagesc(slice_c)
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title('Slice 90 After 1D Blurring - \sigma_x = 20')
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])

subplot(2,2,4)
slice_d = squeeze(rot90(slab_blurred_c,2));
imagesc(slice_d)
colormap(gray);
xlabel('X / voxels')
ylabel('Y / voxels')
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])
title('Slice 90 After 2D Blurring - \sigma_x = 20 - \sigma_y = 20')
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])
suptitle('Testing Asymmetrical Gaussian Blurring on a Single Slice in X-Y Plane')
saveas(gcf,'XY_single_blur','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing YZ plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3) To test symmetrical image blurring and ability of function to handle more
% than one slice

% Take 2 slices far apart on YZ plane
slab = vol(:,[100 250],:);
% assign sigmay and sigmaz values
sigmay = 5; 
sigmaz = 5; 

% blur the two slices in 2D in one call to BlurSlices
slab_blurred = BlurSlices(slab,'Y-Z','',sigmay,sigmaz,'both');

% Plot slice 100 before and after blurring
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
slice_a = squeeze(slab(:,1,:));
imagesc(fliplr(slice_a'));
% NB use transpose and fliplr above for correct image orientation
set(gca,'Xdir','reverse')
% Reverse x-axis values so that 0 is at front of patient
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title('Slice 100 Before Blurring')
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])

subplot(2,2,2)
slice_b = squeeze(slab_blurred(:,1,:));
imagesc(fliplr(slice_b'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title(['Slice 100 After 2D Blurring - \sigma_y = ',num2str(sigmay) ...
    ,' - \sigma_z = ',num2str(sigmaz)])
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])

% Plot slice 250 before and after blurring
subplot(2,2,3)
slice_c = squeeze(slab(:,2,:));
imagesc(fliplr(slice_c'));
% NB use transpose and fliplr above for correct image orientation
set(gca,'Xdir','reverse')
% Reverse x-axis values so that 0 is at front of patient
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title('Slice 250 Before Blurring')
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])

subplot(2,2,4)
slice_d = squeeze(slab_blurred(:,2,:));
imagesc(fliplr(slice_d'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title(['Slice 250 After 2D Blurring - \sigma_y = ',num2str(sigmay) ...
    ,' - \sigma_z = ',num2str(sigmaz)])
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])
suptitle('Testing Symmetrical Gaussian Blurring on Multiple Slices in Y-Z Plane')
saveas(gcf,'YZ_single_blur','jpg');

% 4) Testing ability of BlurSlices to handle just one slice, and checking
% sigmay and sigmaz correspond to the correct axes, and that 1D blurring
% works

% form a single slice slab
slab = vol(:,250,:);
% assign sigmaz value
sigmaz = 20; 
% blur the single slice slab in Z-direction only
slab_blurred_a = BlurSlices(slab,'Y-Z','','',sigmaz,'Z');
% assign sigmay value
sigmay = 20; 
% blur the single slice slab in Y-direction only
slab_blurred_b = BlurSlices(slab,'Y-Z','',sigmay,'','Y');
% blur the single slice slab in both directions
slab_blurred_c = BlurSlices(slab,'Y-Z','',sigmay,sigmaz,'both');

% Plot slice 250 before and after three different blurring permutations
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
slice_a = squeeze(slab);
imagesc(fliplr(slice_a'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title('Slice 250 Before Blurring')
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])

subplot(2,2,2)
slice_b = squeeze(slab_blurred_a);
imagesc(fliplr(slice_b'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title('Slice 250 After 1D Blurring - \sigma_z = 20')
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])

subplot(2,2,3)
slice_c = squeeze(slab_blurred_b);
imagesc(fliplr(slice_c'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title('Slice 250 After 1D Blurring - \sigma_y = 20')
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])

subplot(2,2,4)
slice_d = squeeze(slab_blurred_c);
imagesc(fliplr(slice_d'));
set(gca,'Xdir','reverse')
colormap(gray);
xlabel('Y / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(1) vox_dim(3) vox_dim(2)])
title('Slice 250 After 2D Blurring - \sigma_y = 20 - \sigma_z = 20')
xlim([1,image_dim(1)])
ylim([1,image_dim(3)])
suptitle('Testing Asymmetrical Gaussian Blurring on a Single Slice in Y-Z Plane')
saveas(gcf,'YZ_multiple_blur','jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing XZ plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5) To test symmetrical image blurring and ability of function to handle more
% than one slice

% Take 2 slices from XZ plane
slab = vol([200 300],:,:);
% assign sigmax and sigmaz values
sigmax = 5;
sigmaz = 5; 

% blur the two slices in 2D in one call to BlurSlices
slab_blurred = BlurSlices(slab,'X-Z',sigmax,'',sigmaz,'both');

% Plot slice 200 before and after blurring
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
slice_a = squeeze(slab(1,:,:));
imagesc(fliplr(slice_a'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title('Slice 200 before blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

subplot(2,2,2)
slice_b = squeeze(slab_blurred(1,:,:));
imagesc(fliplr(slice_b'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title(['Slice 40 After 2D Blurring - \sigma_x = ',num2str(sigmax) ...
    ,' - \sigma_z = ',num2str(sigmaz)])
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

% Plot slice 300 before and after blurring
subplot(2,2,3)
slice_c = squeeze(slab(2,:,:));
imagesc(fliplr(slice_c'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title('Slice 300 before blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

subplot(2,2,4)
slice_d = squeeze(slab_blurred(2,:,:));
imagesc(fliplr(slice_d'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title(['Slice 300 After 2D Blurring - \sigma_x = ',num2str(sigmax) ...
    ,' - \sigma_z = ',num2str(sigmaz)])
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])
suptitle('Testing Symmetrical Gaussian Blurring on Multiple Slices in X-Z Plane')
saveas(gcf,'XZ_multiple_blur','jpg');

% 6) Testing ability of BlurSlices to handle just one slice, and checking
% sigmax and sigmaz correspond to the correct axes, and that 1D blurring
% works

% form a single slice slab
slab = vol(250,:,:);
% assign sigmaz value
sigmaz = 20; 
% blur the single slice slab in Z-direction only
slab_blurred_a = BlurSlices(slab,'X-Z','','',sigmaz,'Z');
% assign sigmax value
sigmax = 20; 
% blur the single slice slab in X-direction only
slab_blurred_b = BlurSlices(slab,'X-Z',sigmax,'','','X');
% blur the single slice slab in both directions
slab_blurred_c = BlurSlices(slab,'X-Z',sigmax,'',sigmaz,'both');

% Plot slice 250 before and after three different blurring permutations
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
slice_a = squeeze(slab);
imagesc(fliplr(slice_a'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title('Slice 250 before blurring')
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

subplot(2,2,2)
slice_b = squeeze(slab_blurred_a);
imagesc(fliplr(slice_b'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title('Slice 250 After 1D Blurring - \sigma_z = 20')
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

subplot(2,2,3)
slice_c = squeeze(slab_blurred_b);
imagesc(fliplr(slice_c'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title('Slice 250 After 1D Blurring - \sigma_x = 20')
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])

subplot(2,2,4)
slice_d = squeeze(slab_blurred_c);
imagesc(fliplr(slice_d'));
colormap(gray);
xlabel('X / voxels')
ylabel('Z / voxels')
daspect(1./[vox_dim(2) vox_dim(3) vox_dim(1)])
title('Slice 250 After 2D Blurring - \sigma_x = 20 - \sigma_z = 20')
xlim([1,image_dim(2)])
ylim([1,image_dim(3)])
suptitle('Testing Asymmetrical Gaussian Blurring on a Single Slice in X-Z Plane')
saveas(gcf,'XZ_single_blur','jpg');

catch
% Error message from LoadDICOMVOlume is shown if volume is not loaded
% correctly 
end