clc
close all 
clearvars

% MPHYGB24 - MATLAB coursework assignment 2017/18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Script for Task 9
% Tests functions ComputeProjectionImage and DisplayImage2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NB Script requires 'Bioinformatics Toolbox' for function suptitle.m,
% first used in line 98.

% Load all slices of 3D DICOM volume using function from Task 1:
z_limits(1) = 1;
z_limits(2) = 181;
    
try

Image = LoadDICOMVolume(z_limits);

% To extract image dimensions [no_rows, no_cols, no_slices]
image_dim = size(Image.ImageData);
% calculate minimum value in vol
low_value = min(min(min(Image.ImageData)));
% To extract voxel dimensions in mm (used to scale axes in plots)
vox_dim = Image.VoxelDimensions; % [dy dx dz]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Testing XY plane (same slice thickness, position and orientation as
% figure 5 of question sheet), high resolution sampling and linear interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set orientation of plane & interpolation technique
orientation = 'X-Y';
method = 'linear';

% Set Z slice position of central slice (specified in mm):
central_position = 9;

% Thickness of slab (specified in mm):
slab_thickness_mm = 20;  

% Set slab sampling resolutions to be equal to voxel dimensions in
% respective directions
res_1 = vox_dim(2); % x resolution
res_2 = vox_dim(1); % y resolution
res_3 = vox_dim(3); % z resolution

% define rotation angles for central oblique slice
alpha = 0; % around x-axis
beta = 0; % around y-axis
gamma = 0; % around z-axis 

% compute slab, centre index, no slices in slab, rotated voxel dimensions
% and 3D stencil corresponding to slab
[slab, centre_index, no_slices, rotated_vox_dim, stencil_3D] = ComputeObliqueSlab(Image, ...
    orientation, central_position, res_1, res_2, res_3, method, alpha, beta, ...
    gamma, slab_thickness_mm);

% show central slice
view = 'central';
% compute central slice projection image
projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices);
% open new figure in maximised window
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
% display image with title
DisplayImage2D(projection_image, rotated_vox_dim, orientation)
title('Central Slice Image')

% compute maximum intensity projection image as above
view = 'max';
projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices);
subplot(2,2,2)
DisplayImage2D(projection_image, rotated_vox_dim, orientation)
title('Maximum Intensity Projection Image')

% compute minimum intensity projection as above
subplot(2,2,3)
view = 'min';
projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices);
DisplayImage2D(projection_image, rotated_vox_dim, orientation)
title('Minimum Intensity Projection Image')

% compute median intensity projection as above
subplot(2,2,4)
view = 'median';
projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices);
DisplayImage2D(projection_image, rotated_vox_dim, orientation)
title('Median Intensity Projection Image')
 
% insert title containing all relevant information into current figure
suptitle({['Various Projection Images - XY Plane - Linear Interpolation - Z Position of Slab Centre = ', ...
    num2str(central_position),...
    ' mm - Slab Thickness = ' , num2str(slab_thickness_mm),' mm'];... 
    ['Slab Resolution (x, y, z) = (', num2str(res_1), ' mm, ',num2str(res_2),...
    ' mm, ', num2str(res_3), ' mm)']; ...
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
saveas(gcf,'XY_Test_Project.jpg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Testing YZ plane, nearest neighbour interpolation, low resolution
% sampling and rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set orientation of plane & interpolation technique
orientation = 'Y-Z';
method = 'nearest';

% Set X slice position of central slice (specified in mm):
central_position = 250;

% Thickness of slab (specified in mm):
slab_thickness_mm = 20;  

% Set slab sampling resolutions to be equal to largest voxel dimension x 4
res_1 = 4*vox_dim(3); % y resolution
res_2 = 4*vox_dim(3); % z resolution
res_3 = 4*vox_dim(3); % x resolution

% define rotation angles for central oblique slice
alpha = 0; % around x-axis
beta = 0; % around y-axis
gamma = 45; % around z-axis 

% compute slab, centre index, no slices in slab, rotated voxel dimensions
% and 3D stencil corresponding to slab
[slab, centre_index, no_slices, rotated_vox_dim, stencil_3D] = ComputeObliqueSlab(Image, ...
    orientation, central_position, res_1, res_2, res_3, method, alpha, beta, ...
    gamma, slab_thickness_mm);

% show central slice
view = 'central';
% compute central slice projection image
projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices);
% open new figure in maximised window
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
% display image with title
DisplayImage2D(projection_image, rotated_vox_dim, orientation)
title('Central Slice Image')

% compute maximum intensity projection image as above
view = 'max';
projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices);
subplot(2,2,2)
DisplayImage2D(projection_image, rotated_vox_dim, orientation)
title('Maximum Intensity Projection Image')

% compute minimum intensity projection as above
subplot(2,2,3)
view = 'min';
projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices);
DisplayImage2D(projection_image, rotated_vox_dim, orientation)
title('Minimum Intensity Projection Image')

% compute median intensity projection as above
subplot(2,2,4)
view = 'median';
projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices);
DisplayImage2D(projection_image, rotated_vox_dim, orientation)
title('Median Intensity Projection Image')
 
% insert title containing all relevant information into current figure
suptitle({['Various Projection Images - YZ Plane - Nearest Neighbour Interpolation - X Position of Slab Centre = ', ...
    num2str(central_position),...
    ' mm - Slab Thickness = ' , num2str(slab_thickness_mm),' mm'];... 
    ['Slab Resolution (x, y, z) = (', num2str(res_3), ' mm, ',num2str(res_2),...
    ' mm, ', num2str(res_3), ' mm)']; ...
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
saveas(gcf,'YZ_Test_Project.jpg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Testing XZ plane, cubic spline interpolation, high resolution sampling
% and rotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set orientation of plane & interpolation technique
orientation = 'X-Z';
method = 'spline';

% Set X slice position of central slice (specified in mm):
central_position = 250;

% Thickness of slab (specified in mm):
slab_thickness_mm = 20;

% Set slab sampling resolutions to be equal to respective voxel dimensions
res_1 = vox_dim(1); % x resolution
res_2 = vox_dim(3); % z resolution
res_3 = vox_dim(2); % y resolution


% define rotation angles for central oblique slice
alpha = 0; % around x-axis
beta = 0; % around y-axis
gamma = 45; % around z-axis 

% compute slab, centre index, no slices in slab, rotated voxel dimensions
% and 3D stencil corresponding to slab
[slab, centre_index, no_slices, rotated_vox_dim, stencil_3D] = ComputeObliqueSlab(Image, ...
    orientation, central_position, res_1, res_2, res_3, method, alpha, beta, ...
    gamma, slab_thickness_mm);

% show central slice
view = 'central';
% compute central slice projection image
projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices);
% open new figure in maximised window
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
% display image with title
DisplayImage2D(projection_image, rotated_vox_dim, orientation)
title('Central Slice Image')

% compute maximum intensity projection image as above
view = 'max';
projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices);
subplot(2,2,2)
DisplayImage2D(projection_image, rotated_vox_dim, orientation)
title('Maximum Intensity Projection Image')

% compute minimum intensity projection as above
subplot(2,2,3)
view = 'min';
projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices);
DisplayImage2D(projection_image, rotated_vox_dim, orientation)
title('Minimum Intensity Projection Image')

% compute median intensity projection as above
subplot(2,2,4)
view = 'median';
projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices);
DisplayImage2D(projection_image, rotated_vox_dim, orientation)
title('Median Intensity Projection Image')
 
% insert title containing all relevant information into current figure
suptitle({['Various Projection Images - XZ Plane - Cubic Spline Interpolation - Y Position of Slab Centre = ', ...
    num2str(central_position),...
    ' mm - Slab Thickness = ' , num2str(slab_thickness_mm),' mm'];... 
    ['Slab Resolution (x, y, z) = (', num2str(res_1), ' mm, ',num2str(res_3),...
    ' mm, ', num2str(res_2), ' mm)']; ...
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
saveas(gcf,'XZ_Test_Project.jpg')

catch
% Error message from LoadDICOMVOlume is shown if volume is not loaded
% correctly  
end