clc
close all 
clearvars

% MPHYGB24 - MATLAB coursework assignment 2017/18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Script for Task 8
% Tests functions ComputeObliqueSlab and DisplaySlab3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load all slices of 3D DICOM volume using function from Task 1:
z_limits(1) = 1;
z_limits(2) = 181;

try
    
Image = LoadDICOMVolume(z_limits);

% To extract image dimensions [no_rows, no_cols, no_slices]
image_dim = size(Image.ImageData);
% To extract voxel dimensions in mm (used to scale axes in plots)
vox_dim = Image.VoxelDimensions; % [dy dx dz]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Testing XY plane extraction, Linear Interpolation and High Resolution
% Sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set orientation of plane & interpolation technique
orientation = 'X-Y';
method = 'linear';

% Set Z slice position of central slice (specified in mm):
central_position = 90;

% Thickness of slab (specified in mm):
slab_thickness_mm = 50;

% Set slab sampling resolutions to be equal to voxel dimensions in
% respective directions
res_1 = vox_dim(2); % x resolution
res_2 = vox_dim(1); % y resolution
res_3 = vox_dim(3); % z resolution

% define rotation angles for central oblique slice
alpha = -45; % around x-axis
beta = 0; % around y-axis
gamma = 0; % around z-axis 

% compute slab, centre index and no slices in slab
[slab, centre_index, no_slices] = ComputeObliqueSlab(Image, orientation, ...
    central_position, res_1, res_2, res_3, method, alpha, beta, gamma, slab_thickness_mm);

% show slab in new 3D wire frame plot
DisplaySlab3D(Image, slab, orientation, central_position, alpha, ...
    beta, gamma, centre_index, no_slices,'new')

% insert title containing all relevant information into current figure
title({'Oblique Slab - XY Plane - Linear Interpolation';
    ['Z Position of Slab Centre = ', num2str(central_position),...
    ' mm - Slab Thickness = ' , num2str(slab_thickness_mm),' mm'];... 
    ['Slab Resolution (x, y, z) = (', num2str(res_1), ' mm, ',num2str(res_2),...
    ' mm, ', num2str(res_3), ' mm)']; ...
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
saveas(gcf,'XY_Test_Slab.jpg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Testing YZ plane extraction, Nearest Neighbour Interpolation and Low
% Resolution Sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set orientation of plane & interpolation technique
orientation = 'Y-Z';
method = 'linear';

% Set X slice position of central slice (specified in mm):
central_position = 250;

% Thickness of slab (specified in mm):
slab_thickness_mm = 20;

% Set slab sampling resolutions all to be four times the largest voxel
% dimension
res_1 = 4; % y resolution
res_2 = 4; % z resolution
res_3 = 4; % x resolution

% define rotation angles for central oblique slice
alpha = 0; % around x-axis
beta = 30; % around y-axis
gamma = 0; % around z-axis 

% compute slab, centre index and no slices in slab
[slab, centre_index, no_slices] = ComputeObliqueSlab(Image, orientation, ...
    central_position, res_1, res_2, res_3, method, alpha, beta, gamma, slab_thickness_mm);

% show slab in new 3D wire frame plot
DisplaySlab3D(Image, slab, orientation, central_position, alpha, ...
    beta, gamma, centre_index, no_slices, 'new')

% insert title containing all relevant information into current figure
title({'Oblique Slab - YZ Plane - Nearest Neighbour Interpolation';
    ['X Position of Slab Centre = ', num2str(central_position),...
    ' mm - Slab Thickness = ' , num2str(slab_thickness_mm),' mm'];... 
    ['Slab Resolution (x, y, z) = (', num2str(res_3), ' mm, ',num2str(res_1),...
    ' mm, ', num2str(res_2), ' mm)']; ...
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
saveas(gcf,'YZ_Test_Slab.jpg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Testing XZ plane extraction, Cubic Spline Interpolation and High
% Resolution Sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set orientation of plane & interpolation technique
orientation = 'X-Z';
method = 'spline';

% Set X slice position of central slice (specified in mm):
central_position = 255;

% Thickness of slab (specified in mm):
slab_thickness_mm = 50;

% Set slab sampling resolutions all to be respective pixel dimensions
res_1 = vox_dim(2); % x resolution
res_2 = vox_dim(3); % z resolution
res_3 = vox_dim(1); % y resolution

% define rotation angles for central oblique slice
alpha = 0; % around x-axis
beta = 0; % around y-axis
gamma = 20; % around z-axis 

% compute slab, centre index and no slices in slab
[slab, centre_index, no_slices] = ComputeObliqueSlab(Image, orientation, ...
    central_position, res_1, res_2, res_3, method, alpha, beta, gamma, slab_thickness_mm);

% show slab in new 3D wire frame plot
DisplaySlab3D(Image, slab, orientation, central_position, alpha, ...
    beta, gamma, centre_index, no_slices, 'new')

% insert title containing all relevant information into current figure
title({'Oblique Slab - XZ Plane - Cubic Spline Interpolation';
    ['Y Position of Slab Centre = ', num2str(central_position),...
    ' mm - Slab Thickness = ' , num2str(slab_thickness_mm),' mm'];... 
    ['Slab Resolution (x, y, z) = (', num2str(res_1), ' mm, ',num2str(res_3),...
    ' mm, ', num2str(res_2), ' mm)']; ...
    ['R_x = ', num2str(alpha), '^o - R_y = ', num2str(beta), ...
    '^o - R_z = ', num2str(gamma), '^o']})
saveas(gcf,'XZ_Test_Slab.jpg')

catch 
% Error message from LoadDICOMVOlume is shown if volume is not loaded
% correctly       
end