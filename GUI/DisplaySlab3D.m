function DisplaySlab3D(Image, slab, orientation, central_position, alpha, ...
    beta, gamma, centre_index, no_slices, image, destination) %#ok<INUSD>
% DisplaySlab3D displays a 3D slab in a 3D view of a wire frame plot, to be
% used in right hand view of SliceViewer GUI. 
%
% DESCRIPTION: DisplaySlab3D(Image, slab, orientation, central_position, alpha, ...
%    beta, gamma, centre_index, no_slices)
%
% INPUTS:
%       Image (1 x 1 structure with two fields) - 
%           .ImageData (double matrix) - of dimensions (number of rows, 
%           number of columns,number of slices) containing the voxel grey 
%           level values
%           .VoxelDimensions (double vector) - a 1 by 3 vector containing
%           the (y,x,z) voxel dimensions in mm, respectively
%       
%       orientation (character string) - determines reference slice plane orientation
%           'X-Y' - XY plane, orthogonal to Z axis
%           'Y-Z' - YZ plane, orthogonal to X axis
%           'X-Z' - XZ plane, orthogonal to Y axis
%
%       central_position (double scalar) - determines position of reference 
%       slice in mm along axis orthogonal to the slice plane:
%           For the Z slice position, zero is toward the head of the 
%               patient.
%           For the Y slice position, zero is toward the front of the 
%               patient. 
%           For the X slice position, zero is toward the left of the
%               patient (for this volume this end of the voxel range is
%               where the arm is shown).
%
%       alpha (double scalar) - rotation angle about x axis, alpha can have
%       the following range of values:
%           -180 degrees <= alpha <= +180 degrees.
%
%       beta (double scalar) - rotation angle about y axis, beta can have
%       the following range of values:
%           -180 degrees <= beta <= +180 degrees.
%
%       gamma (double scalar) - rotation angle about z axis, gamma can have
%       the following range of values:
%           -180 degrees <= gamma <= +180 degrees.
%
%       centre_index (integer scalar) - the central index of the number of
%           slices within the returned slab. This is 1 if the slab 
%
%       no_slices (integer scalar) - the number of slices in the slab
%
%       image (character string) - if set to 'new' then figure opens in
%       a maximised window. Set to '' for the figure to plot by default.
%
%       destinaton (axes handle) - handle to the axes on which the projection 
%       image should be shown. This is an optional input argument. 
%
%       NB all three rotation angles are counterclockwise in a RHS
%       cartesian coordinate system, or clockwise in a LHS cartesian
%       coordinate system.
%
%OUTPUTS:
%       NONE.
%
% FUNCTION DEPENDENCIES:
%       NONE.
%
% AUTHOR:
%       Anonymised for MPHYGB24 MATLAB coursework assignment 2017/18

% To assign 3D image intensities to a matrix variable
vol = Image.ImageData;
image_dim = size(vol);
% To extract voxel dimensions in mm (used to scale axes in plots)
vox_dim = Image.VoxelDimensions; % [dy dx dz]

% Assign a 3D matrix of zeros of the same dimensions as the whole image
% volume
V = zeros(image_dim(2),image_dim(1),image_dim(3));

% To carry out orientation specific tasks:

% For XY slice
if strcmp(orientation,'X-Y') == 1
    
    % calculate index of centre slice within vol
    z_centre = round(central_position/vox_dim(3)) + 1;

    % calculate index of start slice within vol
    z_start = z_centre - centre_index + 1;

    % calculate index of end slice within vol
    z_end = z_start + no_slices - 1;

    % insert slab into vol between these boundaries 
    V(:,:,z_start:z_end) = slab;

    % rotate V 180 degrees around Z-axis to match convention in question
    % (figure 1) and previous answers
    V(:,:,z_start:z_end) = rot90(V(:,:,z_start:z_end),2);

    % specify origin of rotation to be the centre of the projection plane of
    % the slab
    origin = [image_dim(2)/2,image_dim(1)/2,z_centre];
    % specify range of slices to show within V in all 3 dimensions
    xslice= [];
    yslice= [];
    zslice= z_start:z_end;

    % open a new maximised window if indicated
    if strcmp(image,'new') == 1
        figure('units','normalized','outerposition',[0 0 1 1])
    end
    
    % retrieve the slab as a surface plot
    % display image in GUI axes, if axes is specified
    if nargin == 12
        destinaton = imagesc(projection_image); 
    else
        %otherwise plot normally
        destinaton = slice(V, xslice, yslice, zslice,'nearest');
    end

    %rotate around x axis
    rotate(destinaton,[1,0,0],alpha,origin)
    %rotate around y axis
    rotate(destinaton,[0,1,0],beta,origin)
    %rotate around z axis
    rotate(destinaton,[0,0,1],gamma,origin)

% For YZ slice
elseif strcmp(orientation,'Y-Z') == 1

    % calculate index of centre slice within vol
    x_centre = round(central_position/vox_dim(2)) + 1;

    % calculate index of start slice within vol
    x_start = x_centre - centre_index + 1;

    % calculate index of end slice within vol
    x_end = x_start + no_slices - 1; 

    % insert slab into vol between these boundaries 
    V(:,x_start:x_end,:) = slab;

    % rotate V 180 degrees around Z-axis to match convention in question
    % (figure 1) and previous answers
    V(:,x_start:x_end,:) = flipud(V(:,x_start:x_end,:));

    % specify origin of rotation to be the centre of the projection plane of
    % the slab
    origin = [x_centre,image_dim(1)/2,image_dim(3)/2];
    % specify range of slices to show within V in all 3 dimensions
    xslice= x_start:x_end;
    yslice= [];
    zslice= [];

    % open a new maximised window if indicated
    if strcmp(image,'new') == 1
        figure('units','normalized','outerposition',[0 0 1 1])
    end
    
    % retrieve the slab as a surface plot
    % display image in GUI axes, if axes is specified
    if nargin == 12
        destinaton = imagesc(projection_image);  
    else
        %otherwise plot normally
        destinaton = slice(V, xslice, yslice, zslice,'nearest');
    end

    %rotate around x axis
    rotate(destinaton,[1,0,0],alpha,origin)
    %rotate around y axis
    rotate(destinaton,[0,1,0],beta,origin)
    %rotate around z axis (NB negative here!)
    rotate(destinaton,[0,0,-1],gamma,origin)

% For XZ slice
elseif strcmp(orientation,'X-Z') == 1

    % calculate index of centre slice within vol
    y_centre = round(central_position/vox_dim(1)) + 1;

    % calculate index of start slice within vol
    y_start = y_centre - centre_index + 1;

    % calculate index of end slice within vol
    y_end = y_start + no_slices - 1; 

    % insert slab into vol between these boundaries 
    V(y_start:y_end,:,:) = slab;

    % rotate V 180 degrees around Z-axis to match convention in question
    % (figure 1) and previous answers
    V(y_start:y_end,:,:) = fliplr(V(y_start:y_end,:,:));

    % specify origin of rotation to be the centre of the projection plane of
    % the slab
    origin = [image_dim(2)/2,y_centre,image_dim(3)/2];
    % specify range of slices to show within V in all 3 dimensions
    xslice= [];
    yslice= y_start:y_end;
    zslice= [];

    % open a new maximised window if indicated
    if strcmp(image,'new') == 1
        figure('units','normalized','outerposition',[0 0 1 1])
    end
    
    % retrieve the slab as a surface plot
    % display image in GUI axes, if axes is specified
    if nargin == 12
        destinaton = imagesc(projection_image);  
    else
        %otherwise plot normally
        destinaton = slice(V, xslice, yslice, zslice,'nearest');
    end

    %rotate around x axis
    rotate(destinaton,[1,0,0],alpha,origin)
    %rotate around y axis
    rotate(destinaton,[0,1,0],beta,origin)
    %rotate around z axis (NB negative here!)
    rotate(destinaton,[0,0,-1],gamma,origin)
    
end

% General non-orientation specific formatting of the plot
colormap(gray)
% set edgecolour to none so detail of slice is visible
set(destinaton,'edgeColor','none')
% reverse the orientation of the Z and Y axes to match convention in
% question
set(gca,'Zdir','reverse','YDir','reverse')
% format plot, view and colourmap to match that of question
set(gca,'Box','on','BoxStyle','full')
set(gca,'View',[-125 30])
% label axes and set appropriate limits
xlabel('X / voxels')
ylabel('Y / voxels')
zlabel('Z / voxels')
xlim([1,image_dim(2)])
ylim([1,image_dim(1)])
zlim([1,image_dim(3)])
% adjust scaling to allow for voxel dimensions
daspect(1./[vox_dim(2) vox_dim(1) vox_dim(3)])

end

