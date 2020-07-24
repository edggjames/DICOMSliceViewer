function DisplayImage2D(projection_image, rotated_vox_dim, orientation, destinaton) %#ok<INUSD>
% DisplayImage2D displays a 2D projection image with the correct image 
% scaling according to orientation and rotated voxel dimensions, to be used
% in left hand axes of SliceViewer GUI
%
% DESCRIPTION:
% DisplayImage2D(projection_image, rotated_vox_dim, orientation). Plots
% image using 'imagesc' function, scales the axes accordingly, turns axes
% 'off', and uses a grey colourmap
%
% INPUTS:       
%       projection_image (double matrix) - a projection image that is
%           either central, maximum, minimum or median intensity projection
%           of a slab. This is returned by the function
%           ComputeProjectionImage
%
%       rotated_vox_dim (double vector) - a 3 x 1 vector containing the 
%           rotated voxel dimensions in the form [dx dy dz]
%
%       orientation (character string) - determines reference slice plane orientation
%           'X-Y' - XY plane, orthogonal to Z axis
%           'Y-Z' - YZ plane, orthogonal to X axis
%           'X-Z' - XZ plane, orthogonal to Y axis
%
%       destinaton (axes handle) - handle to the axes on which the projection 
%       image should be shown. This is an optional input argument. 
%
%OUTPUTS:
%       NONE.
%
% FUNCTION DEPENDENCIES:
%       NONE.
%
% AUTHOR:
%       Anonymised for MPHYGB24 MATLAB coursework assignment 2017/18
 
% To carry out orientation specific tasks:

% For XY slice
if strcmp(orientation,'X-Y') == 1 
    % display image in GUI axes, if axes is specified
    if nargin == 4
        destinaton = imagesc(projection_image); %#ok<NASGU>
    else
    %otherwise plot normally
    imagesc(projection_image)
    end
    % To carry out orientation specific scaling of image to allow for rotated 
    % voxel dimensions:
    daspect(1./[rotated_vox_dim(2) rotated_vox_dim(1) rotated_vox_dim(3)])
    % turn off axis and adjust colormap
    axis off
    colormap (gray)
    
% For YZ slice
elseif strcmp(orientation,'Y-Z') == 1
        % display image in GUI axes, if axes is specified
    if nargin == 4
        destinaton = imagesc(projection_image); %#ok<NASGU>
    else
    %otherwise plot normally
    imagesc(projection_image)
    end
    % To carry out orientation specific scaling of image to allow for rotated 
    % voxel dimensions:
    daspect(1./[rotated_vox_dim(1) rotated_vox_dim(3) rotated_vox_dim(2)])
    % turn off axis and adjust colormap
    axis off
    colormap (gray)
    
% For XZ slice
elseif strcmp(orientation,'X-Z') == 1
        % display image in GUI axes, if axes is specified
    if nargin == 4
        destinaton = imagesc(projection_image); %#ok<NASGU>
    else
    %otherwise plot normally
    imagesc(projection_image)
    end
    % To carry out orientation specific scaling of image to allow for rotated 
    % voxel dimensions:
    daspect(1./[rotated_vox_dim(2) rotated_vox_dim(3) rotated_vox_dim(1)])
    % turn off axis and adjust colormap
    axis off
    colormap (gray)
end

end

