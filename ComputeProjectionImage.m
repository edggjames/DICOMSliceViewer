function projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices)
% Performs projection operations on a slab extracted from an image volume.
% The projection plane is the central slice plane, and the projection angle
% is perpendicular to the projection plane. 

% DESCRIPTION:
% projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
%    centre_index, view, no_slices). The function supports one of four
%    projection operations (detailed below in 'inputs').
%
% INPUTS:       
%       slab (double matrix) - oblique slab from volume, resized to
%           be the same size as the uninterpolated central slice (which is
%           the projection plane). This is a 3D image volume. 
%
%       stencil_3D (double matrix) - a 3D stencil which is a record of if
%       each individual voxel of the slab has been rotated out of the
%       volume or not.
%
%       orientation (character string) - determines reference slice plane orientation
%           'X-Y' - XY plane, orthogonal to Z axis
%           'Y-Z' - YZ plane, orthogonal to X axis
%           'X-Z' - XZ plane, orthogonal to Y axis
%
%       centre_index (integer scalar) - the central index of the number of
%           slices within the returned slab. This is 1 if the slab
%
%       view (character string) - this determined the type of projection
%           image returned. Four types are supported by this function:
%               1) 'central' - central slice plane
%               2) 'max' - maximum intensity projection
%               3) 'min' - minimum intensity projection
%               4) 'median' - median intensity projection
%
%       no_slices (integer scalar) - the number of slices in the slab
%
%OUTPUTS:
%       projection_image (double matrix) - a projection image that is
%           either central, maximum, minimum or median intensity projection
%           of a slab, using the central slice plane as the projection
%           plane.
%
% FUNCTION DEPENDENCIES:
%       'StripBorderObliqueSlice' - removes the border of a slice which has
%           been rotated out of the image volume, by reference to
%           stencil_2D
%
% AUTHOR:
%       Anonymised for MPHYGB24 MATLAB coursework assignment 2017/18

% ensure that slab is the same size as stencil_3D
if no_slices == 1
    slab = imresize(slab, size(stencil_3D));
elseif no_slices > 1
    slab = imresize3(slab, size(stencil_3D));
end
% convert all values in slab that have been rotated out of volume to NaNs
slab = slab.*stencil_3D;

% Perform orientation specific tasks

if strcmp(orientation,'X-Y') == 1  

    % Perform projection view specific tasks
    
    if strcmp(view,'central') == 1
        % return central slice of slab (i.e. projection plane)
        projection_image = squeeze(slab(:,:,centre_index));
        % also form a 2D stencil for central slice view
        stencil_centre = squeeze(stencil_3D(:,:,centre_index));
        % strip low value border from image
        projection_image = StripBorderObliqueSlice(projection_image,stencil_centre);
    else
        % form a 2D stencil which projects all voxels rotated within the volume
        % onto one plane, omitting NaNs that have been rotated out of volume
        stencil_2D = squeeze(max(stencil_3D,[],3,'omitnan'));
        % This is then the stencil for max, min and median projection
        % calculations
    end
    
    if strcmp(view,'max') == 1
        % compute maximum intensity projection of slab along third dimension
        projection_image = squeeze(max(slab,[],3,'omitnan'));
        % strip low value border from image
        projection_image = StripBorderObliqueSlice(projection_image,stencil_2D);
    elseif strcmp(view,'min') == 1
        % compute minimum intensity projection of slab along third dimension
        projection_image = squeeze(min(slab,[],3,'omitnan'));
        % strip low value border from image
        projection_image = StripBorderObliqueSlice(projection_image,stencil_2D);
    elseif strcmp(view,'median') == 1
        % compute median intensity projection of slab along third dimension
        projection_image = squeeze(median(slab,3,'omitnan'));
        % strip low value border from image
        projection_image = StripBorderObliqueSlice(projection_image,stencil_2D);
    end

    
elseif strcmp(orientation,'Y-Z') == 1

    % Perform projection view specific tasks
    
    if strcmp(view,'central') == 1
        if no_slices > 1 
            % return central slice of slab (i.e. projection plane)
            projection_image = squeeze(slab(:,centre_index,:));
            % form a 2D stencil for central slice
            stencil_centre = squeeze(stencil_3D(:,centre_index,:));
            % strip low value border from image
            projection_image = StripBorderObliqueSlice(projection_image,stencil_centre);
        elseif no_slices == 1
            projection_image = StripBorderObliqueSlice(slab,stencil_3D);            
        end
    else
        % form a 2D stencil which projects all voxels rotated within the volume
        % onto one plane, omitting NaNs that have been rotated out of volume
        if no_slices > 1
           stencil_2D = squeeze(max(stencil_3D,[],2,'omitnan'));
        elseif no_slices == 1
           stencil_2D = stencil_3D;
        end
    end
    
    if strcmp(view,'max') == 1
        if no_slices > 1
            % compute maximum intensity projection of slab along second dimension
            projection_image = squeeze(max(slab,[],2,'omitnan'));
        elseif no_slices == 1
            % compute maximum intensity projection of slab along third dimension
            projection_image = squeeze(max(slab,[],3,'omitnan'));
        end
        % strip low value border from image
        projection_image = StripBorderObliqueSlice(projection_image,stencil_2D);
    elseif strcmp(view,'min') == 1
        if no_slices > 1
            % compute minimum intensity projection of slab along second dimension
            projection_image = squeeze(min(slab,[],2,'omitnan'));
        elseif no_slices == 1
            % compute minimum intensity projection of slab along third dimension
            projection_image = squeeze(min(slab,[],3,'omitnan'));
        end
        % strip low value border from image
        projection_image = StripBorderObliqueSlice(projection_image,stencil_2D);
    elseif strcmp(view,'median') == 1
        if no_slices > 1
            % compute median intensity projection of slab along second dimension
            projection_image = squeeze(median(slab,2,'omitnan'));
        elseif no_slices == 1
            % compute median intensity projection of slab along third dimension
            projection_image = squeeze(median(slab,3,'omitnan'));
        end
        % strip low value border from image
        projection_image = StripBorderObliqueSlice(projection_image,stencil_2D);
    end

    
elseif strcmp(orientation,'X-Z') == 1
    
    % Perform projection view specific tasks
    
    if strcmp(view,'central') == 1
        if no_slices > 1
            % return central slice of slab (i.e. projection plane)
            projection_image = squeeze(slab(centre_index,:,:));
            % form a 2D stencil for central slice
            stencil_centre = squeeze(stencil_3D(centre_index,:,:));
            % strip low value border from image
            projection_image = StripBorderObliqueSlice(projection_image,stencil_centre);
        elseif no_slices == 1
            projection_image = StripBorderObliqueSlice(slab,stencil_3D);            
        end
    else
        % form a 2D stencil which projects all voxels rotated within the volume
        % onto one plane, omitting NaNs that have been rotated out of volume
        if no_slices > 1
           stencil_2D = squeeze(max(stencil_3D,[],1,'omitnan'));
        elseif no_slices == 1
           stencil_2D = stencil_3D;
        end
    end
    
    if strcmp(view,'max') == 1
        if no_slices > 1
            % compute maximum intensity projection of slab along first dimension
            projection_image = squeeze(max(slab,[],1,'omitnan'));
        elseif no_slices == 1
            % compute maximum intensity projection of slab along third dimension
            projection_image = squeeze(max(slab,[],3,'omitnan'));
        end
        % strip low value border from image
        projection_image = StripBorderObliqueSlice(projection_image,stencil_2D);
    elseif strcmp(view,'min') == 1
        if no_slices > 1
            % compute minimum intensity projection of slab along first dimension
            projection_image = squeeze(min(slab,[],1,'omitnan'));
        elseif no_slices == 1
            % compute minimum intensity projection of slab along third dimension
            projection_image = squeeze(min(slab,[],3,'omitnan'));
        end
        % strip low value border from image
        projection_image = StripBorderObliqueSlice(projection_image,stencil_2D);
    elseif strcmp(view,'median') == 1
        if no_slices > 1
            % compute median intensity projection of slab along first dimension
            projection_image = squeeze(median(slab,1,'omitnan'));
        elseif no_slices == 1
            % compute median intensity projection of slab along third dimension
            projection_image = squeeze(median(slab,3,'omitnan'));
        end
        % strip low value border from image
        projection_image = StripBorderObliqueSlice(projection_image,stencil_2D);
    end
    
end

% ensure projection_image is 2D by removing singleton dimension if it
% exists
projection_image = squeeze(projection_image);

end

