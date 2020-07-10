function [slab_resize, centre_index, no_slices, rotated_vox_dim, stencil_3D] =  ...
    ComputeObliqueSlab(Image, orientation, central_position, res_1, res_2, ...
    res_3, method, alpha, beta, gamma, slab_thickness_mm)
% ComputeObliqueSlab computes an oblique slab from a 3D DICOM volume. Uses
% Gaussian blurring of all slices in slab when resampling at a lower
% in-plane resolution than rotated original image.
%
% DESCRIPTION: [slab_resize, centre_index, no_slices] = ComputeObliqueSlab ... 
%   (Image, orientation, central_position, res_1, res_2, res_3, method, ...
%   alpha, beta, gamma, slab_thickness_mm)
%       First calculates central slice of slab, then number of slices in
%       rotated slab.
%       Slab is composed of required number of slices and then undergoes 3D
%       interpolation (after Gaussian blurring if required) before being
%       returned. 
%       3D oblique slab is calculated on specified plane orientation, at
%       specified position, according to input 3D rotation angles, central
%       slice position and slab thickness
%       Gaussian blurring of slices is performed to avoid aliasing when sampling
%       at lower in-plane resolutions than original image in any direction.
%       This removes higher spatial frequencies and avoids pixellation and
%       other image artefacts. 
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
%       res1 (double scalar) - defines the resolution in mm in the output
%       slab of the first in-plane dimension, i.e x if XY plane is
%       selected, Y if YZ plane is selected, or x if XZ plane selected. NB
%       This is with reference to the rotated voxel dimension in this
%       function.
%
%       res2 (double scalar) - defines the resolution in mm in the output
%       slab of the second in-plane dimension, i.e y if XY plane is
%       selected, z if YZ plane is selected, or z if XZ plane selected.  NB
%       This is with reference to the rotated voxel dimension in this
%       function.
%       
%       res3 (double scalar) - defines the resolution in mm in the output
%       slab orthogonal to the slab, i.e z if XY plane is selected, x if
%       YZ plane is selected, or y if XZ plane selected.  NB This is with 
%       reference to the rotated voxel dimension in this function. This
%       must be less than or equal to half of the slice thickness. 
%
%       method (character string) - user defined interpolation method. One 
%       of:
%           nearest-neighbour ('nearest'),
%           linear ('linear'), or
%           cubic spline ('spline').
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
%       slab_thickness_mm (double scalar) - the thickness of the output
%       slab as specified by the user in mm. 
%
%       NB all three rotation angles are counterclockwise in a RHS
%       cartesian coordinate system, or clockwise in a LHS cartesian
%       coordinate system.
%
%OUTPUTS:
%       slab_resize (double matrix) - oblique slab from volume, resized to
%           be the same size as the uninterpolated central slice (which is
%           the projection plane). This is a 3D image volume. 
%
%       centre_index (integer scalar) - the central index of the number of
%           slices within the returned slab. This is 1 if the slab 
%
%       no_slices (integer scalar) - the number of slices in the slab
%
%       rotated_vox_dim (double vector) - a 3 x 1 vector of containing the
%           rotated voxel dimensions in the form [dy dx dz].
%
%       stencil_3D (double matrix) - a 3D stencil keeping a record of if
%       each individual voxel of the slab has been rotated out of the
%       volume or not
%
% FUNCTION DEPENDENCIES:
%       'ComputeObliqueSlice' - computes an oblique slice from 3D dicom
%       volume
%       'BlurSlices' - performs 1D or 2D Gaussian blurring on one or more
%       slices of a 3D volume, blurring each slice separately 
%
% AUTHOR:
%       Anonymised for MPHYGB24 MATLAB coursework assignment 2017/18

% To extract voxel dimensions in mm (used to scale axes in plots)
vox_dim = Image.VoxelDimensions; % [dy dx dz]

% compute central oblique slice, rotated voxel dimensions, and rotated
% voxel indices, with no interpolation within invoked function
[central_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, orientation,...
    central_position, res_1, res_2, method, alpha, beta, gamma, 'slab', central_position);

% calculate 2D size of central oblique slice in terms of voxels
size_central_slice = size(central_slice);

% To carry out orientation specific tasks:

% For XY slice
if strcmp(orientation,'X-Y') == 1
    % Check res_3 is less than or equal to slab thickness
    if res_3 > slab_thickness_mm
        error = errordlg({'Z resolution is greater than slab thickness.'; ...
            'Adjust resolution or slice thickness.'},...
            'Resolution Error','modal');
        disp('Oblique slab NOT computed.')
        % block program execution until user has clicked on modal error box
        uiwait(error)
        return
    end

    % thickness of slice is third rotated voxel_dimension
    slice_thickness_mm = rotated_vox_dim(3);

    % no of slices in slab is therefore
    no_slices = slab_thickness_mm / slice_thickness_mm;

    % convert this to an appropriate integer value
    if no_slices < 1
        no_slices = 1;
    else 
        no_slices = round(no_slices);
    end 

    % form a slab to hold all voxel itensities in 
    slab = zeros(size_central_slice(1), size_central_slice(2), no_slices);
    stencil_3D = ones(size_central_slice(1), size_central_slice(2), no_slices)*NaN;

    % assign central slice to centre of this slab
    centre_index = round(no_slices/2);
    slab (:,:,centre_index) = central_slice;
    stencil_3D ( :,:,centre_index) = stencil_2D;

    if no_slices ~= 1
        % assign remaining slices to slab either side of central slice
        temp_position = central_position - (centre_index*vox_dim(3));
        % Open a waitbar dialog box
        h = waitbar(0,'Assigning Slices to Slab...');
        for i = 1:no_slices
            temp_position = temp_position + vox_dim(3);
            %update waitbar dialog box
            waitbar(i / no_slices)
            if i ~=centre_index
                [temp_slice,~,stencil_2D] = ComputeObliqueSlice(Image, orientation, ...
                    temp_position, res_1, res_2, method, alpha, beta, gamma, ...
                    'slab', central_position);
                slab (:,:,i) = temp_slice;
                stencil_3D (:,:,i) = stencil_2D;
            end
        end
        %close waitbar dialog box
        close(h)
    end
    
    % Interpolate this slab in 3D according to res_1, res_2 and res_3
    % To convert resolutions from mm to number of rotated voxel dimensions
    x_res_vox = res_1/rotated_vox_dim(2);
    y_res_vox = res_2/rotated_vox_dim(1);
    z_res_vox = res_3/rotated_vox_dim(3);
    % Convert these into vectors of query points for interpolation
    x_query = 1:x_res_vox:size_central_slice(2);
    y_query = 1:y_res_vox:size_central_slice(1);
    z_query = 1:z_res_vox:no_slices;
    % set up a 3D grid of query points for interpolation
    [x_query,y_query,z_query] = meshgrid(x_query,y_query,z_query);

    % blur the slices of the slab if need be with the following sigma
    % values
    sigmax = ((sqrt(2*log(2)))/pi)*x_res_vox;
    sigmay = ((sqrt(2*log(2)))/pi)*y_res_vox;

    % Threshold for blurring here just over 1, as in-plane resolutions
    % in GUI are given to 3DP. Therefore this adjustment allows for
    % rounding errors. 

    if x_res_vox >= 1.001 && y_res_vox >= 1.001
        dimension = 'both';
    elseif x_res_vox >= 1.001 
        dimension = 'X';
    elseif y_res_vox >= 1.001 
        dimension = 'Y';
    end

    % blur these slice/s if necessary using function from Task 4
    if x_res_vox >= 1.001 || y_res_vox >= 1.001
        slab = BlurSlices(slab,'X-Y',sigmax,sigmay,'',dimension);
    end

    if no_slices == 1
        % replicate slab into 2 layers to allow for 3D interpolation
        slab(:,:,2) = squeeze(slab);
    end

    % perform 3D interpolation
    slab_interp = interp3(slab,x_query,y_query,z_query,method);

% For YZ slice
elseif strcmp(orientation,'Y-Z') == 1
    % Check res_3 is less than slab thickness
    if res_3 > slab_thickness_mm
        error = errordlg({'X resolution is greater than slab thickness.'; ...
            'Adjust resolution or slice thickness.'},...
            'Resolution Error','modal');
        disp('Oblique slab NOT computed.')
        % block program execution until user has clicked on modal error box
        uiwait(error)
        return
    end

    % thickness of slice is third rotated voxel_dimension
    slice_thickness_mm = rotated_vox_dim(2);

    % no of slices in slab is therefore
    no_slices = slab_thickness_mm / slice_thickness_mm;

    % convert this to an appropriate integer value
    if no_slices < 1
        no_slices = 1;
    else 
        no_slices = round(no_slices);
    end 
    
    % form a slab to hold all voxel itensities in 
    slab = zeros(size_central_slice(1), no_slices, size_central_slice(2));
    stencil_3D = ones(size_central_slice(1), no_slices, size_central_slice(2))*NaN;

    % assign central slice to centre of this slab
    centre_index = round(no_slices/2);
    slab (:,centre_index,:) = central_slice;
    stencil_3D (:,centre_index,:) = stencil_2D;

    if no_slices ~= 1
        % assign remaining slices to slab either side of central slice
        temp_position = central_position - (centre_index*vox_dim(2));
        % Open a waitbar dialog box
        h = waitbar(0,'Assigning Slices to Slab...');
        for i = 1:no_slices
            temp_position = temp_position + vox_dim(2);
            %update waitbar dialog box
            waitbar(i / no_slices)
            if i ~=centre_index
                [temp_slice,~,stencil_2D] = ComputeObliqueSlice(Image, orientation, ...
                    temp_position, res_1, res_2, method, alpha, beta, gamma, ...
                    'slab', central_position);
                slab (:,i,:) = temp_slice;
                stencil_3D (:,i,:) = stencil_2D;
            end
        end
        %close waitbar dialog box
        close(h)
    end
    
    % Interpolate this slab in 3D according to res_1, res_2 and res_3
    % To convert resolutions from mm to number of rotated voxel dimensions
    x_res_vox = res_3/rotated_vox_dim(2);
    y_res_vox = res_1/rotated_vox_dim(1);
    z_res_vox = res_2/rotated_vox_dim(3);
    % Convert these into vectors of query points for interpolation
    x_query = 1:x_res_vox:no_slices;
    y_query = 1:y_res_vox:size_central_slice(1);
    z_query = 1:z_res_vox:size_central_slice(2);
    % set up a 3D grid of query points for interpolation
    [x_query,y_query,z_query] = meshgrid(x_query,y_query,z_query);
    
    % blur the slices of the slab if need be with the following sigma
    % values
    sigmay = ((sqrt(2*log(2)))/pi)*y_res_vox;
    sigmaz = ((sqrt(2*log(2)))/pi)*z_res_vox;
    
    % Threshold for blurring here just over 1, as in-plane resolutions
    % in GUI are given to 3DP. Therefore this adjustment allows for
    % rounding errors. 

    if y_res_vox >= 1.001 && z_res_vox >= 1.001
        dimension = 'both';
    elseif y_res_vox >= 1.001 
        dimension = 'Y';
    elseif z_res_vox >= 1.001 
        dimension = 'Z';
    end

    % blur these slices if necessary using function from Task 4
    if y_res_vox >= 1.001 || z_res_vox >= 1.001
        slab = BlurSlices(slab,'Y-Z','',sigmay,sigmaz,dimension);
    end  
    
    if no_slices == 1
        % replicate slab into 2 layers to allow for 3D interpolation
        slab(:,2,:) = squeeze(slab);
    end
    
    % perform 3D interpolation
    slab_interp = interp3(slab,x_query,y_query,z_query,method);

% For XZ slice
elseif strcmp(orientation,'X-Z') == 1 
    
    % Check res_3 is less than or equal to slab thickness
    if res_3 > slab_thickness_mm
        error = errordlg({'Y resolution is greater than slab thickness.'; ...
            'Adjust resolution or slice thickness.'},...
            'Resolution Error','modal');
        disp('Oblique slab NOT computed.')
        % block program execution until user has clicked on modal error box
        uiwait(error)
        return
    end

    % thickness of slice is third rotated voxel_dimension
    slice_thickness_mm = rotated_vox_dim(1);

    % no of slices in slab is therefore
    no_slices = slab_thickness_mm / slice_thickness_mm;

    % convert this to an appropriate integer value
    if no_slices < 1
        no_slices = 1;
    else 
        no_slices = round(no_slices);
    end 
   
    % form a slab to hold all voxel itensities in 
    slab = zeros(no_slices,size_central_slice(1), size_central_slice(2));
    stencil_3D = ones(no_slices,size_central_slice(1), size_central_slice(2))*NaN;

    % assign central slice to centre of this slab
    centre_index = round(no_slices/2);
    slab (centre_index,:,:) = central_slice;
    stencil_3D (centre_index,:,:) = stencil_2D;
    
    if no_slices ~= 1
        % assign remaining slices to slab either side of central slice
        temp_position = central_position - (centre_index*vox_dim(1));
        % Open a waitbar dialog box
        h = waitbar(0,'Assigning Slices to Slab...');
        for i = 1:no_slices
            temp_position = temp_position + vox_dim(1);
            %update waitbar dialog box
            waitbar(i / no_slices)
            if i ~=centre_index
                [temp_slice,~,stencil_2D] = ComputeObliqueSlice(Image, orientation, ...
                    temp_position, res_1, res_2, method, alpha, beta, gamma, ...
                    'slab', central_position);
                slab (i,:,:) = temp_slice;
                stencil_3D (i,:,:) = stencil_2D;
            end
        end
        %close waitbar dialog box
        close(h)
    end
    
    % Interpolate the slab

    % interpolate this slab in 3D according to res_1, res_2 and res_3
    % To convert resolutions from mm to number of rotated voxel dimensions
    x_res_vox = res_1/rotated_vox_dim(2);
    y_res_vox = res_3/rotated_vox_dim(1);
    z_res_vox = res_2/rotated_vox_dim(3);
    % Convert these into vectors of query points for interpolation
    x_query = 1:x_res_vox:size_central_slice(1);
    y_query = 1:y_res_vox:no_slices;
    z_query = 1:z_res_vox:size_central_slice(2);
    % set up a 3D grid of query points for interpolation
    [x_query,y_query,z_query] = meshgrid(x_query,y_query,z_query);
    
    if no_slices == 1
        % replicate slab into 2 layers to allow for 3D interpolation
        slab(2,:,:) = squeeze(slab);
    end
      
    % blur the slices of the slab if need be with the following sigma
    % values
    sigmax = ((sqrt(2*log(2)))/pi)*x_res_vox;
    sigmaz = ((sqrt(2*log(2)))/pi)*z_res_vox;
    
    % Threshold for blurring here just over 1, as in-plane resolutions
    % in GUI are given to 3DP. Therefore this adjustment allows for
    % rounding errors. 

    if x_res_vox >= 1.001 && z_res_vox >= 1.001
        dimension = 'both';
    elseif x_res_vox >= 1.001 
        dimension = 'X';
    elseif z_res_vox >= 1.001 
        dimension = 'Z';
    end

    % blur these slices if necessary using function from Task 4
    if y_res_vox >= 1.001 || z_res_vox >= 1.001
        slab = BlurSlices(slab,'X-Z',sigmax,'',sigmaz,dimension);
    end  
    
    % perform 3D interpolation
    slab_interp = interp3(slab,x_query,y_query,z_query,method);
 
end
    
% resize slab_interp to original slab dimensions
% if slab_interp is 3D
if no_slices > 1
    slab_resize = imresize3(slab_interp,size(slab));
% if slab_interp is 2D
elseif no_slices == 1
    slab_interp = squeeze(slab_interp);
    slab_resize = imresize(slab_interp,size(central_slice));
    % also ensure stencil_3D is only 2D if there is only 1 slice
    stencil_3D = squeeze(stencil_3D);
end

end

