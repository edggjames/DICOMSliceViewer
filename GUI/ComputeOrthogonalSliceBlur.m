function slice = ComputeOrthogonalSliceBlur(Image, orientation, position, res_1, res_2, method)
% ComputeOrthogonalSlice computes an orthogonal slice from a 3D DICOM
% volume. Uses Gaussian blurring of slices to when resampling at a lower
% in-plane resolution. 
%
% DESCRIPTION: slice = ComputeOrthogonalSlice(Image, orientation, ...
%       position, res_1, res_2, method)
%       Calculates volume dimensions in terms of both voxels and mm
%       Extracts voxel dimensions
%       2D slice is calculated on specified plane orientation and at
%       specified position using 3D interpolation (performed using user
%       specified interpolation method)
%       Gaussian blurring of slices is performed to avoid aliasing when sampling
%       at lower in-plane resolutions than original image in any direction.
%       This removes higher spatial frequencies and avoid pixellation and
%       other image artefacts. 
%       2D slice has in-plane pixel resolutions as specified in function
%       inputs
%
% INPUTS:
%       Image (1 x 1 structure with two fields) - 
%           .ImageData (double matrix) - of dimensions (number of rows, 
%           number of columns,number of slices) containing the voxel grey 
%           level values
%           .VoxelDimensions (double vector) - a 1 by 3 vector containing
%           the (y,x,z) voxel dimensions in mm, respectively
%       
%       orientation (character string) - determines slice plane orientation
%           'X-Y' - XY plane, orthogonal to Z axis
%           'Y-Z' - YZ plane, orthogonal to X axis
%           'X-Z' - XZ plane, orthogonal to Y axis
%
%       position (double scalar) - determines position of slice in mm along
%       axis orthogonal to the slice plane:
%           For the Z slice position, zero is toward the head of the 
%               patient.
%           For the Y slice position, zero is toward the front of the 
%               patient. 
%           For the X slice position, zero is toward the left of the
%               patient (for this volume this end of the voxel range is
%               where the arm is shown).
%
%       res1 (double scalar) - defines the resolution in mm in the output
%       slice of the first in-plane dimension, i.e x if XY plane is
%       selected, Y if YZ plane is selected, or x if XZ plane selected. NB
%       decreasing resolution below original respective voxel dimension
%       will not improve results.
%
%       res2 (double scalar) - defines the resolution in mm in the output
%       slice of the second in-plane dimension, i.e y if XY plane is
%       selected, z if YZ plane is selected, or z if XZ plane selected.  NB
%       decreasing resolution below original respective voxel dimension
%       will not improve results.
%       
%       method (character string) - user defined interpolation method. One 
%       of:
%           nearest-neighbour ('nearest'),
%           linear ('linear'), or
%           cubic spline ('spline').
%
%OUTPUTS:
%       Slice (double matrix) - orthogonal slice from volume. Slice 
%       dimensions match original volumne in-plane dimensions due to use of 
%       'imresize' at end of function (even if sampled at lower resolution)
%
% FUNCTION DEPENDENCIES:
%       'BlurSlices' - performs 1D or 2D Gaussian blurring on one or more
%       slices of a 3D volume, blurring each slice separately 
%
% AUTHOR:
%       Anonymised for MPHYGB24 MATLAB coursework assignment 2017/18

% To assign 3D image intensities to a matrix variable
vol = Image.ImageData;

% To get size of image as 3D vector in form of (y,x,z) i.e.
% (number rows, number of columns, number of slices)
image_dim = size(Image.ImageData);
% Then assign as follows
no_rows = image_dim(1);
no_cols = image_dim(2);
no_slices = image_dim(3);

% To get size of voxel as 3D vector in mm in form of (y,x,z) i.e.
% (voxel height, voxel width, voxel length)
vox_dim = Image.VoxelDimensions;
% Then assign as follows
vox_height = vox_dim(1);
vox_width = vox_dim(2);
vox_length = vox_dim(3);

% To get size of image as 3D vector in mm in form of (y,x,z) i.e.
% (image height, image width, image length)
image_size = (image_dim - [1 1 1]).*vox_dim;
% Assume first slice is centred around 0 mm in respective direction
% Then assign as follows
image_height = image_size(1);
image_width = image_size(2);
image_length = image_size(3);

% To carry out orientation specific tasks:

% For XY slice
if strcmp(orientation,'X-Y') == 1
    
    % 1) To convert resolutions from mm to number of voxel dimensions for
    % 3D interpolation for in-plane dimensions
    x_res_vox = res_1/vox_width;
    y_res_vox = res_2/vox_height;
    % Convert these into vectors of query points for interpolation
    x_query = 1:x_res_vox:no_cols;
    y_query = 1:y_res_vox:no_rows;
     
    % 2) Form expressions for sigmax, sigmay and the dimensions to perform
    % Gaussian blurring in
    sigmax = ((sqrt(2*log(2)))/pi)*x_res_vox;
    sigmay = ((sqrt(2*log(2)))/pi)*y_res_vox;

    % Threshold for blurring here is just over 1, as in-plane resolutions
    % in GUI are given to 3DP. Therefore this adjustment allows for
    % rounding errors. 

    if x_res_vox >= 1.001 && y_res_vox >= 1.001
        dimension = 'both';
    elseif x_res_vox >= 1.001 
        dimension = 'X';
    elseif y_res_vox >= 1.001 
        dimension = 'Y';
    end
    
    % 3) Check position is within bounds
    if position < 0 || position > image_length 
        % out of bound measurement has been requested
        error = errordlg({'Cannot display X-Y slice outside volume dimensions'; ...
            ['Z Slice Position cannot be less than 0.000 or greater than ' ...
            , sprintf('%.3f',image_length)]},'Interpolation Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
        return
    end
    
    % 4) convert position to a slice number to index original volume with
    % Assuming slice 1 is centred around 0 mm
    z_query = (position/vox_length) + 1;
     
    % 5) Interpolate according to string contained in method
    if strcmp(method,'nearest') == 1
        % To make this as fast as possible - just use the two slices that
        % neighbour the z_query point, and then perform 3D nearest
        % neighbour interpolation between these two slices (this requires
        % at least 2 sample points in each dimension)
        if round(z_query) == image_dim(3)
            % sample the last 2 slices and z_query point is second slice
            slice_vector = [round(z_query)-1, round(z_query)];
            z_query_nearest = 2;
        else
            % sample slice either side and z_query point is first slice
            slice_vector = [round(z_query), round(z_query)+1];
            z_query_nearest = 1;
        end
        % sample volume by this slice vector
        reduced_vol = vol(:,:,slice_vector);
        % blur these slices if necessary using function from Task 4
        if x_res_vox >= 1.001 || y_res_vox >= 1.001
            reduced_vol = BlurSlices(reduced_vol,'X-Y',sigmax,sigmay,'',dimension);
        end    
        % set up a 3D grid of query points for interpolation
        [x_query,y_query,z_query_nearest] = meshgrid(x_query,y_query,z_query_nearest);
        % perform 3D nearest neighbour interpolation
        slice = interp3(reduced_vol,x_query,y_query,z_query_nearest,method);
        % ensure that this is a 2D slice with no singleton dimensions
        slice = squeeze(slice);
        % resize the slice to the original dimensions (this compensates for
        % sampling at lower resolutions than original voxel dimensions)
        slice = imresize(slice,[no_rows, no_cols]);
        
    elseif strcmp(method,'linear') == 1
        % To make this as fast as possible - just use the two slices that
        % neighbour the z_query point, and then perform 3D linear
        % interpolation between these two slices
        if floor(z_query) == image_dim(3)
            % sample the last 2 slices and z_query point is second slice
            slice_vector = [floor(z_query)-1, floor(z_query)];
            z_query_linear = 2;
        else
            % sample slice either side and z_query point is appropriately
            % spaced between these two slices
            slice_vector = [floor(z_query), floor(z_query)+1];
            z_query_linear = z_query - floor(z_query) + 1;
        end
        % sample volume by this slice vector
        reduced_vol = vol(:,:,slice_vector);
        % blur these slices if necessary using function from Task 4
        if x_res_vox >= 1.001 || y_res_vox >= 1.001
            reduced_vol = BlurSlices(reduced_vol,'X-Y',sigmax,sigmay,'',dimension);
        end 
        % set up a 3D grid of query points for interpolation
        [x_query,y_query,z_query_linear] = meshgrid(x_query,y_query,z_query_linear);
        % perform 3D linear interpolation
        slice = interp3(reduced_vol,x_query,y_query,z_query_linear,method);
        % ensure that this is a 2D slice with no singleton dimensions
        slice = squeeze(slice);
        % resize the slice to the original dimensions (this compensates for
        % sampling at lower resolutions than original voxel dimensions)
    	slice = imresize(slice,[no_rows, no_cols]);
        
    elseif strcmp(method,'spline') == 1
        % The MATLAB cubic spline interpolation algorithm requires 4 grid
        % points in each dimension. Therefore, to make this as fast as 
        % possible just use the two slices either side of the z_query 
        % point, and then perform interpolation on these 4 slices
        if floor(z_query) == 1
            % sample the first 4 slices and z_query point is unchanged
            slice_vector = [1, 2, 3, 4];
            z_query_spline = z_query;
        elseif floor(z_query) == image_dim(3) - 1
            % sample the last 4 slices and z_query point is between last 
            % two slices
            slice_vector = image_dim(3) - 3 : image_dim(3);
            z_query_spline = z_query - floor(z_query) + 3;
        elseif floor(z_query) == image_dim(3)
            % sample the last 4 slices and z_query point is last slice
            slice_vector = image_dim(3) - 3 : image_dim(3);
            z_query_spline = 4;
        else
            % sample 2 slices either side and z_query point is 
            % appropriately spaced between these 4 slices
            slice_vector = floor(z_query)-1 : floor(z_query)+2;
            z_query_spline = z_query - floor(z_query) + 2;
        end
        % sample volume by this slice vector
        reduced_vol = vol(:,:,slice_vector);
        % blur these slices if necessary using function from Task 4
        if x_res_vox >= 1.001 || y_res_vox >= 1.001
            reduced_vol = BlurSlices(reduced_vol,'X-Y',sigmax,sigmay,'',dimension);
        end 
        % set up a 3D grid of query points for interpolation
        [x_query,y_query,z_query_spline] = meshgrid(x_query,y_query,z_query_spline);
        % perform 3D linear interpolation
        slice = interp3(reduced_vol,x_query,y_query,z_query_spline,method);
        % ensure that this is a 2D slice with no singleton dimensions
        slice = squeeze(slice);
        % resize the slice to the original dimensions (this compensates for
        % sampling at lower resolutions than original voxel dimensions)
    	slice = imresize(slice,[no_rows, no_cols]);
    end
    
% For YZ slice
elseif strcmp(orientation,'Y-Z') == 1
        
    % 1) To convert resolutions from mm to number of voxel dimensions for
    % 3D interpolation for in-plane dimensions
    y_res_vox = res_1/vox_height;
    z_res_vox = res_2/vox_length;
    % Convert these into vectors of query points for interpolation
    y_query = 1:y_res_vox:no_rows;
    z_query = 1:z_res_vox:no_slices;
    
    % 2) Form expressions for sigmay, sigmaz and the dimensions to perform
    % Gaussian blurring in
    % is delta here in terms of ratio of new to original voxel dimension??
    sigmay = ((sqrt(2*log(2)))/pi)*y_res_vox;
    sigmaz = ((sqrt(2*log(2)))/pi)*z_res_vox;

    if y_res_vox >= 1.001 && z_res_vox >= 1.001
        dimension = 'both';
    elseif y_res_vox >= 1.001 
        dimension = 'Y';
    elseif z_res_vox >= 1.001 
        dimension = 'Z';
    end

    % 3) Check position is within bounds
    if position < 0 || position > image_width 
        % out of bound measurement has been requested
        error = errordlg({'Cannot display Y-Z slice outside volume dimensions'; ...
            ['X Slice Position cannot be less than 0.000 or greater than ' ...
            , sprintf('%.3f',image_width)]},'Interpolation Error','modal');
        %block program execution until user has clicked on modal error box
        uiwait(error)
        return
    end
    
    % 4) convert position to a slice number to index original volume with
    x_query = (position/vox_width) + 1;
    % flip position around midpoint for correct orientation in X direction
    % (i.e. assumption 0 mm occurs at centre of last slice is made in this 
    % coursework)
    x_query = image_dim(2) - x_query + 1;
 
    % 5) Interpolate according to string contained in method
    if strcmp(method,'nearest') == 1
        % To make this as fast as possible - just use the two slices that
        % neighbour the x_query point, and then perform 3D nearest
        % neighbour interpolation between these two slices
        if round(x_query) == image_dim(2)
            % sample the last 2 slices and x_query point is second slice
            slice_vector = [round(x_query)-1, round(x_query)];
            x_query_nearest = 2;
        else
            % sample slice either side and x_query point is first slice
            slice_vector = [round(x_query), round(x_query)+1];
            x_query_nearest = 1;
        end
        % sample volume by this slice vector
        reduced_vol = vol(:,slice_vector,:);
        % blur these slices if necessary using function from Task 4
        if y_res_vox >= 1.001 || z_res_vox >= 1.001
            reduced_vol = BlurSlices(reduced_vol,'Y-Z','',sigmay,sigmaz,dimension);
        end 
        % set up a 3D grid of query points for interpolation
        [x_query_nearest,y_query,z_query] = meshgrid(x_query_nearest,y_query,z_query);
        % perform 3D nearest neighbour interpolation
        slice = interp3(reduced_vol,x_query_nearest,y_query,z_query,method);
        % ensure that this is a 2D slice with no singleton dimensions
        slice = squeeze(slice);
        % resize the slice to the original dimensions (this compensates for
        % sampling at lower resolutions than original voxel dimensions)
        slice = imresize(slice,[no_rows, no_slices]);
        
    elseif strcmp(method,'linear') == 1
        % To make this as fast as possible - just use the two slices that
        % neighbour the x_query point, and then perform 3D linear
        % interpolation between these two slices
        if floor(x_query) == image_dim(2)
            % sample the last 2 slices and x_query point is second slice
            slice_vector = [floor(x_query)-1, floor(x_query)];
            x_query_linear = 2;
        else
            % sample slice either side and x_query point is appropriately
            % spaced between these two slices
            slice_vector = [floor(x_query), floor(x_query)+1];
            x_query_linear = x_query - floor(x_query) + 1;
        end
        % sample volume by this slice vector
        reduced_vol = vol(:,slice_vector,:);
        % blur these slices if necessary using function from Task 4
        if y_res_vox >= 1.001 || z_res_vox >= 1.001
            reduced_vol = BlurSlices(reduced_vol,'Y-Z','',sigmay,sigmaz,dimension);
        end 
        % set up a 3D grid of query points for interpolation
        [x_query_linear,y_query,z_query] = meshgrid(x_query_linear,y_query,z_query);
        % perform 3D linear interpolation
        slice = interp3(reduced_vol,x_query_linear,y_query,z_query,method);
        % ensure that this is a 2D slice with no singleton dimensions
        slice = squeeze(slice);
        % resize the slice to the original dimensions (this compensates for
        % sampling at lower resolutions than original voxel dimensions)
    	slice = imresize(slice,[no_rows, no_slices]);
        
    elseif strcmp(method,'spline') == 1
        % The MATLAB cubic spline interpolation algorithm requires 4 grid
        % points in each dimension. Therefore, to make this as fast as 
        % possible just use the two slices either side of the x_query 
        % point, and then perform interpolation on these 4 slices
        if floor(x_query) == 1
            % sample the first 4 slices and x_query point is unchanged
            slice_vector = [1, 2, 3, 4];
            x_query_spline = x_query;
        elseif floor(x_query) == image_dim(2) - 1
            % sample the last 4 slices and x_query point is between last 
            % two slices
            slice_vector = image_dim(2) - 3 : image_dim(2);
            x_query_spline = x_query - floor(x_query) + 3;
        elseif floor(x_query) == image_dim(2)
            % sample the last 4 slices and x_query point is last slice
            slice_vector = image_dim(2) - 3 : image_dim(2);
            x_query_spline = 4;
        else
            % sample 2 slices either side and x_query point is 
            % appropriately spaced between these 4 slices
            slice_vector = floor(x_query)-1 : floor(x_query)+2;
            x_query_spline = x_query - floor(x_query) + 2;
        end
        % sample volume by this slice vector
        reduced_vol = vol(:,slice_vector,:);
        % blur these slices if necessary using function from Task 4
        if y_res_vox >= 1.001 || z_res_vox >= 1.001
            reduced_vol = BlurSlices(reduced_vol,'Y-Z','',sigmay,sigmaz,dimension);
        end 
        % set up a 3D grid of query points for interpolation
        [x_query_spline,y_query,z_query] = meshgrid(x_query_spline,y_query,z_query);
        % perform 3D linear interpolation
        slice = interp3(reduced_vol,x_query_spline,y_query,z_query,method);
        % ensure that this is a 2D slice with no singleton dimensions
        slice = squeeze(slice);
        % resize the slice to the original dimensions (this compensates for
        % sampling at lower resolutions than original voxel dimensions)
    	slice = imresize(slice,[no_rows, no_slices]);
    end
    
% For XZ slice
elseif strcmp(orientation,'X-Z') == 1
    
    % 1) To convert resolutions from mm to number of voxel dimensions for
    % 3D interpolation for in-plane dimensions
    x_res_vox = res_1/vox_width;
    z_res_vox = res_2/vox_length;
    % Convert these into vectors of query points for interpolation
    x_query = 1:x_res_vox:no_cols;
    z_query = 1:z_res_vox:no_slices;
    
    % 2) Form expressions for sigmax, sigmaz and the dimensions to perform
    % Gaussian blurring in
    % is delta here in terms of ratio of new to original voxel dimension??
    sigmax = ((sqrt(2*log(2)))/pi)*x_res_vox;
    sigmaz = ((sqrt(2*log(2)))/pi)*z_res_vox;

    if x_res_vox >= 1.001 && z_res_vox >= 1.001
        dimension = 'both';
    elseif x_res_vox >= 1.001 
        dimension = 'X';
    elseif z_res_vox >= 1.001 
        dimension = 'Z';
    end
    
    % 3) Check position is within bounds
    if position < 0 || position > image_height 
        % out of bound measurement has been requested
        error = errordlg({'Cannot display X-Z slice outside volume dimensions'; ...
            ['Y Slice Position cannot be less than 0.000 or greater than ' ...
            , sprintf('%.3f',image_height)]},'Interpolation Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
        return
    end
    
    % 4) convert position to a slice number to index original volume with
    y_query = (position/vox_height) + 1;
    % flip position around midpoint for correct orientation in Y direction
    % (i.e. assumption 0 mm occurs at centre of last slice is made in this 
    % coursework)
    y_query = image_dim(1) - y_query + 1;

    % 5) Interpolate according to string contained in method
    if strcmp(method,'nearest') == 1
        % To make this as fast as possible - just use the two slices that
        % neighbour the y_query point, and then perform 3D nearest
        % neighbour interpolation between these two slices
        if round(y_query) == image_dim(1)
            % sample the last 2 slices and y_query point is second slice
            slice_vector = [round(y_query)-1, round(y_query)];
            y_query_nearest = 2;
        else
            % sample slice either side and y_query point is first slice
            slice_vector = [round(y_query), round(y_query)+1];
            y_query_nearest = 1;
        end
        % sample volume by this slice vector
        reduced_vol = vol(slice_vector,:,:);
        % blur these slices if necessary using function from Task 4
        if x_res_vox >= 1.001 || z_res_vox >= 1.001
            reduced_vol = BlurSlices(reduced_vol,'X-Z',sigmax,'',sigmaz,dimension);
        end 
        % set up a 3D grid of query points for interpolation
        [x_query,y_query_nearest,z_query] = meshgrid(x_query,y_query_nearest,z_query);
        % perform 3D nearest neighbour interpolation
        slice = interp3(reduced_vol,x_query,y_query_nearest,z_query,method);
        % ensure that this is a 2D slice with no singleton dimensions
        slice = squeeze(slice);
        % resize the slice to the original dimensions (this compensates for
        % sampling at lower resolutions than original voxel dimensions)
        slice = imresize(slice,[no_cols, no_slices]);
        
    elseif strcmp(method,'linear') == 1
        % To make this as fast as possible - just use the two slices that
        % neighbour the y_query point, and then perform 3D linear
        % interpolation between these two slices
        if floor(y_query) == image_dim(1)
            % sample the last 2 slices and y_query point is second slice
            slice_vector = [floor(y_query)-1, floor(y_query)];
            y_query_linear = 2;
        else
            % sample slice either side and y_query point is appropriately
            % spaced between these two slices
            slice_vector = [floor(y_query), floor(y_query)+1];
            y_query_linear = y_query - floor(y_query) + 1;
        end
        % sample volume by this slice vector
        reduced_vol = vol(slice_vector,:,:);
        % blur these slices if necessary using function from Task 4
        if x_res_vox >= 1.001 || z_res_vox >= 1.001
            reduced_vol = BlurSlices(reduced_vol,'X-Z',sigmax,'',sigmaz,dimension);
        end 
        % set up a 3D grid of query points for interpolation
        [x_query,y_query_linear,z_query] = meshgrid(x_query,y_query_linear,z_query);
        % perform 3D linear interpolation
        slice = interp3(reduced_vol,x_query,y_query_linear,z_query,method);
        % ensure that this is a 2D slice with no singleton dimensions
        slice = squeeze(slice);
        % resize the slice to the original dimensions (this compensates for
        % sampling at lower resolutions than original voxel dimensions)
    	slice = imresize(slice,[no_cols, no_slices]);
        
    elseif strcmp(method,'spline') == 1
        % The MATLAB cubic spline interpolation algorithm requires 4 grid
        % points in each dimension. Therefore, to make this as fast as 
        % possible just use the two slices either side of the y_query 
        % point, and then perform interpolation on these 4 slices
        if floor(y_query) == 1
            % sample the first 4 slices and y_query point is unchanged
            slice_vector = [1, 2, 3, 4];
            y_query_spline = y_query;
        elseif floor(y_query) == image_dim(1) - 1
            % sample the last 4 slices and y_query point is between last 
            % two slices
            slice_vector = image_dim(1) - 3 : image_dim(1);
            y_query_spline = y_query - floor(y_query) + 3;
        elseif floor(y_query) == image_dim(1)
            % sample the last 4 slices and y_query point is last slice
            slice_vector = image_dim(1) - 3 : image_dim(1);
            y_query_spline = 4;
        else
            % sample 2 slices either side and y_query point is 
            % appropriately spaced between these 4 slices
            slice_vector = floor(y_query)-1 : floor(y_query)+2;
            y_query_spline = y_query - floor(y_query) + 2;
        end
        % sample volume by this slice vector
        reduced_vol = vol(slice_vector,:,:);
        % blur these slices if necessary using function from Task 4
        if x_res_vox >= 1.001 || z_res_vox >= 1.001
            reduced_vol = BlurSlices(reduced_vol,'X-Z',sigmax,'',sigmaz,dimension);
        end 
        % set up a 3D grid of query points for interpolation
        [x_query,y_query_spline,z_query] = meshgrid(x_query,y_query_spline,z_query);
        % perform 3D linear interpolation
        slice = interp3(reduced_vol,x_query,y_query_spline,z_query,method);
        % ensure that this is a 2D slice with no singleton dimensions
        slice = squeeze(slice);
        % resize the slice to the original dimensions (this compensates for
        % sampling at lower resolutions than original voxel dimensions)
    	slice = imresize(slice,[no_cols, no_slices]);
    end
    
end

end