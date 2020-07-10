function [output_slice, rotated_vox_dim, stencil_2D] = ComputeObliqueSlice(Image, ...
  orientation, position, res_1, res_2, method, alpha, beta, gamma, thickness, central_position)
% ComputeObliqueSlice computes an oblique slice from a 3D DICOM
% volume. Uses Gaussian blurring of slices when resampling at a lower
% in-plane resolution than rotated original image
%
% DESCRIPTION: [output_slice, rotated_vox_dim] = ComputeObliqueSlice(Image, ...
%     orientation, position, res_1, res_2, method, alpha, beta, gamma)
%       Calculates volume dimensions in terms of both voxels and mm
%       Extracts voxel dimensions
%       2D oblique slice is calculated on specified plane orientation, at
%       specified position, according to input 3D rotation angles,
%       using 3D interpolation (performed using user specified 
%       interpolation method).
%       Gaussian blurring of slices is performed to avoid aliasing when sampling
%       at lower in-plane resolutions than original image in any direction.
%       This removes higher spatial frequencies and avoids pixellation and
%       other image artefacts. 
%       2D slice has in-plane pixel resolutions as specified in function
%       inputs.
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
%       position (double scalar) - determines position of reference slice in mm along
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
%       This is with reference to the rotated voxel dimension in this
%       function.
%
%       res2 (double scalar) - defines the resolution in mm in the output
%       slice of the second in-plane dimension, i.e y if XY plane is
%       selected, z if YZ plane is selected, or z if XZ plane selected.  NB
%       This is with reference to the rotated voxel dimension in this
%       function.
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
%       thickness (character string) - if set to 'slice' then interpolation
%       is performed within the function. If set to 'slab' then only the
%       intensities of the nearest oblique slice are returned, and interpolation
%       will occur within the invoking function. 
%
%       central_position (double scalar) - the position of central slice of
%       slab specified in mm. Use this input if calculating a slab, to
%       specify the centre of rotation of each slice in the slab. This is
%       an optional input argument. 
%
%       NB all three rotation angles are counterclockwise in a RHS
%       cartesian coordinate system, or clockwise in a LHS cartesian
%       coordinate system.
%
%OUTPUTS:
%       output_slice (double matrix) - oblique slice from volume, according
%           to function input. 
%
%       rotated_vox_dim (double vector) - a 3 x 1 vector of containing the
%           rotated voxel dimensions in the form [dy dx dz].
%
%       stencil_2D (double matrix) - a matrix of dimensions mathcing the
%           in-plane dimensions of the unrotated slice. This keeps a record
%           of if a voxel has been rotated out of the volume. If so then a
%           'NaN' is recorded, if not then a '1' is recorded.
%
% FUNCTION DEPENDENCIES:
%       'BlurSlices' - performs 1D or 2D Gaussian blurring on one or more
%       slices of a 3D volume, blurring each slice separately 
%       'ComputeOrthogonalSliceBlur' - computes an orthogonal slice if all
%       three rotation parameters are set to zero
%
% AUTHOR:
%       Anonymised for MPHYGB24 MATLAB coursework assignment 2017/18

% To assign 3D image intensities to a matrix variable
vol = Image.ImageData;
% calculate minimum value in vol
low_value = min(min(min(vol)));
% To extract voxel dimensions in mm (used to scale axes in plots)
vox_dim = Image.VoxelDimensions; % [dy dx dz]
% To extract image dimensions in terms of number of voxels
image_dim = size(Image.ImageData); % [rows cols slices]
    
% To get size of image as 3D vector in mm in form of (y,x,z) i.e.
% (image height, image width, image length)
image_size = (image_dim - [1 1 1]).*vox_dim;
% Assume first slice is centred around 0 mm in respective direction
% Then assign as follows
image_height = image_size(1);
image_width = image_size(2);
image_length = image_size(3);

% Calculate the three RHS counterclockwise / LHS clockwise rotation 
% matrices

% For the first dimension (y)
% -180 degrees <= beta <= +180 degrees
Ry = [ 1, 0,          0
       0, cosd(beta), -sind(beta)
       0, sind(beta), cosd(beta) ];

% For the second dimension (x) 
% -180 degrees <= alpha <= +180 degrees
Rx = [ cosd(alpha),  0, sind(alpha)
       0,            1, 0
       -sind(alpha), 0, cosd(alpha) ];

% For the third dimension (z)
% -180 degrees <= gamma <= +180 degrees
Rz = [ cosd(gamma), -sind(gamma), 0
       sind(gamma), cosd(gamma),  0
       0,           0,            1 ];

% Multiply the three rotation matrices together

R = Ry * Rx * Rz;

% Scale rotated image voxels correctly according to rotation matrix
rotated_vox_dim = (abs(R) * [vox_dim(1); vox_dim(2); vox_dim(3)]);

% To carry out orientation specific tasks:

% For XY slice
if strcmp(orientation,'X-Y') == 1
    
    % Check position is within bounds
    if position < 0 || position > image_length 
        % out of bound measurement has been requested
        error = errordlg({'Reference X-Y slice is outside volume dimensions'; ...
            ['Z Slice Position cannot be less than 0.000 or greater than ' ...
            , sprintf('%.3f',image_length)]},'Interpolation Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
        return
    end
    
    % Convert position to a slice number to index original volume with, 
    % assuming slice 1 is centred around 0 mm
    z_query = (position/vox_dim(3)) + 1;
    
    % initialise a matrix for stencil_2D
    stencil_2D = ones(image_dim(1),image_dim(2))*NaN;
    
    if strcmp(thickness,'slice') == 1
    
    % Assign a slab of 2 or 4 matrices of zeros to hold rotated slices in, which 
    % are each the same size as the padded unrotated slice. These will be 
    % interpolated later
    
    if strcmp(method,'nearest') == 1 || strcmp(method,'linear') == 1
        rotated_slab = ones(image_dim(1),image_dim(2),2)*low_value;
        first_slice = 0;  last_slice = 1;
    elseif strcmp(method,'spline') == 1
        rotated_slab = ones(image_dim(1),image_dim(2),4)*low_value;
        first_slice = -1;  last_slice = 2;
    end
    
    % Multiply by each voxel location by R in turn in 2 for loops for each
    % of the 2 or 4 rotated slices

    for k = first_slice:last_slice
    
        % looping through y values
        for i = 1:image_dim(1)
    
            % looping through x values
            for j = 1:image_dim(2)
        
                % Transform the coordinates so the centre of rotation is at the
                % centre of the slice
                Transformed = [i - (image_dim(1)/2);
                               j - (image_dim(2)/2);
                               0 - k];
        
                % Rotate the coordinates by pre-multiplying by R
                Rotated = R * Transformed;
            
                % Transform back to original coordinate system
                Transformed_rotated = Rotated + [image_dim(1)/2;
                                                 image_dim(2)/2;
                                                 floor(z_query)+k];
           
                % Extract the rotated coordinates (rounded to the nearest 
                % coordinate)                             
                y = round(Transformed_rotated(1));
                x = round(Transformed_rotated(2));
                z = round(Transformed_rotated(3));
        
                % Assign relevant image intensity to rotated slice (if it
                % exists)

                if  y >= 1 && y <= image_dim(1) && ...
                    x >= 1 && x <= image_dim(2) && ...
                    z >= 1 && z <= image_dim(3)
           
                    if strcmp(method,'nearest') == 1 || strcmp(method,'linear') == 1
                        rotated_slab(i,j,k+1) = vol(y,x,z);
                    elseif strcmp(method,'spline') == 1
                        rotated_slab(i,j,k+2) = vol(y,x,z);
                    end
                    if k == 0
                        stencil_2D(i,j) = 1;
                    end
                end
            end
        end
    end      
        
    % To convert resolutions from mm to number of rotated voxel dimensions for
    % 3D interpolation for in-plane dimensions
    x_res_vox = res_1/rotated_vox_dim(2);
    y_res_vox = res_2/rotated_vox_dim(1);
    % Convert these into vectors of query points for interpolation
    x_query = 1:x_res_vox:image_dim(2);
    y_query = 1:y_res_vox:image_dim(1);

    % blur the slices if need be here using similar code to
    % ComputeOrthogonalSliceBlur. Form expressions for sigmax, sigmay and the 
    % dimensions to perform Gaussian blurring in

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

     % blur these slices if necessary using function from Task 4
     if x_res_vox >= 1.001 || y_res_vox >= 1.001
        rotated_slab = BlurSlices(rotated_slab,'X-Y',sigmax,sigmay,'',dimension);
     end  

    % Interpolate between these slices
    if strcmp(method,'nearest') == 1
        % define a z_query point for interpolation (i.e. first or second slice)
        if floor(z_query) == round(z_query)
            z_query_interp = 1;
        else
            z_query_interp = 2;
        end
        % set up a 3D grid of query points for interpolation
        [x_query,y_query,z_query_interp] = meshgrid(x_query,y_query,z_query_interp);
        % perform 3D interpolation
        output_slice = interp3(rotated_slab,x_query,y_query,z_query_interp,method);
    
    elseif strcmp(method,'linear') == 1
        % define a z_query point for interpolation (somewhere between first and
        % second slice)
        z_query_interp = z_query - floor(z_query) + 1;
        % set up a 3D grid of query points for interpolation
        [x_query,y_query,z_query_interp] = meshgrid(x_query,y_query,z_query_interp);
        % perform 3D interpolation
        output_slice = interp3(rotated_slab,x_query,y_query,z_query_interp,method);

    elseif strcmp(method,'spline') == 1
        % define a z_query point for interpolation (somewhere between second
        % and third slice)
        z_query_interp = z_query - floor(z_query) + 2;
        % set up a 3D grid of query points for interpolation
        [x_query,y_query,z_query_interp] = meshgrid(x_query,y_query,z_query_interp);
        % perform 3D interpolation
        output_slice = interp3(rotated_slab,x_query,y_query,z_query_interp,method);
    end
    
    elseif strcmp(thickness,'slab') == 1
    % Assign a matrix to hold returned slice in
    output_slice = ones(image_dim(1),image_dim(2))*low_value;    
    % Calculate the difference between this slice and the centre of slab
    % (so that the slice is rotated around the centre of the slab, not the
    % centre of each slice in the slab!)
    delta_slab = central_position - position;
    
    % looping through y values
        for i = 1:image_dim(1)
    
            % looping through x values
            for j = 1:image_dim(2)
        
                % Transform the coordinates so the centre of rotation is at the
                % centre of the slab
                Transformed = [i - (image_dim(1)/2);
                               j - (image_dim(2)/2);
                               0 - delta_slab];
        
                % Rotate the coordinates by pre-multiplying by R
                Rotated = R * Transformed;
            
                % Transform back to original coordinate system
                Transformed_rotated = Rotated + [image_dim(1)/2;
                                                 image_dim(2)/2;
                                                 z_query + delta_slab];
           
                % Extract the rotated coordinates (rounded to the nearest 
                % coordinate)                             
                y = round(Transformed_rotated(1));
                x = round(Transformed_rotated(2));
                z = round(Transformed_rotated(3));
        
                % Assign relevant image intensity to rotated slice (if it
                % exists)

                if  y >= 1 && y <= image_dim(1) && ...
                    x >= 1 && x <= image_dim(2) && ...
                    z >= 1 && z <= image_dim(3)
                        output_slice(i,j) = vol(y,x,z);
                        stencil_2D(i,j) = 1;
                end
            end
        end
    end
    
% For YZ slice
elseif strcmp(orientation,'Y-Z') == 1

    % Check position is within bounds
    if position < 0 || position > image_width 
        % out of bound measurement has been requested
        error = errordlg({'Reference Y-Z slice is outside volume dimensions'; ...
            ['X Slice Position cannot be less than 0.000 or greater than ' ...
            , sprintf('%.3f',image_width)]},'Interpolation Error','modal');
        %block program execution until user has clicked on modal error box
        uiwait(error)
        return
    end
    
    % Convert position to a slice number to index original volume with, 
    x_query = (position/vox_dim(2)) + 1;
    % assuming last slice is centred around 0 mm
    x_query = image_dim(2) - x_query + 1;
    
    % initialise a matrix for stencil_2D
    stencil_2D = ones(image_dim(1),image_dim(3))*NaN;
    
    if strcmp(thickness,'slice') == 1
    
    % Assign a slab of 2 or 4 matrices of zeros to hold rotated slices in, which 
    % are each the same size as the padded unrotated slice. These will be 
    % interpolated later
    if strcmp(method,'nearest') == 1 || strcmp(method,'linear') == 1
        rotated_slab = ones(image_dim(1),2,image_dim(3))*low_value;
        first_slice = 0;  last_slice = 1;
    elseif strcmp(method,'spline') == 1
        rotated_slab = ones(image_dim(1),4,image_dim(3))*low_value;
        first_slice = -1;  last_slice = 2;
    end
    
    % Multiply by each voxel location by R in turn in 2 for loops for each
    % of the 2 or 4 rotated slices     

    for k = first_slice:last_slice
    
        % looping through y values
        for i = 1:image_dim(1)
    
            % looping through z values
            for j = 1:image_dim(3)
        
                % Transform the coordinates so the centre of rotation is at the
                % centre of the slice
                Transformed = [i - (image_dim(1)/2);
                               0 - k
                               j - (image_dim(3)/2)];
        
                % Rotate the coordinates by pre-multiplying by R
                Rotated = R * Transformed;

                % Transform coordinates back to that of original volume

                Transformed_rotated = Rotated + [image_dim(1)/2;                                                 
                                                 floor(x_query)+k
                                                 image_dim(3)/2;];
           
                % Extract the rotated coordinates (rounded to the nearest 
                % coordinate)                             
                y = round(Transformed_rotated(1));
                x = round(Transformed_rotated(2));
                z = round(Transformed_rotated(3));
        
                % Assign relevant image intensity to rotated slice (if it
                % exists)

                if  y >= 1 && y <= image_dim(1) && ...
                    x >= 1 && x <= image_dim(2) && ...
                    z >= 1 && z <= image_dim(3)
                        if strcmp(method,'nearest') == 1 || strcmp(method,'linear') == 1
                            rotated_slab(i,k+1,j) = vol(y,x,z);
                        elseif strcmp(method,'spline') == 1
                            rotated_slab(i,k+2,j) = vol(y,x,z);
                        end
                        if k == 0
                            stencil_2D(i,j) = 1;
                        end
                end
            end
        end
    end
    
    % To convert resolutions from mm to number of rotated voxel dimensions for
    % 3D interpolation for in-plane dimensions
    y_res_vox = res_1/rotated_vox_dim(1);
    z_res_vox = res_2/rotated_vox_dim(3);
    % Convert these into vectors of query points for interpolation
    y_query = 1:y_res_vox:image_dim(1);
    z_query = 1:z_res_vox:image_dim(3);

    % blur the slices if need be here using similar code to
    % ComputeOrthogonalSliceBlur. Form expressions for sigmay, sigmaz and the 
    % dimensions to perform Gaussian blurring in

    sigmay = ((sqrt(2*log(2)))/pi)*y_res_vox;
    sigmaz = ((sqrt(2*log(2)))/pi)*z_res_vox;

    if y_res_vox >= 1.001 && z_res_vox >= 1.001
        dimension = 'both';
    elseif y_res_vox >= 1.001 
        dimension = 'Y';
    elseif z_res_vox >= 1.001 
        dimension = 'Z';
    end

    % Blur these slices if necessary using function from Task 4
    if y_res_vox >= 1.001 || z_res_vox > 1.001
        rotated_slab = BlurSlices(rotated_slab,'Y-Z','',sigmay,sigmaz,dimension);
    end 

    % Interpolate between these slices
    if strcmp(method,'nearest') == 1
        % define an x_query point for interpolation (i.e. first or second slice)
        if floor(x_query) == round(x_query)
            x_query_interp = 1;
        else
            x_query_interp = 2;
        end
        % set up a 3D grid of query points for interpolation
        [x_query_interp,y_query,z_query] = meshgrid(x_query_interp,y_query,z_query);
        % perform 3D interpolation
        output_slice = interp3(rotated_slab,x_query_interp,y_query,z_query,method);
    
    elseif strcmp(method,'linear') == 1
        % define an x_query point for interpolation (somewhere between first and
        % second slice)
        x_query_interp = x_query - floor(x_query) + 1;
        % set up a 3D grid of query points for interpolation
        [x_query_interp,y_query,z_query] = meshgrid(x_query_interp,y_query,z_query);
        % perform 3D interpolation
        output_slice = interp3(rotated_slab,x_query_interp,y_query,z_query,method);

    elseif strcmp(method,'spline') == 1
        % define an x_query point for interpolation (somewhere between second
        % and third slice)
        x_query_interp = x_query - floor(x_query) + 2;
        % set up a 3D grid of query points for interpolation
        [x_query_interp,y_query,z_query] = meshgrid(x_query_interp,y_query,z_query);
        % perform 3D interpolation
        output_slice = interp3(rotated_slab,x_query_interp,y_query,z_query,method);
    end
    
    elseif strcmp(thickness,'slab') == 1
        
    % Also assign a matrix to hold position indices of rotated voxels in
    output_slice = ones(image_dim(1),image_dim(3))*low_value;
    % Calculate the difference between this slice and the centre of slab
    % (so that the slice is rotated around the centre of the slab, not the
    % centre of each slice in the slab!)
    delta_slab = central_position - position;
        
    % looping through y values
    for i = 1:image_dim(1)
    
        % looping through z values
        for j = 1:image_dim(3)
        
            % Transform the coordinates so the centre of rotation is at the
            % centre of the slab
            Transformed = [i - (image_dim(1)/2)
                           0 - delta_slab
                           j - (image_dim(3)/2)];
        
             % Rotate the coordinates by pre-multiplying by R
             Rotated = R * Transformed;

             % Transform coordinates back to that of original volume

             Transformed_rotated = Rotated + [image_dim(1)/2                                                 
                                              x_query + delta_slab
                                              image_dim(3)/2];
           
             % Extract the rotated coordinates (rounded to the nearest 
             % coordinate)                             
             y = round(Transformed_rotated(1));
             x = round(Transformed_rotated(2));
             z = round(Transformed_rotated(3));
        
             % Assign relevant image intensity to rotated slice (if it
             % exists)

             if  y >= 1 && y <= image_dim(1) && ...
                 x >= 1 && x <= image_dim(2) && ...
                 z >= 1 && z <= image_dim(3)
                    output_slice(i,j) = vol(y,x,z);
                    stencil_2D(i,j) = 1;
             end
        end
    end    
    end
   
% For XZ slice
elseif strcmp(orientation,'X-Z') == 1

    % Check position is within bounds
    if position < 0 || position > image_height
        % out of bound measurement has been requested
        error = errordlg({'Reference X-Z slice is outside volume dimensions'; ...
            ['Y Slice Position cannot be less than 0.000 or greater than ' ...
            , sprintf('%.3f',image_height)]},'Interpolation Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
        return
    end
    
    % Convert position to a slice number to index original volume with, 
    y_query = (position/vox_dim(1)) + 1;
    % assuming last slice is centred around 0 mm
    y_query = image_dim(1) - y_query + 1;
    
    % initialise a matrix for stencil_2D
    stencil_2D = ones(image_dim(2),image_dim(3))*NaN;
    
    if strcmp(thickness,'slice') == 1
    
    % Assign a slab of 2 or 4 matrices of zeros to hold rotated slices in, which 
    % are each the same size as the padded unrotated slice. These will be 
    % interpolated later
    if strcmp(method,'nearest') == 1 || strcmp(method,'linear') == 1
        rotated_slab = ones(2,image_dim(2),image_dim(3))*low_value;
        first_slice = 0;  last_slice = 1;
    elseif strcmp(method,'spline') == 1
        rotated_slab = ones(4,image_dim(2),image_dim(3))*low_value;
        first_slice = -1; last_slice = 2;
    end
    
    % Multiply by each voxel location by R in turn in 2 for loops for each
    % of the 2 or 4 rotated slices
    
    for k = first_slice:last_slice

        % looping through x values
        for i = 1:image_dim(2)
    
            % looping through z values
            for j = 1:image_dim(3)
        
                % Transform the coordinates so the centre of rotation is at the
                % centre of the slice
                Transformed = [0 - k
                               i - (image_dim(2)/2)
                               j - (image_dim(3)/2)];
        
                % Rotate the coordinates by pre-multiplying by R
                Rotated = R * Transformed;

                % Transform coordinates back to that of original volume

                Transformed_rotated = Rotated + [floor(y_query)+k
                                                 image_dim(2)/2
                                                 image_dim(3)/2];
           
                % Extract the rotated coordinates (rounded to the nearest 
                % coordinate)                             
                y = round(Transformed_rotated(1));
                x = round(Transformed_rotated(2));
                z = round(Transformed_rotated(3));
        
                % Assign relevant image intensity to rotated slice (if it
                % exists)

                if  y >= 1 && y <= image_dim(1) && ...
                    x >= 1 && x <= image_dim(2) && ...
                    z >= 1 && z <= image_dim(3)
           
                    if strcmp(method,'nearest') == 1 || strcmp(method,'linear') == 1
                        rotated_slab(k+1,i,j) = vol(y,x,z);
                    elseif strcmp(method,'spline') == 1
                        rotated_slab(k+2,i,j) = vol(y,x,z);
                    end
                    if k == 0
                        stencil_2D(i,j) = 1;
                    end  
                end
            end
        end
    end
 
    % To convert resolutions from mm to number of rotated voxel dimensions for
    % 3D interpolation for in-plane dimensions
    x_res_vox = res_1/rotated_vox_dim(2);
    z_res_vox = res_2/rotated_vox_dim(3);
    % Convert these into vectors of query points for interpolation
    x_query = 1:x_res_vox:image_dim(2);
    z_query = 1:z_res_vox:image_dim(3);

    % blur the slices if need be here using similar code to
    % ComputeOrthogonalSliceBlur. Form expressions for sigmax, sigmaz and the 
    % dimensions to perform Gaussian blurring in

    sigmax = ((sqrt(2*log(2)))/pi)*x_res_vox;
    sigmaz = ((sqrt(2*log(2)))/pi)*z_res_vox;

    if x_res_vox >= 1.001 && z_res_vox >= 1.001
        dimension = 'both';
    elseif x_res_vox >= 1.001 
        dimension = 'X';
    elseif z_res_vox >= 1.001 
        dimension = 'Z';
    end

    % Blur these slices if necessary using function from Task 4
    if x_res_vox >= 1.001 || z_res_vox >= 1.001
        rotated_slab = BlurSlices(rotated_slab,'X-Z',sigmax,'',sigmaz,dimension);
    end 

    % Interpolate between these slices
    if strcmp(method,'nearest') == 1
        % define a y_query point for interpolation (i.e. first or second slice)
        if floor(y_query) == round(y_query)
            y_query_interp = 1;
        else
            y_query_interp = 2;
        end
        % set up a 3D grid of query points for interpolation
        [x_query,y_query_interp,z_query] = meshgrid(x_query,y_query_interp,z_query);
        % perform 3D interpolation
        output_slice = interp3(rotated_slab,x_query,y_query_interp,z_query,method);
    
    elseif strcmp(method,'linear') == 1
        % define a z_query point for interpolation (somewhere between first and
        % second slice)
        y_query_interp = y_query - floor(y_query) + 1;
        % set up a 3D grid of query points for interpolation
        [x_query,y_query_interp,z_query] = meshgrid(x_query,y_query_interp,z_query);
        % perform 3D interpolation
        output_slice = interp3(rotated_slab,x_query,y_query_interp,z_query,method);

    elseif strcmp(method,'spline') == 1
        % define a z_query point for interpolation (somewhere between second
        % and third slice)
        y_query_interp = y_query - floor(y_query) + 2;
        % set up a 3D grid of query points for interpolation
        [x_query,y_query_interp,z_query] = meshgrid(x_query,y_query_interp,z_query);
        % perform 3D interpolation
        output_slice = interp3(rotated_slab,x_query,y_query_interp,z_query,method);
    end
    
    elseif strcmp(thickness,'slab') == 1 
    
    % Also assign a matrix to hold position indices of rotated voxels in
    output_slice = ones(image_dim(2),image_dim(3))*low_value;  
    % Calculate the difference between this slice and the centre of slab
    % (so that the slice is rotated around the centre of the slab, not the
    % centre of each slice in the slab!)
    delta_slab = central_position - position;
    
    % looping through x values
    for i = 1:image_dim(2)
    
        % looping through z values
        for j = 1:image_dim(3)
        
            % Transform the coordinates so the centre of rotation is at the
            % centre of the slab
            Transformed = [0 - delta_slab
                           i - (image_dim(2)/2)
                           j - (image_dim(3)/2)];
        
            % Rotate the coordinates by pre-multiplying by R
            Rotated = R * Transformed;

            % Transform coordinates back to that of original volume
            Transformed_rotated = Rotated + [y_query + delta_slab
                                             image_dim(2)/2
                                             image_dim(3)/2];
           
            % Extract the rotated coordinates (rounded to the nearest 
            % coordinate)                             
            y = round(Transformed_rotated(1));
            x = round(Transformed_rotated(2));
            z = round(Transformed_rotated(3));
        
            % Assign relevant image intensity to rotated slice (if it
            % exists)

            if  y >= 1 && y <= image_dim(1) && ...
                x >= 1 && x <= image_dim(2) && ...
                z >= 1 && z <= image_dim(3)
                    output_slice(i,j) = vol(y,x,z);    
                    stencil_2D(i,j) = 1;
            end
        end
    end
    end
    
end

% ensure that output slice is a 2D slice with no singleton dimensions
output_slice = squeeze(output_slice);

end

