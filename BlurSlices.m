function slab_blurred = BlurSlices(slab, orientation, sigmax, sigmay, sigmaz, dimension)
% BlurSlices performs 1D or 2D Gaussian blurring on one or more slices of a 3D 
% volume, blurring each slice separately 
%
% DESCRIPTION: slab_blurred = BlurSlices(slab, orientation, sigmax, sigmay, sigmaz)
%       Performs 1D or 2D Gaussian blurring on each slice in input slab in a
%       specified orientation
%       Any one or two of sigmax, sigmay or sigmaz are specified for blurring
%       Slices are returned combined into a slab
%
% INPUTS:
%       slab (double 3D matrix) - holds n 2D m by p slices in a n by m by p
%       matrix 
%           
%       orientation (character string) - determines slice plane orientation
%           'X-Y' - XY plane, orthogonal to Z axis (if chosen then sigmaz
%           should be input as '').
%           'Y-Z' - YZ plane, orthogonal to X axis (if chosen then sigmax
%           should be input as '').
%           'X-Z' - XZ plane, orthogonal to Y axis (if chosen then sigmay
%           should be input as '').
%
%       sigmax (double scalar) - determines the width of the Gaussian 
%       kernel in the x-direction (only used if orientation is set to 'X-Y',
%       or 'X-Z', and x-direction is being blurred, otherwise set to '').
%
%       sigmay (double scalar) - determines the width of the Gaussian 
%       kernel in the y-direction (only used if orientation is set to 'X-Y',
%       or 'Y-Z', and y-direction is being blurred, otherwise set to '').        
%
%       sigmaz (double scalar) - determines the width of the Gaussian 
%       kernel in the z-direction (only used if orientation is set to 'Y-Z',
%       or 'X-Z', and z-direction is being blurred, otherwise set to '').
%
%       dimension (character string) - determines the dimension/s to blur,
%       can be set to either 'X', 'Y', 'Z' or 'both'.
%
%OUTPUTS:
%       slab_blurred (double 3D matrix) - n 2D m by p blurrred slices in a 
%       n by m by p matrix.
%
% FUNCTION DEPENDENCIES:
%       GaussianImfilterSep - computes 2 x 1D convolutions of an image with
%       a Gaussian Kernel, making use of the separability of the Gaussian
%       Kernel
%
% AUTHOR:
%       Anonymised for MPHYGB24 MATLAB coursework assignment 2017/18

% Firstly calculate dimensions of input slab
dim_slab = size(slab);
% Assign first two dimensions of blurred slab 
no_rows = dim_slab(1);
no_cols = dim_slab(2);

% For X-Y slice
if strcmp(orientation,'X-Y') == 1
    
    % assign third dimension of blurred slab
    % if only 1 slice is input
    if ismatrix(slab) == 1
        no_slices = 1;
    % otherwise more than one slice has been input
    else
    no_slices = dim_slab(3);
    end
    
    % preallocate matrix of zeros to hold blurred slab in
    slab_blurred = zeros (no_rows, no_cols, no_slices);
    
    % allocate the dimensions to blur 
    if strcmp(dimension,'X') == 1
        dimension = '1';
    elseif strcmp(dimension,'Y') == 1
        dimension = '2';  
    elseif strcmp(dimension,'both') == 1
        dimension = 'both';
    end
        
    %loop through slices
    for i = 1:no_slices
    % blur each slice and update slab_blurred matrix accordingly
    slice = squeeze(slab(:,:,i)); % reduce to 2D slice
    % NB note order of sigmax and sigmay in function argument
    slab_blurred(:,:,i) = GaussianImfilterSep(slice, dimension, sigmax, sigmay);
    end

% For Y-Z slice
elseif strcmp(orientation,'Y-Z') == 1
    
    % assign third dimension of blurred slab
    no_slices = dim_slab(3);
    
    % preallocate matrix of zeros to hold blurred slab in
    slab_blurred = zeros (no_rows, no_cols, no_slices);
    
    % allocate the dimensions to blur 
    if strcmp(dimension,'Y') == 1
        dimension = '2';
    elseif strcmp(dimension,'Z') == 1
        dimension = '1';  
    elseif strcmp(dimension,'both') == 1
        dimension = 'both';
    end
    
    %loop through columns
    for i = 1:no_cols
    % blur each slice and update slab_blurred matrix accordingly
    slice = squeeze(slab(:,i,:)); % reduce to 2D slice
    % NB note order of sigmaz and sigmay in function argument
    slab_blurred(:,i,:) = GaussianImfilterSep(slice, dimension, sigmaz, sigmay);
    end
    
% For X-Z slice    
elseif strcmp(orientation,'X-Z') == 1
    
    % assign third dimension of blurred slab
    no_slices = dim_slab(3);
    
    % preallocate matrix of zeros to hold blurred slab in
    slab_blurred = zeros (no_rows, no_cols, no_slices);
    
    % allocate the dimensions to blur 
    if strcmp(dimension,'X') == 1
        dimension = '2';
    elseif strcmp(dimension,'Z') == 1
        dimension = '1';  
    elseif strcmp(dimension,'both') == 1
        dimension = 'both';
    end
    
    % loop through rows
    for i = 1:no_rows
        % blur each slice and update slab_blurred matrix accordingly
        slice = squeeze(slab(i,:,:)); % reduce to 2D slice
        % NB note order of sigmaz and sigmax in function argument
        slab_blurred(i,:,:) = GaussianImfilterSep(slice, dimension, sigmaz, sigmax);
    end
    
end


end
