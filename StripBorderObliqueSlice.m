function output_slice = StripBorderObliqueSlice(input_slice, stencil_2D)
% Used before 2D imagesc plotting, but not before 3D slice/slab plotting

% DESCRIPTION: output_slice = StripBorderObliqueSlice(input_slice,low_value)
    % Strips the low value border from a 2D slice. This low_value border
    % is produced from rotation of a slice to be outside original
    % volume. Multiplies input_slice by stencil_2D, then strips complete
    % NaN rows and columns.

% INPUTS:
    % input_slice (double matrix) - 2D image slice
    % stencil_2D (double matrix)  - 2D image which is effectively a record
    % of if a voxel has been rotated out of the volume or not. 

% OUTPUTS:
    % output_slice (double matrix) - 2D image slice, with any rows/columns
    % that have been completely rotated out the volume removed. 

% Function dependencies: None

% Author: Anonymised for MPHYGB24 MATLAB coursework assignment 2017/18

% Resize input_slice to be the same size as stencil_2D
size_stencil_2D = size(stencil_2D);
input_slice = imresize(input_slice, [size_stencil_2D(1) size_stencil_2D(2)]);
% Multiply input_slice by stencil_2D to retain values inside volume
output_slice = stencil_2D.*input_slice;
% Remove rows and columns which contain all NaNs from slice by converting 
% to a NULL row/column
output_slice(:,~any(output_slice,1)) = [];  % rows, first dimension
output_slice(~any(output_slice,2),:) = [];  % columns, second dimension
% NB. Any remaining NaNs will scale as lowest value colour on imagesc plots

% ensure output_slice is 2D by removing singleton dimensions if they
% exist
output_slice = squeeze(output_slice);
end

