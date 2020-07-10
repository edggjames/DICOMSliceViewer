function  image_filtered = GaussianImfilterSep(image, dimension, sigma_1, sigma_2)
% Used by parent function BlurSlices
% GaussianImfilterSep - computes one or two x 1D convolutions of an image with a
% Gaussian Kernel, making use of the separability of the Gaussian Kernel
%
% DESCRIPTION: image_filtered = GaussianImfilterSep(image, sigma_1, sigma_2)
%       Computes a one or two x 1D central convolutions of an image with one 
%       or two 1D Gaussian distributions - this is effectively image blurring/smoothing.
%
% INPUTS:
%       image (double matrix) - image/slice to be filtered/blurred 
%
%       dimension (character string) - the dimension/s to be blurred, can
%       be either set to '1', '2', or 'both' 
%
%       sigma_1 (double scalar)  - Standard deviation of first dimension direction
%       Gaussian Kernel (determines the 'width' of the Gaussian in this 
%       direction)
%       
%       sigma_2 (double scalar)  - Standard deviation of second dimension direction 
%       Gaussian kernel (determines the 'width' of the Gaussian in this 
%       direction)
%
% OUTPUTS:
%       image_filtered (double matrix) - filtered image/slice after convolution
%       with one or two x 1D Gaussian kernels (size of which depends on sigma_1 and
%       sigma_2).
%
% FUNCTION DEPENDENCIES:
%       Dependent upon function 'CompGaussian' to form the Gaussian 1D
%       kernels
%       Also dependent upon 'imfilter.n' - image filtering - part of 
%       MATLAB Image Processing Toolbox)
% 
% AUTHOR:
%        Anonymised for MPHYGB24 MATLAB coursework assignment 2017/18
    
% To filter only in the first dimension
if strcmp (dimension,'1') == 1
    %First define size of Gaussian kernel in first dimension direction using same 
    %statement as inbuilt MATLAB function 'imgaussfilt' 
    N_1 = 2*ceil(2*sigma_1)+1 ; %must be an odd number
    %Generate a distribution of these values around zero
    first = -floor(N_1/2) : floor(N_1/2);

    %Then form 1D Gaussian kernel based on first and sigma_1
    G_1 = CompGaussian(first, sigma_1);  % 1 by N kernel

    %Filter the image using a central 1D convolution of the input image
    %with the 1D Gaussian kernel in the first dimension
    image_filtered = imfilter(image,G_1,'conv','replicate');
    
% To filter only in the second dimension
elseif strcmp (dimension,'2') == 1
    %First define size of Gaussian kernel in second dimension using same 
    %statement as inbuilt MATLAB function 'imgaussfilt' 
    N_2 = 2*ceil(2*sigma_2)+1 ; %must be an odd number
    
    %Generate a distribution of these values around zero
    second = -floor(N_2/2) : floor(N_2/2);

    %Then form 2D Gaussian kernel based on second and sigma_2
    G_2 = CompGaussian(second, sigma_2); % 1 by N kernel
    G_2 = G_2'; % tranpose to a N by 1 kernel
    
    %Filter the image using a central 1D convolution of the input image
    %with the 1D Gaussian kernel in the second dimension
    image_filtered = imfilter(image,G_2,'conv','replicate');

% To filter in both dimensions    
elseif strcmp (dimension,'both') == 1
    
    %First define size of Gaussian kernel in each direction using same 
    %statement as inbuilt MATLAB function 'imgaussfilt' 
    N_1 = 2*ceil(2*sigma_1)+1 ; %must be an odd number
    N_2 = 2*ceil(2*sigma_2)+1 ; %must be an odd number
    %Generate a distribution of these values around zero
    first = -floor(N_1/2) : floor(N_1/2);
    second = -floor(N_2/2) : floor(N_2/2);

    %Then form 2D Gaussian kernel based on first, second, sigma_1 and sigma_2
    G_1 = CompGaussian(first, sigma_1);  % 1 by N kernel
    G_2 = CompGaussian(second, sigma_2); % 1 by N kernel
    G_2 = G_2'; % tranpose to a N by 1 kernel
    
    %Then filter image in both dimensions
    image_filtered = imfilter(image,G_1,'conv','replicate');
    image_filtered = imfilter(image_filtered,G_2,'conv','replicate');
    
end 

end

