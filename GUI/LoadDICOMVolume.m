function Image = LoadDICOMVolume(range, folder)
% LoadDICOMVolume loads a 3D DICOM volume
%
% DESCRIPTION: Image = LoadDICOMVolume(range, folder)
%       Loads all DICOM files (images and metadata) from a specified folder
%       Files are spatially sorted according to Z value
%       Voxel 3D dimensions are extracted
%       The volume is returned within specified Z limits
%
% INPUTS:
%       range (double vector) - a two element vector, where the first
%       element specifies the lower bound of Z slices to return, and the
%       second element specifies the upper bound of Z slices to return.
%
%       folder (character string) - an optional input argument specifying
%       the folder from which to load the DICOM volume. If no input folder
%       is specified then the user is directed to a standard dialog box to
%       select a folder from which to load the DICOM volume.
% 
% OUTPUTS:
%       Image (1 x 1 structure with two fields) - 
%           .ImageData (double matrix) - of dimensions (number of rows, 
%           number of columns,number of slices) containing the voxel grey 
%           level values
%           .VoxelDimensions (double vector) - a 1 by 3 vector containing
%           the (y,x,z) voxel dimensions in mm, respectively
%
% FUNCTION DEPENDENCIES:
%       None
%
% AUTHOR:
%       Anonymised for MPHYGB24 MATLAB coursework assignment 2017/18


% To load DICOM files from a selected folder into a volume between a
% specified range of z values

% if user has not specified a folder as a function argument
if nargin == 1

    % Select folder to load volume from (dialog box opens in the current 
    % working directory)  
    folder = uigetdir('','Select a folder to load volume from');
    if folder == 0
        % User has clicked the cancel button, display an error message
        error = errordlg('User has clicked cancel.','Load Error','modal');
        % block program execution until user has clicked on modal error box
        disp('DICOM volume NOT loaded.')
        uiwait(error)
        Image=0;
        return; 
    end
end

% To get number of slices in volume

% build a pathname for all DICOM files in selected folder
dicom_pathname = fullfile(folder,'*.dcm');
% build a list of all the DICOM files in the selected folder
dicom_list = dir(dicom_pathname); % structure with 6 fields

% if this structure is empty then the folder does not contain any DICOM
% files
if isempty(dicom_list)
    error = errordlg ...
    ('The selected folder does not contain any DICOM images.', ...
    'Load Error','modal');
    disp('DICOM volume NOT loaded.')
    % block program execution until user has clicked on modal error box
    uiwait(error)
    Image=0;
    return; 
end

% assign number of slices
no_slices = length(dicom_list);

% error checking for Z range values

if range(1) < 1 || range(1) > range(2) || range(2) > no_slices
    error = errordlg({'Enter a valid range of slices to load.'; ...
        ['The volume contains ', num2str(no_slices),' slices.']}, ...
        'Z range error','modal');
        disp('DICOM volume NOT loaded.')
    % block program execution until user has clicked on modal error box
    uiwait(error)
    Image=0;
    return;
end

% diplay message to user at this point
disp('Loading DICOM volume ...')

% To loop through all slices of volume
% Open a waitbar dialog box
h = waitbar(0,'Loading DICOM Volume...');

for i = 1:no_slices
    
    %update waitbar dialog box
    waitbar(i / no_slices)
   
    % For each slice build a full filename  
    fullfilename = fullfile(folder, dicom_list(i).name);

    % load in DICOM metadata using dicominfo
    info = dicominfo (fullfilename);  % structure with 62 fields
    % load in image using dicomread
    slice = dicomread (fullfilename); % matrix of type int16
  
    % for first slice only
    
    if i == 1
        
        % To get width and height of volume
        
        no_cols = info.Width; %number of columns in image
        no_rows = info.Height; %number of rows in image
        % Initialise a structure called Image
        % Initialise field matrix to hold pixel grey level values
        Image.ImageData = zeros(no_rows,no_cols,no_slices);
        % To initialise a field vector to hold the voxel dimensions (mm) in 
        Image.VoxelDimensions = zeros(1,3);
        % To extract the in-plane voxel dimensions
        Image.VoxelDimensions(1:2) = info.PixelSpacing(1:2); 
        % To initialise a vector to hold Z value in for each slice
        Z = zeros(no_slices, 1);

    end
    
    % To define vectors u, v and p for each slice 
    u = info.ImageOrientationPatient(1:3); % 3 x 1 vector of type double
    v = info.ImageOrientationPatient(4:6); % 3 x 1 vector of type double
    p = info.ImagePositionPatient; % 3 x 1 vector of type double
    % To define the position Z (in mm) in the volume for individual slice
    % Assign to Z vector
    Z(i) = dot(p,cross(u,v)); % scalar double
        
    % Convert slice to double for greater precision and assign to ImageData
    Image.ImageData(:,:,i) = double(slice); 

end

%close waitbar dialog box
close(h)

% Then sort all slices by Z number

% First sort Z and keep a record of sorts in vector I
[Z,I] = sort(Z);

% Then the orthogonal voxel dimension is the difference between any two 
% consecutive Z values, choose first two here for simplicity
Image.VoxelDimensions(3) = Z(2) - Z(1);

% Then sort the slices in ImageData using I as a logical index
Image.ImageData = Image.ImageData(:,:,I);

% Then curtail the number of slices in ImageData dependent upon specified
% range of slices

Image.ImageData = Image.ImageData(:,:,range(1):range(2));

%Display message to user at this point. 
disp('DICOM volume loaded.')

end

