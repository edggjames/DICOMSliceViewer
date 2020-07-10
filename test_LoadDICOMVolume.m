clearvars
close all 
clc

% MPHYGB24 - MATLAB coursework assignment 2017/18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Script for Task 1
% Tests function LoadDICOMVolume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) To load DICOM volume from a folder specified by user using a standard
% dialog box, with specific range of Z values:

% NB can test error messages here if user clicks cancel or selects a folder
% which contains no DICOM files

z_limits(1) = 1;
z_limits(2) = 181;

Image = LoadDICOMVolume(z_limits);

% Quick check that Image structure has been loaded correctly

if isequal(Image,0)
    % Error message from LoadDICOMVOlume is shown as user has clicked
    % cancel or specified a folder which exists but which contains no DICOM
    % files
else
    % Display fields of Image structure
    disp(Image) 
    % Display all Z slices of volume 
    range = z_limits(2) - z_limits(1) + 1;
    for i = 1:range
        imshow(Image.ImageData(:,:,i),[])
        title(['Showing slice ', num2str(i - z_limits(1) + 1),' of slices ', ...
            num2str(z_limits(1)),' to ' num2str(z_limits(2))])
        pause(0.01)
    end
end

% 2) To load DICOM volume from a pre-specified folder using optional input 
% argument:

% Load name of folder on my PC (which is in current working directory)
folder = fullfile(pwd,'\Pancreas-03130\');
% Please put your folder file path here if different from the above. 

Image = LoadDICOMVolume(z_limits, folder);

% Quick check that Image structure has been loaded correctly

% Display fields of Image structure
disp(Image) 

% Display all Z slices of volume 
figure
range = z_limits(2) - z_limits(1) + 1;
for i = 1:range
    imshow(Image.ImageData(:,:,i),[])
    title(['Showing slice ', num2str(i - z_limits(1) + 1),' of slices ', ...
            num2str(z_limits(1)),' to ' num2str(z_limits(2))])
    pause(0.01)
end

% 3) To test error message for when folder containing no DICOM files is
% selected

% load name of folder containing no DICOM files

LoadDICOMVolume(z_limits, 'folder_wrong');

% 4) To test error message when start slice OF less than 1 is specified
z_limits(1) = 0;

LoadDICOMVolume(z_limits, folder);

% 5) To test error message when end slice of more than upper bound is 
% specified
z_limits(1) = 1;
z_limits(2) = 182;

LoadDICOMVolume(z_limits, folder);

% 6) To test error message when end slice less than start slice is
% specified
z_limits(1) = 10;
z_limits(2) = 9;

LoadDICOMVolume(z_limits, folder);

% 7) To check that one slice can be loaded
z_limits(1) = 10;
z_limits(2) = 10;

% Display structure and slice to test
Image = LoadDICOMVolume(z_limits, folder);
disp(Image) 
figure
imshow(Image.ImageData,[])
title(['Showing slice ', num2str(z_limits(1)), ' only'])