% MPHYGB24 - MATLAB coursework assignment 2017/18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUI Script for Tasks 3 & 5 & 7 & 10
% To be used in conjunction with SliceViewer.fig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = SliceViewer(varargin)
% SLICEVIEWER MATLAB code for SliceViewer.fig
%      SLICEVIEWER, by itself, creates a new SLICEVIEWER or raises the existing
%      singleton*.
%
%      H = SLICEVIEWER returns the handle to a new SLICEVIEWER or the handle to
%      the existing singleton*.
%
%      SLICEVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SLICEVIEWER.M with the given input arguments.
%
%      SLICEVIEWER('Property','Value',...) creates a new SLICEVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SliceViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SliceViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SliceViewer

% Last Modified by GUIDE v2.5 23-Jan-2018 20:05:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SliceViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @SliceViewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before SliceViewer is made visible.
function SliceViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SliceViewer (see VARARGIN)

% Choose default command line output for SliceViewer
handles.output = hObject;
% Set 'DICOM' to empty (i.e. a record of if a DICOM is loaded or not)
set(handles.Load_Pushbutton,'UserData',[]);
% invoke function to set all values to default
set_default_values(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SliceViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = SliceViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in Plane_Pop_Up_Menu.
function Plane_Pop_Up_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Plane_Pop_Up_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Plane_Pop_Up_Menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Plane_Pop_Up_Menu

% If a volume is loaded, update max, min and step values for both sliders and 
% default values for resolution edit boxes
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    % display a warning dialog box
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
else
    Image = get(handles.Load_Pushbutton,'UserData');
    % 1) To extract voxel dimensions in mm
    vox_dim = Image.VoxelDimensions; % [dy dx dz]
    % 2) To get size of image as 3D vector in form of (y,x,z)
    image_dim = size(Image.ImageData);
    % 3) To get size of image as 3D vector in mm in form of (y,x,z)
    image_size = (image_dim - [1 1 1]).*vox_dim;
end

% update static texts, sliders and image resolutions relating to X, Y and Z 
% when a different plane is selected.
switch get(hObject,'Value')
    case 1
        % orientation = 'X-Y'
        set(handles.Res_1_Text, 'String', 'X (mm):');
        set(handles.Res_2_Text, 'String', 'Y (mm):');
        set(handles.Res_3_Text, 'String', 'Z (mm):');
        set(handles.Slice_Position_Text, 'String', 'Central Slice Position (Z / mm):');
        if empty ~= true
            % DICOM volume is loaded
            position = num2str(image_size(3)/2,'%.3f');
            set(handles.Slice_Position_Edit_Box,'String', position);
            set(handles.Slice_Position_Slider,'Min',0,'Max',image_size(3), ...
                'Value', image_size(3)/2);
            % set minor slider step to 1 slice, and major slider step to 10
            % slices
            set(handles.Slice_Position_Slider,'SliderStep',[1, 10] / (image_dim(3)-1));
            res_1 = num2str (vox_dim(2),'%.3f');
            res_2 = num2str (vox_dim(1),'%.3f');
            res_3 = num2str (vox_dim(3),'%.3f');
            set(handles.Res_1_Edit_Box,'String', res_1);
            set(handles.Res_2_Edit_Box,'String', res_2);
            set(handles.Res_3_Edit_Box,'String', res_3);
            set(handles.slab_thickness_slider,'Min',vox_dim(3),'Max',...
                image_size(3),'Value', vox_dim(3));
            % For slab thickness slider, set minor slider step to add 1 slice to slab,
            % and major slider step to add 10 slices to slab
            set(handles.slab_thickness_slider,'SliderStep',[1, 10] / (image_dim(3)-2));
            set(handles.slab_thickness_edit_box,'String', res_3);
        end
    case 2
        % orientation = 'Y-Z'
        set(handles.Res_1_Text, 'String', 'Y (mm):');
        set(handles.Res_2_Text, 'String', 'Z (mm):');
        set(handles.Res_3_Text, 'String', 'X (mm):');
        set(handles.Slice_Position_Text, 'String', 'Central Slice Position (X / mm):');
        if empty ~= true
            position = num2str(image_size(2)/2,'%.3f');
            set(handles.Slice_Position_Edit_Box,'String', position);
            set(handles.Slice_Position_Slider,'Min',0,'Max',image_size(2), ...
                'Value', image_size(2)/2);
            % set minor slider step to 1 slice, and major slider step to 10
            % slices
            set(handles.Slice_Position_Slider,'SliderStep',[1, 10] / (image_dim(2)-1));
            res_1 = num2str (vox_dim(1),'%.3f');
            res_2 = num2str (vox_dim(3),'%.3f');
            res_3 = num2str (vox_dim(2),'%.3f');
            set(handles.Res_1_Edit_Box,'String', res_1);
            set(handles.Res_2_Edit_Box,'String', res_2);
            set(handles.Res_3_Edit_Box,'String', res_3);
            % For slab thickness slider, set minor slider step to add 1 slice to slab,
            % and major slider step to add 10 slices to slab
            set(handles.slab_thickness_slider,'Min',vox_dim(2),'Max',...
                image_size(2),'Value', vox_dim(2));
            % For slab thickness slider, set minor slider step to add 1 slice to slab,
            % and major slider step to add 10 slices to slab
            set(handles.slab_thickness_slider,'SliderStep',[1, 10] / (image_dim(2)-2));
            set(handles.slab_thickness_edit_box,'String', res_3)
        end
        
    case 3
        % orientation = 'X-Z'
        set(handles.Res_1_Text, 'String', 'X (mm):');
        set(handles.Res_2_Text, 'String', 'Z (mm):');
        set(handles.Res_3_Text, 'String', 'Y (mm):');
        set(handles.Slice_Position_Text, 'String', 'Central Slice Position (Y / mm):');
        if empty ~= true
            position = num2str(image_size(1)/2,'%.3f');
            set(handles.Slice_Position_Edit_Box,'String', position);
            set(handles.Slice_Position_Slider,'Min',0,'Max',image_size(1), ...
                'Value', image_size(1)/2);
            % set minor slider step to 1 slice, and major slider step to 10
            % slices
            set(handles.Slice_Position_Slider,'SliderStep',[1, 10] / (image_dim(1)-1));
            res_1 = num2str (vox_dim(2),'%.3f');
            res_2 = num2str (vox_dim(3),'%.3f');
            res_3 = num2str (vox_dim(1),'%.3f');
            set(handles.Res_1_Edit_Box,'String', res_1);
            set(handles.Res_2_Edit_Box,'String', res_2);
            set(handles.Res_3_Edit_Box,'String', res_3);
            % For slab thickness slider, set minor slider step to add 1 slice to slab,
            % and major slider step to add 10 slices to slab
            set(handles.slab_thickness_slider,'Min',vox_dim(1),'Max',...
                image_size(1),'Value', vox_dim(1));
            set(handles.slab_thickness_slider,'SliderStep',[1, 10] / ((image_dim(1)-2)));
            set(handles.slab_thickness_edit_box,'String', res_3);
        end
end

% also reset all rotational parameters to zero
% alpha edit box and slider
set(handles.alpha_edit_box,'String', '0.0');
set(handles.alpha_slider,'Value', 0);
% beta edit box and slider
set(handles.beta_edit_box,'String', '0.0');
set(handles.beta_slider,'Value', 0);
% gamma edit box and slider
set(handles.gamma_edit_box,'String', '0.0');
set(handles.gamma_slider,'Value', 0);

% call the plot_slab_and_image function to update image display if a volume is loaded
if empty ~= true
    plot_slab_and_image(handles)
end

% --- Executes during object creation, after setting all properties.
function Plane_Pop_Up_Menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Plane_Pop_Up_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Res_1_Edit_Box_Callback(hObject, eventdata, handles)
% hObject    handle to Res_1_Edit_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Res_1_Edit_Box as text
%        str2double(get(hObject,'String')) returns contents of Res_1_Edit_Box as a double

% To show all values to 3dp
res_1_string = get(hObject,'String');
% Convert string to double
res_1_number = str2double(res_1_string);

% check to see if a volume is loaded
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    if res_1_number >= 0.5 && res_1_number <= 20
        % Convert to 3dp string
        res_1_string_3dp = sprintf('%.3f',res_1_number);
        % set value to this updated string
        set(hObject,'String', res_1_string_3dp);
    elseif res_1_number < 0.5
        % resolution value is too low
        set(hObject,'String', '0.500');
    elseif res_1_number > 20
        % resolution value is too high
        set(hObject,'String', '20.000');
    end
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
else
    % Only accept sensible resolution values between 0.5 and 20 mm
    if res_1_number >= 0.5 && res_1_number <= 20
        % Convert to 3dp string
        res_1_string_3dp = sprintf('%.3f',res_1_number);
        % set value to this updated string
        set(hObject,'String', res_1_string_3dp);
    elseif res_1_number < 0.5
        % resolution value is too low
        error = errordlg('Minimum resolution value is 0.500','Resolution Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
        set(hObject,'String', '0.500');
    elseif res_1_number > 20
        % resolution value is too high
        error = errordlg('Maximum resolution value is 20.000','Resolution Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
        set(hObject,'String', '20.000');
    end
end

% --- Executes during object creation, after setting all properties.
function Res_1_Edit_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Res_1_Edit_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function Slice_Position_Slider_Callback(hObject, eventdata, handles)
% hObject    handle to Slice_Position_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% To link slice position slider value to slice position edit box
position = get(hObject,'Value');
% Convert number to string with 3dp precision
position = num2str (position,'%.3f');
% Set edit box to slider value
set(handles.Slice_Position_Edit_Box,'String', position)

% Check if a DICOM volume is loaded
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
end

% --- Executes during object creation, after setting all properties.
function Slice_Position_Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slice_Position_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in Load_Pushbutton.
function Load_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load DICOM volume from a folder specified by user using a standard
% dialog box, with specific range of Z values (here chosen to be all slices
% with prior knowledge of volume z dimensions):
z_limits(1) = 1;
z_limits(2) = 181;

Loaded_Image = LoadDICOMVolume(z_limits);
if isequal(Loaded_Image,0)
    % if user clicks cancel or selects a folder with no DICOM files in then
    % LoadDICOMVolume shows appropriate error message. Loaded_Image
    % structure remainds unaffected. 
else
    % Volume has been loaded. Store the loaded image in UserData field of 
    % handles.Load_Pushbutton
    set(hObject,'UserData',Loaded_Image)
    % set variables to default loaded image values
    set_loaded_image_values(handles)
    % call the plot_slab_and_image function to update image display
    plot_slab_and_image(handles)
end

% --- Executes on button press in Reset_Pushbutton.
function Reset_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Reset_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check if image is loaded
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    % If not then set default values for when image is not loaded
    set_default_values(handles);
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
else
    % If image is loaded then set default values for loaded image
    set_loaded_image_values(handles);
    % update display
    plot_slab_and_image(handles)
end

function Slice_Position_Edit_Box_Callback(hObject, eventdata, handles)
% hObject    handle to Slice_Position_Edit_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Slice_Position_Edit_Box as text
%        str2double(get(hObject,'String')) returns contents of Slice_Position_Edit_Box as a double

% To show all values to 3dp
position_string = get(hObject,'String');
% Convert string to double
position_number = str2double(position_string);

empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
    if position_number >= 0 && position_number <=100
        % Convert to 3dp string
        position_string_3dp = sprintf('%.3f', position_number);
        % set value to this updated string
        set(hObject,'String', position_string_3dp);
        % To link slice position edit box value to slice position slider value
        set(handles.Slice_Position_Slider,'Value', position_number)
    elseif position_number < 0
        set(hObject,'String', '0.000');
        set(handles.Slice_Position_Slider,'Value', 0)
    elseif position_number > 100
        set(hObject,'String', '100.000');
        set(handles.Slice_Position_Slider,'Value', 100)   
    end
else
    % use function to get orthogonal image length
    length = get_image_length(handles);
    if position_number >= 0 && position_number <=length
        % Convert to 3dp string
        position_string_3dp = sprintf('%.3f', position_number);
        % set value to this updated string
        set(hObject,'String', position_string_3dp);
        % To link slice position edit box value to slice position slider value
        set(handles.Slice_Position_Slider,'Value', position_number)
    elseif position_number < 0
        set(hObject,'String', '0.000');
        set(handles.Slice_Position_Slider,'Value', 0)
        error = errordlg('Central slice position cannot have a negative value.', ...
            'Slab Error','modal');
        disp('Oblique slab NOT computed.')
        % block program execution until user has clicked on modal error box
        uiwait(error)
    elseif position_number > length
        position_string_3dp = sprintf('%.3f', length);
        set(hObject,'String', position_string_3dp);
        set(handles.Slice_Position_Slider,'Value', length)
        error = errordlg({'Central slice position exceeds volume dimensions.' ...
            ['Orthogonal image volume dimension is ', sprintf('%.3f',length),' mm.']} ...
            ,'Slab Error','modal');
        disp('Oblique slab NOT computed.')
        % block program execution until user has clicked on modal error box
        uiwait(error)
    end
end

% --- Executes during object creation, after setting all properties.
function Slice_Position_Edit_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slice_Position_Edit_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Res_2_Edit_Box_Callback(hObject, eventdata, handles)
% hObject    handle to Res_2_Edit_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Res_2_Edit_Box as text
%        str2double(get(hObject,'String')) returns contents of Res_2_Edit_Box as a double

% To show all values to 3dp
res_2_string = get(hObject,'String');
% Convert string to double
res_2_number = str2double(res_2_string);

% check to see if a volume is loaded
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    if res_2_number >= 0.5 && res_2_number <= 20
        % Convert to 3dp string
        res_2_string_3dp = sprintf('%.3f',res_2_number);
        % set value to this updated string
        set(hObject,'String', res_2_string_3dp);
    elseif res_2_number < 0.5
        % resolution value is too low
        set(hObject,'String', '0.500');
    elseif res_2_number > 20
        % resolution value is too high
        set(hObject,'String', '20.000');
    end
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
else
    % Only accept sensible resolution values
    if res_2_number >= 0.5 && res_2_number <= 20
        % Convert to 3dp string
        res_2_string_3dp = sprintf('%.3f',res_2_number);
        % set value to this updated string
        set(hObject,'String', res_2_string_3dp);
    elseif res_2_number < 0.5
        % resolution value is too low
        error = errordlg('Minimum resolution value is 0.500','Resolution Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
        set(hObject,'String', '0.500');
    elseif res_2_number > 20
        % resolution value is too high
        error = errordlg('Maximum resolution value is 20.000','Resolution Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
        set(hObject,'String', '20.000');
    end
end

% --- Executes during object creation, after setting all properties.
function Res_2_Edit_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Res_2_Edit_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% plot slab and image function invoked by load image, changing view plane,
% reset, or recalculate slab . NB not invoked by changing projection view
% method alone.
function plot_slab_and_image(handles)
% retrieve loaded image from load dicom push button 
Image = get(handles.Load_Pushbutton,'UserData');
% retrieve all relevant variables from GUI fields:
[orientation, central_position, res_1, res_2, res_3, method, alpha, ...
    beta, gamma, slab_thickness_mm, view] = get_variables(handles);

% check that unrotated slab is fully contained within original image volume
% before plotting, get image length:
image_length = get_image_length(handles);
% retrieve unrotated orthogonal voxel dimension
orth_vox_dim = get_orthogonal_vox_dim(handles);
% set tolerance for 3DP accuracy
tol = 0.001;
% set upper and lower bounds
lower = central_position - slab_thickness_mm/2 + orth_vox_dim/2 + tol;
upper = central_position + slab_thickness_mm/2 - orth_vox_dim/2 - tol;
    try
if lower <  0 || upper > image_length
    % out of bound measurement has been requested
    error = errordlg({'Slab must be completely contained by image volume.' ...
        'Adjust central slice position and/or slab thickness parameters.',...
        ['Orthogonal image volume dimension is ', sprintf('%.3f',image_length),' mm.']} ...
        ,'Slab Error','modal');
        disp('Oblique slab NOT computed.')  
        % block program execution until user has clicked on modal error box
        uiwait(error)
else
% compute slab, centre index, no slices in slab, rotated voxel dimensions
% and 3D stencil corresponding to slab
disp('Computing oblique slab ...')
[slab, centre_index, no_slices, rotated_vox_dim, stencil_3D] = ComputeObliqueSlab(Image, ...
    orientation, central_position, res_1, res_2, res_3, method, alpha, beta, ...
    gamma, slab_thickness_mm);
disp('Slab computed OK.')

disp('Computing projection image ...')
% compute central slice projection image
projection_image = ComputeProjectionImage(slab, stencil_3D, orientation, ...
    centre_index, view, no_slices);
disp('Projection image computed OK.')

% Plot slab and projection_image in GUI

% 1) Display oblique slice in 2D in Left plot
axes(handles.Left)
% Call DisplayImage2D to plot image into left axes
disp('Updating left hand display ...')
DisplayImage2D(projection_image, rotated_vox_dim, orientation, handles.Left)
disp('Left hand display updated OK.')

% 2) Display slab in 3D in Right plot
axes(handles.Right)
disp('Updating right hand display ...')
DisplaySlab3D(Image, slab, orientation, central_position, alpha, ...
    beta, gamma, centre_index, no_slices,'', handles.Right)
disp('Right hand display updated OK.')

% store variables so can be accessed by plot_image function (for when just
% need to change left display, but not right display)
% Create a Data structure with 6 fields
Data.slab = slab;
Data.centre_index = centre_index;
Data.no_slices = no_slices;
Data.rotated_vox_dim = rotated_vox_dim;
Data.orientation = orientation;
Data.stencil_3D = stencil_3D;
% Store this structure in current 'memory'
set(handles.projection_image_drop_down,'UserData',Data);
end
    catch
    end

% plot image function invoked by changing projection image view only.
function plot_image(handles)
% retrieve slab, centre index, no slices in slab, rotated voxel dimensions
% and 3D stencil corresponding to slab
Data = get(handles.projection_image_drop_down,'UserData');
% retrieve updated projection image view information
switch get(handles.projection_image_drop_down,'Value')
    case 1
        view = 'central';
    case 2
        view = 'max';
    case 3
        view = 'min';
    case 4
        view = 'median';
end
% compute central slice projection image
disp('Computing projection image ...')
projection_image = ComputeProjectionImage(Data.slab, Data.stencil_3D, ...
    Data.orientation, Data.centre_index, view, Data.no_slices);
disp('Projection image computed OK.')

% Plot projection_image in GUI
axes(handles.Left)
disp('Updating left hand display ...')
% Call DisplayImage2D to plot image into left axes
DisplayImage2D(projection_image, Data.rotated_vox_dim, Data.orientation, handles.Left)
disp('Left hand display updated OK.')

% --- Executes on selection change in Interpolation_Method_Popupmenu.
function Interpolation_Method_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Interpolation_Method_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Interpolation_Method_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Interpolation_Method_Popupmenu

% Check if a volume has been loaded or not
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
end

% --- Executes during object creation, after setting all properties.
function Interpolation_Method_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Interpolation_Method_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function beta_edit_box_Callback(hObject, eventdata, handles)
% hObject    handle to beta_edit_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta_edit_box as text
%        str2double(get(hObject,'String')) returns contents of beta_edit_box as a double

% check to see if volume is loaded or not
empty = isempty(get(handles.Load_Pushbutton,'UserData'));

% To show all values to 1dp
beta_string = get(hObject,'String');
% Convert string to double
beta_number = str2double(beta_string);

% Only accept rotation angles within range +-180 degrees
if beta_number >= -180 && beta_number <= 180
    % Convert to 1dp string
    beta_string_1dp = sprintf('%.1f',beta_number);
    % set value to this updated string
    set(hObject,'String', beta_string_1dp);
    % To link beta edit box value to beta slider value
    set(handles.beta_slider,'Value', beta_number)
    if empty == true
        warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
    end
elseif beta_number < -180
    % rotation value is too low
    set(hObject,'String', '-180.0');
    % To link beta edit box value to beta slider value
    set(handles.beta_slider,'Value', -180)
    if empty == true
        warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
    else
        error = errordlg('Minimum rotation value is -180.0','Rotation Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
    end
elseif beta_number > 180
    % rotation value is too high
    set(hObject,'String', '180.0');
    % To link beta edit box value to beta slider value
    set(handles.beta_slider,'Value', 180)
    if empty == true
        warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
    else
        error = errordlg('Maximum rotation value is 180.0','Rotation Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
    end
end

% --- Executes during object creation, after setting all properties.
function beta_edit_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_edit_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gamma_edit_box_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_edit_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gamma_edit_box as text
%        str2double(get(hObject,'String')) returns contents of gamma_edit_box as a double

% check to see if volume is loaded or not
empty = isempty(get(handles.Load_Pushbutton,'UserData'));

% To show all values to 1dp
gamma_string = get(hObject,'String');
% Convert string to double
gamma_number = str2double(gamma_string);

% Only accept rotation angles within range +-180 degrees
if gamma_number >= -180 && gamma_number <= 180
    % Convert to 1dp string
    gamma_string_1dp = sprintf('%.1f',gamma_number);
    % set value to this updated string
    set(hObject,'String', gamma_string_1dp);
    % To link gamma edit box value to gamma slider value
    set(handles.gamma_slider,'Value', gamma_number)
    if empty == true
        warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
    end
elseif gamma_number < -180
    % rotation value is too low
    set(hObject,'String', '-180.0');
    % To link beta gamma box value to gamma slider value
    set(handles.gamma_slider,'Value', -180)
    if empty == true
        warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
    else
        error = errordlg('Minimum rotation value is -180.0','Rotation Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
    end
elseif gamma_number > 180
    % rotation value is too high
    set(hObject,'String', '180.0');
    % To link gamma edit box value to gamma slider value
    set(handles.gamma_slider,'Value', 180)
    if empty == true
        warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
    else
        error = errordlg('Maximum rotation value is 180.0','Rotation Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
    end
end

% --- Executes during object creation, after setting all properties.
function gamma_edit_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_edit_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function alpha_slider_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% To link alpha slider value to alpha edit box
alpha = get(hObject,'Value');
% Convert number to string with 1dp precision
alpha = num2str (alpha,'%.1f');
set(handles.alpha_edit_box,'String', alpha)

% Check if a volume is loaded or not
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
end

% --- Executes during object creation, after setting all properties.
function alpha_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function beta_slider_Callback(hObject, eventdata, handles)
% hObject    handle to beta_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% To link beta slider value to beta edit box
beta = get(hObject,'Value');
% Convert number to string with 1dp precision
beta = num2str (beta,'%.1f');
set(handles.beta_edit_box,'String', beta)

% Check if a volume is loaded or not
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
end

% --- Executes during object creation, after setting all properties.
function beta_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function gamma_slider_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% To link gamma slider value to gamma edit box
gamma = get(hObject,'Value');
% Convert number to string with 1dp precision
gamma = num2str (gamma,'%.1f');
set(handles.gamma_edit_box,'String', gamma)

% Check if a volume is loaded or not
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
end

% --- Executes during object creation, after setting all properties.
function gamma_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function alpha_edit_box_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_edit_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_edit_box as text
%        str2double(get(hObject,'String')) returns contents of alpha_edit_box as a double

% check to see if volume is loaded or not
empty = isempty(get(handles.Load_Pushbutton,'UserData'));

% To show all values to 1dp
alpha_string = get(hObject,'String');
% Convert string to double
alpha_number = str2double(alpha_string);

% Only accept rotation angles within range +-180 degrees
if alpha_number >= -180 && alpha_number <= 180
    % Convert to 1dp string
    alpha_string_1dp = sprintf('%.1f',alpha_number);
    % set value to this updated string
    set(hObject,'String', alpha_string_1dp);
    % To link alpha edit box value to alpha slider value
    set(handles.alpha_slider,'Value', alpha_number)
    if empty == true
        warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
    end
elseif alpha_number < -180
    % rotation value is too low
    set(hObject,'String', '-180.0');
    % To link alpha edit box value to alpha slider value
    set(handles.alpha_slider,'Value', -180)
    if empty == true
        warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
    else
        error = errordlg('Minimum rotation value is -180.0','Rotation Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
    end
elseif alpha_number > 180
    % rotation value is too high
    set(hObject,'String', '180.0');
    % To link alpha edit box value to alpha slider value
    set(handles.alpha_slider,'Value', 180)
    if empty == true
        warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
    else
        error = errordlg('Maximum rotation value is 180.0','Rotation Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
    end
end

% --- Executes during object creation, after setting all properties.
function alpha_edit_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_edit_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in projection_image_drop_down.
function projection_image_drop_down_Callback(hObject, eventdata, handles)
% hObject    handle to projection_image_drop_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns projection_image_drop_down contents as cell array
%        contents{get(hObject,'Value')} returns selected item from projection_image_drop_down

% Check if a volume is loaded or not
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');    
else
    % Update left hadn display
    plot_image(handles)
end

% --- Executes during object creation, after setting all properties.
function projection_image_drop_down_CreateFcn(hObject, eventdata, handles)
% hObject    handle to projection_image_drop_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function slab_thickness_slider_Callback(hObject, eventdata, handles)
% hObject    handle to slab_thickness_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% To link slab thickness slider value to slab thickness edit box
thickness = get(hObject,'Value');
% Convert number to string with 3dp precision
thickness = num2str (thickness,'%.3f');
set(handles.slab_thickness_edit_box,'String', thickness)

% Check to see if volume is loaded or not
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
end

% --- Executes during object creation, after setting all properties.
function slab_thickness_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slab_thickness_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slab_thickness_edit_box_Callback(hObject, eventdata, handles)
% hObject    handle to slab_thickness_edit_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slab_thickness_edit_box as text
%        str2double(get(hObject,'String')) returns contents of slab_thickness_edit_box as a double

% To show all values to 3dp
thickness_string = get(hObject,'String');
% Convert string to double
thickness_number = str2double(thickness_string);

% check to see if volume is loaded or not
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    % retrieve slice position from slice position edit box
    if thickness_number < 1
        set(hObject,'String', '1.000');
        % To link slab thickness edit box value to slab thickness slider value
        set(handles.slab_thickness_slider,'Value', 1)
    elseif thickness_number > 100
        set(hObject,'String', '100.000');
        % To link slab thickness edit box value to slab thickness slider value
        set(handles.slab_thickness_slider,'Value', 100)
    else
        thickness_string_3dp = sprintf('%.3f', thickness_number);
        set(hObject,'String', thickness_string_3dp);
        % To link slab thickness edit box value to slab thickness slider value
        set(handles.slab_thickness_slider,'Value', thickness_number)
    end
        warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
else
    %get upper and lower bounds for slice thickness
    min = get(handles.slab_thickness_slider,'Min');
    max = get(handles.slab_thickness_slider,'Max');
    if thickness_number >= min && thickness_number <= max
        % thickness number is ok - convert to 3dp string
        thickness_string_3dp = sprintf('%.3f',thickness_number);
        % set value to this updated string
        set(hObject,'String', thickness_string_3dp);
        % To link slab thickness edit box value to slab thickness slider value
        set(handles.slab_thickness_slider,'Value', thickness_number)
    elseif thickness_number < min
        % thickness number is too low
        thickness_string_3dp = sprintf('%.3f', min);
        set(hObject,'String', thickness_string_3dp);
        % To link slab thickness edit box value to slab thickness slider value
        set(handles.slab_thickness_slider,'Value', min)
        error = errordlg...
            ({'Minimum slab thickness is the unrotated orthogonal voxel';...
            ['dimension (',sprintf('%.3f', min),' mm).'];...
            'Slab thickness reverted to minimum value.'}, ...
            'Slab Thickness Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
    elseif thickness_number > max
        % set value of edit box to max 
        thickness_string_3dp = sprintf('%.3f', max);
        set(hObject,'String', thickness_string_3dp);
        % To link slab thickness edit box value to slab thickness slider value
        set(handles.slab_thickness_slider,'Value', max)
        % slice thickness value is too high
        error = errordlg({['Maximum slab thickness is the orthogonal volume length (',...
            sprintf('%.3f', max),' mm).'];...
            'Slab thickness reverted to maximum value.'},...
            'Slab Thickness Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
    end
end

% --- Executes during object creation, after setting all properties.
function slab_thickness_edit_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slab_thickness_edit_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Res_3_Edit_Box_Callback(hObject, eventdata, handles)
% hObject    handle to Res_3_Edit_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Res_3_Edit_Box as text
%        str2double(get(hObject,'String')) returns contents of Res_3_Edit_Box as a double

% To show all values to 3dp
res_3_string = get(hObject,'String');
% Convert string to double
res_3_number = str2double(res_3_string);

% check to see if a volume is loaded
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    if res_3_number >= 0.5 && res_3_number <= 20
        % Convert to 3dp string
        res_3_string_3dp = sprintf('%.3f',res_3_number);
        % set value to this updated string
        set(hObject,'String', res_3_string_3dp);
    elseif res_3_number < 0.5
        % resolution value is too low
        set(hObject,'String', '0.500');
    elseif res_3_number > 20
        % resolution value is too high
        set(hObject,'String', '20.000');
    end
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
else
    % Only accept sensible resolution values
    if res_3_number >= 0.5 && res_3_number <= 20
        % Convert to 3dp string
        res_3_string_3dp = sprintf('%.3f',res_3_number);
        % set value to this updated string
        set(hObject,'String', res_3_string_3dp);    
    elseif res_3_number < 0.5
        error = errordlg('Minimum resolution value is 0.500','Resolution Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
        set(hObject,'String', '0.500');
    elseif res_3_number > 20
        % resolution value is too high
        error = errordlg('Maximum resolution value is 20.000','Resolution Error','modal');
        % block program execution until user has clicked on modal error box
        uiwait(error)
        set(hObject,'String', '20.000');
    end
end

% --- Executes during object creation, after setting all properties.
function Res_3_Edit_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Res_3_Edit_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_default_values(handles)
% Set default values for all variables when starting function or when
% clicking reset with no volume loaded
set(handles.Slice_Position_Edit_Box,'String', '50.000');
set(handles.Slice_Position_Slider,'Min',0,'Max',100,'Value', 50);
% set minor slider step to 1 unit, and major slider step to 10 units
set(handles.Slice_Position_Slider,'SliderStep',[1, 10] / 100);
set(handles.Res_1_Edit_Box,'String', '1.000');
set(handles.Res_2_Edit_Box,'String', '1.000');
set(handles.Res_3_Edit_Box,'String', '1.000');
set(handles.Plane_Pop_Up_Menu,'Value', 1);
set(handles.Interpolation_Method_Popupmenu,'Value', 1);
% alpha edit box and slider
set(handles.alpha_edit_box,'String', '0.0');
set(handles.alpha_slider,'Min',-180,'Max',180,'Value', 0);
% set minor slider step to 1 degree, and major slider step to 10 degrees
set(handles.alpha_slider,'SliderStep',[1, 10] / (360));
% beta edit box and slider
set(handles.beta_edit_box,'String', '0.0');
set(handles.beta_slider,'Min',-180,'Max',180,'Value', 0);
set(handles.beta_slider,'SliderStep',[1, 10] / (360));
% gamma edit box and slider
set(handles.gamma_edit_box,'String', '0.0');
set(handles.gamma_slider,'Min',-180,'Max',180,'Value', 0);
set(handles.gamma_slider,'SliderStep',[1, 10] / (360));
% slab thickness and slider
set(handles.slab_thickness_edit_box,'String', '1.000');
set(handles.slab_thickness_slider,'Min',1,'Max',100,'Value', 1);
% set minor slider step to 1 unit, and major slider step to 10 units
set(handles.slab_thickness_slider,'SliderStep',[1, 10] / 99);
% projection image popupmenu
set(handles.projection_image_drop_down,'Value', 1);
set(handles.Res_1_Text, 'String', 'X (mm):');
set(handles.Res_2_Text, 'String', 'Y (mm):');
set(handles.Res_3_Text, 'String', 'Z (mm):');
set(handles.Slice_Position_Text, 'String', 'Central Slice Position (Z / mm):');

function set_loaded_image_values(handles)
% Set default values when first loading an image and when clicking reset
% when an image is loaded
Loaded_Image = get(handles.Load_Pushbutton,'UserData');
% 1) To extract voxel dimensions in mm (used to scale axes in plots, for
% default resolutions, and for slider increment values)
vox_dim = Loaded_Image.VoxelDimensions; % [dy dx dz]
% 2) To get size of image as 3D vector in form of (y,x,z) i.e.
% (number rows, number of columns, number of slices)
image_dim = size(Loaded_Image.ImageData);
% 3) To get size of image as 3D vector in mm in form of (y,x,z) i.e.
% (image height, image width, image length) (used as maximum range for
% slider)
image_size = (image_dim - [1 1 1]).*vox_dim;
position = num2str(image_size(3)/2,'%.3f');
% set position of slice position edit box and slider
set(handles.Slice_Position_Edit_Box,'String', position);
set(handles.Slice_Position_Slider,'Min',0,'Max',image_size(3),'Value', image_size(3)/2);
% set minor slider step to 1 slice, and major slider step to 10 slices
set(handles.Slice_Position_Slider,'SliderStep',[1, 10] / (image_dim(3)-1));
res_1 = num2str (vox_dim(2),'%.3f');
res_2 = num2str (vox_dim(1),'%.3f');
res_3 = num2str (vox_dim(3),'%.3f');
set(handles.Res_1_Edit_Box,'String', res_1);
set(handles.Res_2_Edit_Box,'String', res_2);
set(handles.Res_3_Edit_Box,'String', res_3);
% alpha edit box and slider
set(handles.alpha_edit_box,'String', '0.0');
set(handles.alpha_slider,'Value', 0);
% beta edit box and slider
set(handles.beta_edit_box,'String', '0.0');
set(handles.beta_slider,'Value', 0)
% gamma edit box and slider
set(handles.gamma_edit_box,'String', '0.0');
set(handles.gamma_slider,'Value', 0);
% projection drop_down menu
set(handles.projection_image_drop_down,'Value', 1);
% interpolation drop_down menu
set(handles.Interpolation_Method_Popupmenu,'Value', 1);
% set view plane menu
set(handles.Plane_Pop_Up_Menu,'Value', 1);
set(handles.Res_1_Text, 'String', 'X (mm):');
set(handles.Res_2_Text, 'String', 'Y (mm):');
set(handles.Res_3_Text, 'String', 'Z (mm):');
set(handles.Slice_Position_Text, 'String', 'Central Slice Position (Z / mm):');
% slab_thickness
set(handles.slab_thickness_slider,'Min',vox_dim(3),'Max',image_size(3));
% set minor slider step to add 1 slice to slab, and major slider step to
% add 10 slices to slab
set(handles.slab_thickness_slider,'SliderStep',[1, 10] / (image_dim(3)-2));
set(handles.slab_thickness_edit_box,'String', res_3);
set(handles.slab_thickness_slider,'Value', vox_dim(3))

function [orientation, central_position, res_1, res_2, res_3, method, alpha, ...
    beta, gamma, slab_thickness_mm, view] = get_variables(handles)
%Function to retrieve all variables when either plotting both left and
%right views or just updating left view

% retrieve slice orientation from view plane drop down menu
switch get(handles.Plane_Pop_Up_Menu,'Value')
    case 1
        orientation = 'X-Y';
    case 2
        orientation = 'Y-Z';
    case 3
        orientation = 'X-Z';
end
% retrieve central slice position from central slice position edit box
central_position = str2double(get(handles.Slice_Position_Edit_Box,'String'));
% retrieve resolution 1 from res_1 edit box
res_1 = str2double(get(handles.Res_1_Edit_Box,'String'));
% retrieve resolution 2 from res_2 edit box
res_2 = str2double(get(handles.Res_2_Edit_Box,'String'));
% retrieve resolution 3 from res_3 edit box
res_3 = str2double(get(handles.Res_3_Edit_Box,'String'));
% retrieve interpolation method from drop down menu
switch get(handles.Interpolation_Method_Popupmenu,'Value')
    case 1
        method = 'nearest';
    case 2
        method = 'linear';
    case 3
        method = 'spline';
end
% retrieve three rotation parameters
alpha = str2double(get(handles.alpha_edit_box,'String'));
beta = str2double(get(handles.beta_edit_box,'String'));
gamma = str2double(get(handles.gamma_edit_box,'String'));
% retrieve slab thickness
slab_thickness_mm = str2double(get(handles.slab_thickness_edit_box,'String'));
% retrieve projection image view 
switch get(handles.projection_image_drop_down,'Value')
    case 1
        view = 'central';
    case 2
        view = 'max';
    case 3
        view = 'min';
    case 4
        view = 'median';
end

function length = get_image_length(handles)
% retrieve loaded image from load dicom push button 
Image = get(handles.Load_Pushbutton,'UserData');
% 1) To extract voxel dimensions in mm
vox_dim = Image.VoxelDimensions; % [dy dx dz]
% 2) To get size of image as 3D vector in form of (y,x,z)
image_dim = size(Image.ImageData);
% 3) To get size of image as 3D vector in mm in form of (y,x,z)
image_size = (image_dim - [1 1 1]).*vox_dim;
% retrieve slice orientation from view plane drop down menu
switch get(handles.Plane_Pop_Up_Menu,'Value')
    case 1 
        % XY plane
        length = image_size(3);
    case 2
        % YZ plane
        length = image_size(2);
    case 3
        % XZ plane
        length = image_size(1);
end

function orth_vox_dim = get_orthogonal_vox_dim(handles)
% retrieve loaded image from load dicom push button 
Image = get(handles.Load_Pushbutton,'UserData');
% 1) To extract voxel dimensions in mm
vox_dim = Image.VoxelDimensions; % [dy dx dz]
% retrieve slice orientation from view plane drop down menu
switch get(handles.Plane_Pop_Up_Menu,'Value')
    case 1 
        % XY plane
        orth_vox_dim = vox_dim(3);
    case 2
        % YZ plane
        orth_vox_dim = vox_dim(2);
    case 3
        % XZ plane
        orth_vox_dim = vox_dim(1);
end

% --- Executes on button press in refresh_button.
function refresh_button_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check to see if a volume is loaded
empty = isempty(get(handles.Load_Pushbutton,'UserData'));
if empty == true
    warndlg('Please load a DICOM volume into the GUI.','No DICOM Loaded','replace');
else
    try
    % update display of GUI
    plot_slab_and_image(handles)
    catch
    % error message shown in case of orthogonal resolution error    
    end
end
