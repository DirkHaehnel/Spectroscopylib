function varargout = query_metadata(varargin)
% QUERY_METADATA M-file for query_metadata.fig
%      QUERY_METADATA by itself, creates a new QUERY_METADATA or raises the
%      existing singleton*.
%
%      H = QUERY_METADATA returns the handle to a new QUERY_METADATA or the handle to
%      the existing singleton*.
%
%      QUERY_METADATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUERY_METADATA.M with the given input arguments.
%
%      QUERY_METADATA('Property','Value',...) creates a new QUERY_METADATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before query_metadata_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to query_metadata_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help query_metadata

% Last Modified by GUIDE v2.5 19-Dec-2007 08:55:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @query_metadata_OpeningFcn, ...
                   'gui_OutputFcn',  @query_metadata_OutputFcn, ...
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

% --- Executes just before query_metadata is made visible.
function query_metadata_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to query_metadata (see VARARGIN)

% Choose default command line output for query_metadata
handles.output = struct();

% Map betweeen stimulating wavelength and coresponding distance
handles.pinholes = {
    '200' 200;
    '150' 150;
    '100' 100;
    '75' 75;
    '50' 50
};

handles.wavelength_distances = {
    '470' 470 432 '510';
%    '532' 532 389;
    '637' 637 470 '675'
    };

wavelengths = strvcat(handles.wavelength_distances{:,1});

set(handles.pinhole, 'String', strvcat(handles.pinholes{:,1}));
set(handles.lambdaex_menu1, 'String', wavelengths);
set(handles.lambdaex_menu2, 'String', wavelengths);
set(handles.distance, 'String', handles.wavelength_distances{1, 3});
set(handles.distance2, 'String', handles.wavelength_distances{1, 3});

handles.output.LambdaEm1 = str2double(get(handles.lambda_em1, 'String'));
handles.output.LambdaEm2 = str2double(get(handles.lambda_em2, 'String'));

if size(varargin, 2) > 1
    if varargin{1,2}.Temperature == 0
         varargin{1,2}.Temperature = 293.15;
    end
    defaults = varargin{2};
    handles.output = defaults;
%    if isfield(defaults, 'Temperature')
        set(handles.temperature, 'String', defaults.Temperature);
%    end
    
    if isfield(defaults, 'Pinhole')
        pinhole_index = find(handles.pinholes == defaults.Pinhole);
        if size(pinhole_index) > 0
            set(handles.pinhole, 'Value', pinhole_index(2));
        end
    else
        set(handles.pinhole, 'Value', 2);
    end
    if isfield(defaults, 'LambdaEx1')
        index = find(handles.wavelength_distances{:,2} == defaults.LambdaEx1);
        if size(index) > 0
            set(handles.lambdaex_menu1, 'Value', index);
            set(handles.distance, 'String', handles.wavelength_distances{index,3});
        end
    else
        set(handles.lambdaex_menu1, 'Value', 1);
        set(handles.distance, 'String', handles.wavelength_distances{1,3});
    end
    if isfield(defaults, 'LambdaEx2')
        index = find(handles.wavelength_distances{:,2} == defaults.LambdaEx2);
        if size(index) > 0
            set(handles.lambdaex_menu2, 'Value', index);
            set(handles.distance2, 'String', handles.wavelength_distances{index,3});
        end
    else
        set(handles.lambdaex_menu2, 'Value', 2);
        set(handles.distance2, 'String', handles.wavelength_distances{2,3});
    end
    if isfield(defaults, 'LambdaEm1')
        set(handles.lambda_em1, 'String', defaults.LambdaEm1);
    else
        set(handles.lambda_em1, 'String', handles.wavelength_distances{1,4});
    end
    if isfield(defaults, 'LambdaEm2')
        set(handles.lambda_em2, 'String', handles.wavelength_distances{2,4});
    else
        set(handles.lambda_em2, 'String', '675');
    end
    if size(varargin, 2) > 2
        if varargin{3} == 1
           set(handles.lambdaex_menu2, 'Visible', 'off');
           set(handles.lambda_em2, 'Visible', 'off');
           set(handles.distance2, 'Visible', 'off');
           set(handles.distance2_label, 'Visible', 'off');
           set(handles.lambdaem2_label, 'Visible', 'off');
           set(handles.lambdaex2_label, 'Visible', 'off');
           set(handles.wavelength2_label, 'Visible', 'off');
        else
           set(handles.lambdaex_menu2, 'Visible', 'on');
           set(handles.lambda_em2, 'Visible', 'on');
           set(handles.distance2, 'Visible', 'on');
           set(handles.distance2_label, 'Visible', 'on');
           set(handles.lambdaem2_label, 'Visible', 'on');
           set(handles.lambdaex2_label, 'Visible', 'on');
           set(handles.wavelength2_label, 'Visible', 'on');
        end
    end
end

% Update handles structure
guidata(hObject, handles);

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);

    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);

% Make the GUI modal
set(handles.query_metadata,'WindowStyle','modal')

% UIWAIT makes query_metadata wait for user response (see UIRESUME)
uiwait(handles.query_metadata);

% --- Outputs from this function are returned to the command line.
function varargout = query_metadata_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.query_metadata);

% --- Executes when user attempts to close query_metadata.
function query_metadata_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to query_metadata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(handles.query_metadata, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.query_metadata);
else
    % The GUI is no longer waiting, just close it
    delete(handles.query_metadata);
end


% --- Executes on key press over query_metadata with no controls selected.
function query_metadata_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to query_metadata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for "enter" or "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    % User said no by hitting escape
    handles.output = struct;
    
    % Update handles structure
    guidata(hObject, handles);
    
    uiresume(handles.query_metadata);
end    
    
if isequal(get(hObject,'CurrentKey'),'return')
    uiresume(handles.query_metadata);
end    

function lambda_ex1_Callback(hObject, eventdata, handles)
function lambda_ex1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lambda_ex2_Callback(hObject, eventdata, handles)
function lambda_ex2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lambda_em1_Callback(hObject, eventdata, handles)
function lambda_em1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lambda_em2_Callback(hObject, eventdata, handles)
function lambda_em2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pinhole_Callback(hObject, eventdata, handles)
function pinhole_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function distance_Callback(hObject, eventdata, handles)
function distance_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function accept_metadata_Callback(hObject, eventdata, handles)
handles.output = struct();
index1 = get(handles.lambdaex_menu1, 'Value');
index2 = get(handles.lambdaex_menu2, 'Value');
handles.output.LambdaEx1 = handles.wavelength_distances{index1, 2};
handles.output.LambdaEx2 = handles.wavelength_distances{index2, 2};
handles.output.LambdaEm1 = str2double(get(handles.lambda_em1, 'String'));
handles.output.LambdaEm2 = str2double(get(handles.lambda_em2, 'String'));
handles.output.Distance1 = handles.wavelength_distances{index1, 3};
handles.output.Distance2 = handles.wavelength_distances{index2, 3};
pinhole_index = get(handles.pinhole, 'Value');
handles.output.Pinhole = handles.pinholes{pinhole_index, 2};
handles.output.Temperature = str2double(get(handles.temperature, 'String'));

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.query_metadata);

function distance2_Callback(hObject, eventdata, handles)
function distance2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lambdaex_menu1_Callback(hObject, eventdata, handles)
index = get(hObject, 'Value');
set(handles.distance, 'String', handles.wavelength_distances{index, 3});
set(handles.lambda_em1, 'String', handles.wavelength_distances{index,4});

function lambdaex_menu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lambdaex_menu2_Callback(hObject, eventdata, handles)
index = get(hObject, 'Value');
set(handles.distance2, 'String', handles.wavelength_distances{index, 3});
set(handles.lambda_em2, 'String', handles.wavelength_distances{index,4});

function lambdaex_menu2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function temperature_Callback(hObject, eventdata, handles)
function temperature_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
