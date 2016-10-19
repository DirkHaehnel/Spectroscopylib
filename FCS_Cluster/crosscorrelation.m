function varargout = crosscorrelation(varargin)
% CROSSCORRELATION M-file for crosscorrelation.fig
%      CROSSCORRELATION, by itself, creates a new CROSSCORRELATION or raises the existing
%      singleton*.
%
%      H = CROSSCORRELATION returns the handle to a new CROSSCORRELATION or the handle to
%      the existing singleton*.
%
%      CROSSCORRELATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CROSSCORRELATION.M with the given input arguments.
%
%      CROSSCORRELATION('Property','Value',...) creates a new CROSSCORRELATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crosscorrelation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crosscorrelation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help crosscorrelation

% Last Modified by GUIDE v2.5 09-Nov-2010 17:08:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @crosscorrelation_OpeningFcn, ...
                   'gui_OutputFcn',  @crosscorrelation_OutputFcn, ...
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

% --- Executes just before crosscorrelation is made visible.
function crosscorrelation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crosscorrelation (see VARARGIN)

% Choose default command line output for crosscorrelation
handles.output = hObject;
handles.pt3fn = '';
handles.last_path = '';
handles.mcs_zoom_x = 'auto';
handles.mcs_zoom_y = 'auto';
handles.lifetime_zoom_x = 'auto';
handles.lifetime_zoom_y = 'auto';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crosscorrelation wait for user response (see UIRESUME)
% uiwait(handles.crosscorrelation);


% --- Outputs from this function are returned to the command line.
function varargout = crosscorrelation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function mcs_binning_edit_Callback(hObject, eventdata, handles)
handles = update_display_mcs(handles);
handles = update_graph(handles);
guidata(hObject, handles);

function mcs_binning_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function zoom_area_Callback(hObject, eventdata, handles)
old_context = gca;
axes(handles.mcs_axes);
[x y] = ginput(2);
x=sort(x);
y=sort(y);
xlim(x);
ylim(y);
handles.mcs_zoom_x = x;
handles.mcs_zoom_y = y;
axes(old_context);
guidata(hObject, handles);

function reset_zoom_Callback(hObject, eventdata, handles)
old_context = gca;
axes(handles.mcs_axes);
xlim('auto');
ylim('auto');
handles.mcs_zoom_x = 'auto';
handles.mcs_zoom_y = 'auto';
axes(old_context);
guidata(hObject, handles);

function mcs_timegate1_edit_Callback(hObject, eventdata, handles)
handles = update_graph(handles);
guidata(hObject, handles);

function mcs_timegate1_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mcs_timegate2_edit_Callback(hObject, eventdata, handles)
handles = update_graph(handles);
guidata(hObject, handles);

function mcs_timegate2_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function reset_tg_button_Callback(hObject, eventdata, handles)
set(handles.mcs_timegate1_edit, 'String', 0);
set(handles.mcs_timegate2_edit, 'String', Inf);
handles = update_graph(handles);
guidata(hObject, handles);

function select_tg_area_button_Callback(hObject, eventdata, handles)
[x y] = ginput(2);
x=round(sort(x)*100)/100;
set(handles.mcs_timegate1_edit, 'String', x(1));
set(handles.mcs_timegate2_edit, 'String', x(2));
handles = update_graph(handles);
guidata(hObject, handles);

function channel1_checkbox_Callback(hObject, eventdata, handles)
handles = update_display_mcs(handles);
handles = update_display_lifetimedata(handles);
handles = update_graph(handles);
guidata(hObject, handles);

function channel2_checkbox_Callback(hObject, eventdata, handles)
handles = update_display_mcs(handles);
handles = update_display_lifetimedata(handles);
handles = update_graph(handles);
guidata(hObject, handles);

function channel3_checkbox_Callback(hObject, eventdata, handles)
handles = update_display_mcs(handles);
handles = update_display_lifetimedata(handles);
handles = update_graph(handles);
guidata(hObject, handles);

function channel4_checkbox_Callback(hObject, eventdata, handles)
handles = update_display_mcs(handles);
handles = update_display_lifetimedata(handles);
handles = update_graph(handles);
guidata(hObject, handles);

function select_all_channels_Callback(hObject, eventdata, handles)
set(handles.channel1_checkbox, 'Value', 1.0);
set(handles.channel2_checkbox, 'Value', 1.0);
set(handles.channel3_checkbox, 'Value', 1.0);
set(handles.channel4_checkbox, 'Value', 1.0);
handles = update_display_mcs(handles);
handles = update_display_lifetimedata(handles);
handles = update_graph(handles);
guidata(hObject, handles);

function unselect_all_channels_Callback(hObject, eventdata, handles)
set(handles.channel1_checkbox, 'Value', 0);
set(handles.channel2_checkbox, 'Value', 0);
set(handles.channel3_checkbox, 'Value', 0);
set(handles.channel4_checkbox, 'Value', 0);
handles = update_display_mcs(handles);
handles = update_display_lifetimedata(handles);
handles = update_graph(handles);
guidata(hObject, handles);

function import_pt3_Callback(hObject, eventdata, handles)
[filename, filepath] = uigetfile('*.ht3', 'Read HT3 File', handles.last_path);
if ~isfloat(filename)
    handles.pt3fn = fullfile(filepath, filename);
    handles.last_path = filepath();
    tg.tg1 = [0 Inf];
    [handles.mcs handles.lifetimedata handles.head] = part_mcs(handles.pt3fn, 3000000, tg, handles.status_text);
    set(handles.mcs_timegate1_edit, 'String', 0);
    set(handles.mcs_timegate2_edit, 'String', Inf);
    set(handles.status_text, 'String', '');
    handles = update_display_mcs(handles);
    handles = recalculate_lifetimedata(handles, 1);
    handles = update_display_lifetimedata(handles);
    handles = update_graph(handles);
end
guidata(hObject, handles);

function handles = update_graph(handles)
handles = update_mcs_graph(handles);
handles = update_lifetimedata_graph(handles);

function handles = update_mcs_graph(handles)
if ~ isfield(handles, 'display_mcs')
    return
end
old_context = gca;
axes(handles.mcs_axes);
tg1 = str2double(get(handles.mcs_timegate1_edit, 'String'));
tg2 = str2double(get(handles.mcs_timegate2_edit, 'String'));
entry_count = 15;
if (tg1 ~= Inf) && (tg1 > 0)
    entry_count = entry_count + 3;
end
if (tg2 ~= Inf) && (tg2 > 0)
    entry_count = entry_count + 3;
end

plotparam = cell(1, entry_count);
plotparam{1} = handles.display_mcs(:,1);
plotparam{2} = handles.display_mcs(:,6);
plotparam{3} = 'b';
plotparam{4} = handles.display_mcs(:,1);
plotparam{5} = handles.display_mcs(:,2);
plotparam{6} = 'r';
plotparam{7} = handles.display_mcs(:,1);
plotparam{8} = handles.display_mcs(:,3);
plotparam{9} = 'g';
plotparam{10} = handles.display_mcs(:,1);
plotparam{11} = handles.display_mcs(:,4);
plotparam{12} = 'c';
plotparam{13} = handles.display_mcs(:,1);
plotparam{14} = handles.display_mcs(:,5);
plotparam{15} = 'm';
cur_entry = 16;
if (tg1 ~= Inf) && (tg1 > 0)
    plotparam{cur_entry} = [tg1 tg1];
    plotparam{cur_entry + 1} = [1 max(handles.display_mcs(:,6))*1.1];
    plotparam{cur_entry + 2} = '--rs';
    cur_entry = cur_entry + 3;
end
if (tg2 ~= Inf) && (tg2 > 0)
    plotparam{cur_entry} = [tg2 tg2];
    plotparam{cur_entry + 1} = [1 max(handles.display_mcs(:,6))*1.1];
    plotparam{cur_entry + 2} = '--rs';
    cur_entry = cur_entry + 3;
end

semilogy(plotparam{:});

xlabel('time [s]');
ylabel('counts');

xlim(handles.mcs_zoom_x);
ylim(handles.mcs_zoom_y);

axes(old_context);

function handles = update_display_mcs(handles)
if ~ isfield(handles, 'mcs')
    return
end
binning = round(str2double(get(handles.mcs_binning_edit, 'String')));
set(handles.mcs_binning_edit, 'String', binning);

if binning ~= 1
    rows = ceil(size(handles.mcs, 1) / binning);
    handles.display_mcs = zeros(rows, 6);
    ok = zeros(1,4);
    ok(1) = get(handles.channel1_checkbox, 'Value') == 1.0;
    ok(2) = get(handles.channel2_checkbox, 'Value') == 1.0;
    ok(3) = get(handles.channel3_checkbox, 'Value') == 1.0;
    ok(4) = get(handles.channel4_checkbox, 'Value') == 1.0;
    for i=1:rows
        for j=2:5
            if ok(j-1)
                handles.display_mcs(i,j)=sum(handles.mcs(((i-1)*binning+1):(i*binning),j));
            end
        end
        handles.display_mcs(i,1)=handles.mcs(i*binning,1);
    end
else
    handles.display_mcs = zeros(size(handles.mcs));
    handles.display_mcs(:,1) = handles.mcs(:,1);
    if get(handles.channel1_checkbox, 'Value') == 1.0
        handles.display_mcs(:,2) = handles.mcs(:,2);
    end
    if get(handles.channel2_checkbox, 'Value') == 1.0
        handles.display_mcs(:,3) = handles.mcs(:,3);
    end
    if get(handles.channel3_checkbox, 'Value') == 1.0
        handles.display_mcs(:,4) = handles.mcs(:,4);
    end
    if get(handles.channel4_checkbox, 'Value') == 1.0
        handles.display_mcs(:,5) = handles.mcs(:,5);
    end   
end
handles.display_mcs(:,6) = handles.display_mcs(:,2) + handles.display_mcs(:,3) + handles.display_mcs(:,4) +handles.display_mcs(:,5);

function handles = recalculate_lifetimedata(handles, dont_read)
if exist('dont_read', 'var') && dont_read > 0
    % well we use the existing data (usefull mostly when freshly importing a
    % file
else
    tg1 = str2double(get(handles.mcs_timegate1_edit, 'String'));
    tg2 = str2double(get(handles.mcs_timegate2_edit, 'String'));
    [handles.lifetimedata, handles.head] = part_lifetimedata(handles.pt3fn, 3000000, [tg1 tg2], handles.status_text);
    set(handles.status_text, 'String', '');
end

function handles = update_display_lifetimedata(handles)
if ~ isfield(handles, 'lifetimedata')
    return
end
handles.display_lifetimedata = zeros(size(handles.lifetimedata));
if get(handles.display_unit_slot, 'Value') == 1.0
    handles.display_lifetimedata(:,1) = handles.lifetimedata(:,1) / handles.head.Resolution;
else
    handles.display_lifetimedata(:,1) = handles.lifetimedata(:,1);
end
if get(handles.channel1_checkbox, 'Value') == 1.0
    handles.display_lifetimedata(:,2) = handles.lifetimedata(:,2);
end
if get(handles.channel2_checkbox, 'Value') == 1.0
    handles.display_lifetimedata(:,3) = handles.lifetimedata(:,3);
end
if get(handles.channel3_checkbox, 'Value') == 1.0
    handles.display_lifetimedata(:,4) = handles.lifetimedata(:,4);
end
if get(handles.channel4_checkbox, 'Value') == 1.0
    handles.display_lifetimedata(:,5) = handles.lifetimedata(:,5);
end   
handles.display_lifetimedata(:,6) =  handles.display_lifetimedata(:,2) + handles.display_lifetimedata(:,3) + handles.display_lifetimedata(:,4) +handles.display_lifetimedata(:,5);

function handles = update_lifetimedata_graph(handles)
if ~ isfield(handles, 'lifetimedata')
    return
end
old_context = gca;
axes(handles.lifetime_axes);
plotparam = cell(1, 15);
plotparam{1} = handles.display_lifetimedata(:,1);
plotparam{2} = handles.display_lifetimedata(:,6);
plotparam{3} = 'b';
plotparam{4} = handles.display_lifetimedata(:,1);
plotparam{5} = handles.display_lifetimedata(:,2);
plotparam{6} = 'r';
plotparam{7} = handles.display_lifetimedata(:,1);
plotparam{8} = handles.display_lifetimedata(:,3);
plotparam{9} = 'g';
plotparam{10} = handles.display_lifetimedata(:,1);
plotparam{11} = handles.display_lifetimedata(:,4);
plotparam{12} = 'c';
plotparam{13} = handles.display_lifetimedata(:,1);
plotparam{14} = handles.display_lifetimedata(:,5);
plotparam{15} = 'm';
semilogy(plotparam{:});
if get(handles.display_unit_slot, 'Value') == 1.0
    xlabel('slots');
else
    xlabel('ns');
end
ylabel('count');
hold on;
for i=1:4
    if get(handles.(sprintf('tg%i_show', i)), 'Value') ~= 1.0
        continue
    end
    tg_start = str2double(get(handles.(sprintf('tg%i_start', i)), 'String'));
    tg_end = str2double(get(handles.(sprintf('tg%i_end', i)), 'String'));
    if get(handles.display_unit_slot, 'Value') ~= 1.0
        tg_start = tg_start * handles.head.Resolution;
        tg_end = tg_end * handles.head.Resolution;
    end
    y = max(handles.display_lifetimedata(:,6)) * 1.1;
    color = {'c--', 'm--', 'r--', 'g--'};
    semilogy([tg_start tg_start], [1 y], color{i}, [tg_end tg_end], [1 y], color{i});
end
xlim(handles.lifetime_zoom_x);
ylim(handles.lifetime_zoom_y);
hold off;
axes(old_context);

function zoom_area_lifetime_Callback(hObject, eventdata, handles)
old_context = gca;
axes(handles.lifetime_axes);
[x y] = ginput(2);
x=sort(x);
y=sort(y);
xlim(x);
ylim(y);
handles.lifetime_zoom_x = x;
handles.lifetime_zoom_y = y;
axes(old_context);
guidata(hObject, handles);

function reset_zoom_lifetime_Callback(hObject, eventdata, handles)
old_context = gca;
axes(handles.lifetime_axes);
xlim('auto');
ylim('auto');
handles.lifetime_zoom_x = 'auto';
handles.lifetime_zoom_y = 'auto';
axes(old_context);
guidata(hObject, handles);

function update_lifetime_calculation_Callback(hObject, eventdata, handles)
handles = recalculate_lifetimedata(handles);
handles = update_display_lifetimedata(handles);
handles = update_lifetimedata_graph(handles);
guidata(hObject, handles);

function display_unit_select_SelectionChangeFcn(hObject, eventdata, handles)
handles = update_display_lifetimedata(handles);
handles = update_lifetimedata_graph(handles);
guidata(hObject, handles);

function tg1_start_Callback(hObject, eventdata, handles)

function tg1_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tg2_start_Callback(hObject, eventdata, handles)

function tg2_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tg3_start_Callback(hObject, eventdata, handles)

function tg3_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tg4_start_Callback(hObject, eventdata, handles)

function tg4_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tg1_end_Callback(hObject, eventdata, handles)

function tg1_end_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tg2_end_Callback(hObject, eventdata, handles)

function tg2_end_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tg3_end_Callback(hObject, eventdata, handles)


function tg3_end_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tg4_end_Callback(hObject, eventdata, handles)

function tg4_end_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tg4_select_Callback(hObject, eventdata, handles)
[x y] = ginput(2);
x=sort(x);
y=sort(y);
if get(handles.display_unit_slot, 'Value') ~= 1.0
    x = x / handles.head.Resolution;
end
x=round(x);
set(handles.tg4_start, 'String', x(1));
set(handles.tg4_end, 'String', x(2));
handles = update_lifetimedata_graph(handles);
guidata(hObject, handles);

function tg3_select_Callback(hObject, eventdata, handles)
[x y] = ginput(2);
x=sort(x);
y=sort(y);
if get(handles.display_unit_slot, 'Value') ~= 1.0
    x = x / handles.head.Resolution;
end
x=round(x);
set(handles.tg3_start, 'String', x(1));
set(handles.tg3_end, 'String', x(2));
handles = update_lifetimedata_graph(handles);
guidata(hObject, handles);

function tg2_select_Callback(hObject, eventdata, handles)
[x y] = ginput(2);
x=sort(x);
y=sort(y);
if get(handles.display_unit_slot, 'Value') ~= 1.0
    x = x / handles.head.Resolution;
end
x=round(x);
set(handles.tg2_start, 'String', x(1));
set(handles.tg2_end, 'String', x(2));
handles = update_lifetimedata_graph(handles);
guidata(hObject, handles);

function tg1_select_Callback(hObject, eventdata, handles)
[x y] = ginput(2);
x=sort(x);
y=sort(y);
if get(handles.display_unit_slot, 'Value') ~= 1.0
    x = x / handles.head.Resolution;
end
x=round(x);
set(handles.tg1_start, 'String', x(1));
set(handles.tg1_end, 'String', x(2));
handles = update_lifetimedata_graph(handles);
guidata(hObject, handles);

function tg1_show_Callback(hObject, eventdata, handles)
handles = update_lifetimedata_graph(handles);
guidata(hObject, handles);

function tg2_show_Callback(hObject, eventdata, handles)
handles = update_lifetimedata_graph(handles);
guidata(hObject, handles);

function tg3_show_Callback(hObject, eventdata, handles)
handles = update_lifetimedata_graph(handles);
guidata(hObject, handles);

function tg4_show_Callback(hObject, eventdata, handles)
handles = update_lifetimedata_graph(handles);
guidata(hObject, handles);

function [corrtype] = fetch_corrtype(handles)
    if get(handles.calc_select_one1, 'Value') > 0
        corrtype=1;
    elseif get(handles.calc_select_one2, 'Value') > 0
        corrtype=2;
    else
        corrtype=3;
    end 
            
function [tg1, tg2, channel] = fetch_timegates(handles)
    tg1 = zeros(1,2);
    tg2 = zeros(16,2);
    channel = zeros(16,1);
    mcs_tg1 = str2double(get(handles.mcs_timegate1_edit, 'String'));
    mcs_tg2 = str2double(get(handles.mcs_timegate2_edit, 'String'));
    tg1(:) = [mcs_tg1 mcs_tg2];
    for i=1:4
        name1 = sprintf('tg%i_start', i);
        name2 = sprintf('tg%i_end', i);
        lifetime_gate_start = str2double(get(handles.(name1), 'String'));
        lifetime_gate_end = str2double(get(handles.(name2), 'String'));

        for j=1:4
            index = (i-1) * 4 + j;
            tg2(index,:) = [lifetime_gate_start lifetime_gate_end];
            channel(index) = j;
        end
    end
    if get(handles.calc_select_one1, 'Value') > 0 || get(handles.calc_select_one2, 'Value') > 0
        % This should quater the calculation time, as only the upper left
        % quadrant of the matrix has to be calculated
        %
        % By setting it to -1 we make sure there is no tcspc value, that can
        % match that    
        tg2(9:16, :) = zeros(8,2)-1;
    end 
    
function [auto, autotime, ltd, head, res] = spawn_part_crosscorrelation(work_fn, handles, tg1_supplied)
% This is a thin wrapper around the calls to part_crosscorraltion
% that were previously done for mass and single calculation

    [tg1, tg2, channel] = fetch_timegates(handles);
    
    if exist('tg1_supplied', 'var')
        tg1 = tg1_supplied;
    end

    options = struct();
    options.statusfield = handles.status_text;
    calc_dc_url = get(handles.calc_dc_url, 'String');
    if not(isempty(calc_dc_url))
        options.dc.lookupurl = calc_dc_url;
    end
    calc_dc_jmname = get(handles.calc_dc_jmname, 'String');
    if not(isempty(calc_dc_jmname))
        options.dc.jmname = calc_dc_jmname;
    end

    calc_dc_config = get(handles.calc_dc_config, 'String');
    if not(isempty(calc_dc_config))
        options.dc.config = calc_dc_config;
    end

    package_size = str2double(get(handles.calc_packagesize, 'String'));
    overlap = str2double(get(handles.calc_overlap, 'String'));
    calc_reverse = get(handles.calc_reverse, 'Value');

    calc_options.Nsub = str2double(get(handles.nsub, 'String'));
    calc_options.max_t = str2double(get(handles.t_max, 'String'));
    
    [auto, autotime, ltd, head, res] = part_crosscorrelation( ...
                                             work_fn, ...
                                             tg1, ...
                                             tg2, ...
                                             channel, ...
                                             package_size, ...
                                             overlap, ...
                                             calc_reverse, ...
                                             options, calc_options);
                                         
function mass_calculation_Callback(hObject, eventdata, handles)
    tic;

    [input_fn, input_path] = uigetfile('*.ht3', 'Select HT3 to process', handles.last_path,'MultiSelect', 'on');

    if isequal(input_fn,0)
        return
    end

    if ~ iscell(input_fn)
        input_fn = {input_fn};
    end

    handles.last_path = input_path;

%    base_path = uigetdir('Destination for correlation data', handles.last_path);
%    if isequal(base_path,0)
%        return
%    end
%
%    handles.last_path = base_path;

    defaults = struct();
    head = ht3v2read_cl(fullfile(input_path, input_fn{1}), [0 0]);
    defaults.Temperature = head.ParamStart(1);
    
    corrtype = fetch_corrtype(handles);
    if corrtype == 3
        metadata = query_metadata(0, defaults, 2);
    else
        metadata = query_metadata(0, defaults, 1);
    end

    if metadata.Temperature < head.ParamStart(1) + 0.1 && metadata.Temperature > head.ParamStart(1) - 0.1
        changedTemp = 0;
    else
        changedTemp = 1;
    end
    
    guidata(hObject, handles);

    files_per_run = str2double(get(handles.files_per_run, 'String'));
    start_index = 1;

    while start_index <= size(input_fn, 2)
        last_index = start_index + files_per_run - 1;
        if last_index > size(input_fn, 2)
            last_index = size(input_fn, 2);
        end

        subset_fn = {input_fn{start_index:last_index}};

        work_fn = cell(size(subset_fn));
        for i=1:size(work_fn,2)
            work_fn{i} = fullfile(input_path, subset_fn{i});
        end

        [auto, atime, ltd, head, res] = spawn_part_crosscorrelation(work_fn, handles);

        for i=1:size(subset_fn,2)
            in_fn = subset_fn{i};
            out_fn = in_fn(1:end-4);
            res = res{i};
            head = head{i};
            res.autotime = atime{:,:};
            [tg1, tg2] = fetch_timegates(handles);
            res.metadata = metadata;
            res.metadata.timegate_mcs = tg1;
            res.metadata.timegate_lt = tg2;
            res.metadata.corrtype = corrtype;
            res.metadata.correlation = {auto{i,:}};

            save(fullfile(input_path, out_fn), 'head', 'res', 'informations');
        end
        start_index = last_index + 1;
    end
    set(handles.status_text, 'String', 'Work done!');
    times = int32(toc);
    hh = num2str(fix(times/3600));
    min = num2str(fix((times-(fix(times/3600)*3600))/60));
    if size(min) == 1
       min = ['0' min]; 
    end
    sec = num2str(fix(times-((fix(times/3600)*3600)+(fix((times-(fix(times/3600)*3600))/60)*60))));
    if size(sec) == 1
       sec = ['0' sec]; 
    end 
    zeit = [num2str(hh) ':' num2str(min) ':' num2str(sec)];
    fprintf(['Calculation took: ' zeit ' hh:mm:ss' '\n']);
            
function single_calculation_Callback(hObject, eventdata, handles)
    tic;
    
    %[base_fn, base_path] = uiputfile('*.xcd; *.ltd2', 'Destination for correlation data', handles.last_path);
    %files are stored in same folder than ht3-files
    [base_path, base_fn] = fileparts(handles.pt3fn);
    
    if isequal(base_fn,0) || isequal(base_path,0)
        return
    end

    if strcmp(base_fn(1,end-3:end), '.xcd')
       base_fn =base_fn(1,1:end-4);
    end

    defaults = struct();
    [head] = ht3v2read_cl(handles.pt3fn, [0 0]);
    defaults.Temperature = head.ParamStart(1);
    
    corrtype = fetch_corrtype(handles);
    if corrtype == 3
        metadata = query_metadata(0, defaults, 2);
    else
        metadata = query_metadata(0, defaults, 1);
    end
    
    handles.last_path = base_path;

    guidata(hObject, handles);

    [auto, atime, ltd, head, res] = spawn_part_crosscorrelation({handles.pt3fn}, handles);
   
    corrtype = fetch_corrtype(handles);
    [tg1, tg2] = fetch_timegates(handles);
    res = res{1};
    res.autotime = atime{:,:};
    head = head{1};
    res.metadata = metadata;
    res.metadata.timegate_mcs = tg1;
    res.metadata.timegate_lt = tg2;
    res.metadata.corrtype = corrtype;
    res.metadata.correlation = {auto{1,:}};
    save(fullfile(base_path, base_fn), 'head', 'res');

    set(handles.status_text, 'String', 'Work done!');
    times = int32(toc);
    hh = num2str(fix(times/3600));
    min = num2str(fix((times-(fix(times/3600)*3600))/60));
    if size(min) == 1
       min = ['0' min]; 
    end
    sec = num2str(fix(times-((fix(times/3600)*3600)+(fix((times-(fix(times/3600)*3600))/60)*60))));
    if size(sec) == 1
       sec = ['0' sec]; 
    end 
    zeit = [num2str(hh) ':' num2str(min) ':' num2str(sec)];
    fprintf(['Calculation took: ' zeit ' hh:mm:ss' '\n']);

function calc_packagesize_Callback(hObject, eventdata, handles)
function calc_packagesize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function calc_overlap_Callback(hObject, eventdata, handles)
function calc_overlap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function calc_reverse_Callback(hObject, eventdata, handles)

function t_max_Callback(hObject, eventdata, handles)
function t_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nsub_Callback(hObject, eventdata, handles)
function nsub_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function files_per_run_Callback(hObject, eventdata, handles)
function files_per_run_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function calc_dc_config_Callback(hObject, eventdata, handles)
function calc_dc_config_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function calc_dc_url_Callback(hObject, eventdata, handles)
function calc_dc_url_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function calc_dc_jmname_Callback(hObject, eventdata, handles)
function calc_dc_jmname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tg_autodetect_4_Callback(hObject, eventdata, handles)
tg_autodetect(hObject, handles, 4);

function tg_autodetect_2_Callback(hObject, eventdata, handles)
tg_autodetect(hObject, handles, 2);

function tg_autodetect(hObject, handles, count)
bounds = detect_peaks(handles.display_lifetimedata(:,1), handles.display_lifetimedata(:,6), count);
if get(handles.display_unit_slot, 'Value') ~= 1.0
    bounds = bounds / handles.head.Resolution;
end
if size(bounds,1) < 4
    for i=size(bounds,1)+1:4
        bounds(i,1) = 0;
        bounds(i,2) = 65536;
    end
end
for i=1:4
    set(handles.(sprintf('tg%d_start', i)), 'String', bounds(i, 1));
    set(handles.(sprintf('tg%d_end', i)), 'String', bounds(i,2));
end
handles = update_lifetimedata_graph(handles);
guidata(hObject, handles);

function multiple_mcs_Callback(hObject, eventdata, handles)
    tic;
    
    mcs_answer = inputdlg('Specify ranges for MCS in seconds. One pair per line. Two values seperated by space.', 'MCS Selection', 25);
    
    if size(mcs_answer) == 0
        return
    end
    
    mcs_answer = mcs_answer{1};
    mcs_selection = zeros(0,2);
    
    for i=1:size(mcs_answer, 1)
        bounds = sscanf(mcs_answer(i, :), '%f');
        if size(bounds, 1) == 2
            mcs_selection(size(mcs_selection, 1)+1, :) = bounds';
        end
    end
    
    [base_fn, base_path] = uiputfile('*.xcd; *.ltd2', 'Destination for correlation data', handles.last_path);
    if isequal(base_fn,0) || isequal(base_path,0)
        return
    end

    if strcmp(base_fn(1,end-3:end), '.xcd')
        base_fn =base_fn(1,1:end-4);
    end

    defaults = struct();
    [head] = ht3v2read_cl(handles.pt3fn, [0 0]);
    defaults.Temperature = head.ParamStart(1);
    
    corrtype = fetch_corrtype(handles);
    if corrtype == 3
        metadata = query_metadata(0, defaults, 2);
    else
        metadata = query_metadata(0, defaults, 1);
    end
    
    handles.last_path = base_path;

    guidata(hObject, handles);

    for i=1:size(mcs_selection, 1)
        corrtype = fetch_corrtype(handles);
        [tg1, tg2] = fetch_timegates(handles);
        tg1 = mcs_selection(i, :);
        [auto, atime, ltd, head] = spawn_part_crosscorrelation({handles.pt3fn}, handles, tg1);
        lifetime_data = ltd{1};
        head_data = head{1};
        autotime = atime{1};
        correlation = {auto{1,:}};
        fn = sprintf('%s-mmcs-%i', base_fn, i);
        save(fullfile(base_path, fn), 'lifetime_data', 'head_data', 'metadata', 'tg1', 'tg2', 'correlation', 'autotime', 'corrtype');
        tg2_export = zeros(4,2);
        for j=1:4:size(tg2,1)
            tg2_export(floor(j/4)+1,:) = tg2(j,:);
        end
        export_lifetimedata(base_path, fn, head_data, metadata, tg1, tg2_export, lifetime_data);
        export_crosscorrelation(base_path, fn, head_data, metadata, tg1, tg2_export, autotime, correlation, corrtype);
    end
    set(handles.status_text, 'String', 'Work done!');
    times = int32(toc);
    hh = num2str(fix(times/3600));
    min = num2str(fix((times-(fix(times/3600)*3600))/60));
    if size(min) == 1
       min = ['0' min]; 
    end
    sec = num2str(fix(times-((fix(times/3600)*3600)+(fix((times-(fix(times/3600)*3600))/60)*60))));
    if size(sec) == 1
       sec = ['0' sec]; 
    end 
    zeit = [num2str(hh) ':' num2str(min) ':' num2str(sec)];
    fprintf(['Calculation took: ' zeit ' hh:mm:ss' '\n']);
