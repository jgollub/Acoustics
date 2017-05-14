function varargout = Newmark2D_AcousticScanGUI(varargin)
% NEWMARK2D_NEARFIELDSCANGUI MATLAB code for Newmark2D_NearFieldScanGUI.fig
%      NEWMARK2D_NEARFIELDSCANGUI, by itself, creates a new NEWMARK2D_NEARFIELDSCANGUI or raises the existing
%      singleton*.
%
%      H = NEWMARK2D_NEARFIELDSCANGUI returns the handle to a new NEWMARK2D_NEARFIELDSCANGUI or the handle to
%      the existing singleton*.
%
%      NEWMARK2D_NEARFIELDSCANGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEWMARK2D_NEARFIELDSCANGUI.M with the given input arguments.
%
%      NEWMARK2D_NEARFIELDSCANGUI('Property','Value',...) creates a new NEWMARK2D_NEARFIELDSCANGUI or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Newmark2D_NearFieldScanGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Newmark2D_NearFieldScanGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Newmark2D_NearFieldScanGUI

% Last Modified by GUIDE v2.5 24-Aug-2013 00:32:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Newmark2D_NearFieldScanGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Newmark2D_NearFieldScanGUI_OutputFcn, ...
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

% --- Executes just before Newmark2D_NearFieldScanGUI is made visible.
function Newmark2D_NearFieldScanGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Newmark2D_NearFieldScanGUI (see VARARGIN)

% Choose default command line output for Newmark2D_NearFieldScanGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes Newmark2D_NearFieldScanGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the abortpushbutton flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to abortpushbutton the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

%initialize stage communication
handles.Newmark2D_stage=Newmark2D_stage_start;
handles.steptomm=5000;

%set default stage range
handles.xStageRange=1000.6; %in mm
set(handles.tagStageRangeX,'String',num2str(handles.xStageRange))
handles.yStageRange=999.1672;
set(handles.tagStageRangeY,'String',num2str(handles.yStageRange))

%arbitrary but should be roughly in the center
handles.defZeroInXsteps=2490013; 
handles.defZeroInYsteps=2490013;

%set other variables
handles.metricdata.moverelativeXPos = 0;
handles.metricdata.moverelativeYPos = 0;

handles.metricdata.motorspeed = 25;
set(handles.motorspeed,'String',num2str(handles.metricdata.motorspeed))

handles.metricdata.xrangemin = -(11+2)*25.4;
handles.metricdata.xrangemax = (11+2)*25.4;
handles.metricdata.yrangemin = -(8+2)*25.4;
handles.metricdata.yrangemax = (8+2)*25.4;
handles.metricdata.stepsize = 5;
handles.metricdata.fsamples = 401;
handles.metricdata.fstart = 17.5;
handles.metricdata.fstop = 26.5;
handles.metricdata.sports = 0;
handles.metricdata.buttonState=1;

handles.metricdata.IFbandInHZ=1000;

%handles.metricdata.calfile='SOLT_Cal_S31_S32_S21_p3Coupler_Off_6ft_13dBm_07082013';
%handles.metricdata.calfile='ECal_3port_3ftpanA_4ftpanB_9ftHorn082413';
handles.metricdata.calfile='ECAL_P123_Coupler_off_9ftToHorn_-5dBm_09172013'
% handles.metricdata.calfile='Ecal_3port_3ftpanA_4ftPanB_9ftULLHorn083013';
%handles.metricdata.calfile='ResponsThru_Cal_S31andS32_B_3ftpanA_9ftULLHorn_10dBm_08282013';

% [currPosX, currPosY] = Newmark2D_stage_getPosition(handles.Newmark2D_stage,handles.defZeroInXsteps,handles.defZeroInYsteps);
%  
% set(handles.currentXpos, 'String', num2str(currPosX));
% set(handles.currentYpos, 'String', num2str(currPosY));
set(handles.currentXpos, 'String', 'NaN');
set(handles.currentYpos, 'String', 'NaN');

set(handles.xaxismovedistance, 'String', '0');
set(handles.yaxismovedistance, 'String', '0');
set(handles.xrangemax, 'String', num2str(handles.metricdata.xrangemax));
set(handles.xrangemin, 'String', num2str(handles.metricdata.xrangemin));
set(handles.yrangemax, 'String', num2str(handles.metricdata.yrangemax));
set(handles.yrangemin, 'String', num2str(handles.metricdata.yrangemin));
set(handles.stepsize, 'String', num2str(handles.metricdata.stepsize));
set(handles.fsamples, 'String', num2str(handles.metricdata.fsamples));
set(handles.fstop, 'String', num2str(handles.metricdata.fstop));
set(handles.fstart, 'String', num2str(handles.metricdata.fstart));
set(handles.switchports, 'String', num2str(handles.metricdata.sports));

set(handles.IFbandwidthTag, 'String', num2str(handles.metricdata.IFbandInHZ));
set(handles.Calibration_File, 'String', handles.metricdata.calfile);


% Update handles structure
guidata(handles.figure1, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Newmark2D_NearFieldScanGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% handles.metricdata.X = X;
% handles.metricdata.Y = Y;
% handles.metricdata.f = f;
% handles.metricdata.measurements = measurements;
varargout{1} = handles.output;


% --- Executes on button press in startscanpushbutton.
function startscanpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to startscanpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stagespeed = handles.metricdata.motorspeed;
xmin = handles.metricdata.xrangemin;
xmax = handles.metricdata.xrangemax;
ymin = handles.metricdata.yrangemin;
ymax = handles.metricdata.yrangemax;
dstep = handles.metricdata.stepsize;
frequencysamples = handles.metricdata.fsamples;
fstart = handles.metricdata.fstart;
fstop = handles.metricdata.fstop;
switches = handles.metricdata.sports;
IFBandwidth=handles.metricdata.IFbandInHZ;
calfile=handles.metricdata.calfile;
measureSparam=handles.metricdata.buttonState;


%COMMENTED OUT BY TIM - also moved "data" block and save inside _fcn
[data] = Newmark2D_acoustic_fcn(handles.Newmark2D_stage,stagespeed,...
    handles.defZeroInXsteps,handles.defZeroInYsteps,xmin,xmax,ymin,ymax,dstep,frequencysamples,...
    fstart,fstop,switches,IFBandwidth,calfile,measureSparam);
%COMMENTED OUT BY TIM - also moved "data" block and save inside _fcn


% --- Executes on button press in abortpushbutton.
function abortpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to abortpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Newmark2D_stage_StopStage(handles.Newmark2D_stage)
initialize_gui(gcbf, handles, true);

% --- Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (hObject == handles.english)
    set(handles.text4, 'String', 'lb/cu.in');
    set(handles.text5, 'String', 'cu.in');
    set(handles.text6, 'String', 'lb');
else
    set(handles.text4, 'String', 'kg/cu.m');
    set(handles.text5, 'String', 'cu.m');
    set(handles.text6, 'String', 'kg');
end

% --- Executes on button press in moverelativepushbutton.
function moverelativepushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to moverelativepushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Newmark2D_stage_moveToRelative(handles.Newmark2D_stage,handles.metricdata.motorspeed,...
     handles.defZeroInXsteps, handles.defZeroInYsteps,handles.metricdata.moverelativeXPos,...
     handles.metricdata.moverelativeYPos);
[x, y] = Newmark2D_stage_getPosition(handles.Newmark2D_stage,handles.defZeroInXsteps,handles.defZeroInYsteps);
set(handles.currentXpos, 'String', num2str(x));
set(handles.currentYpos, 'String', num2str(y));


%show range of stage available
XegdeNeg=-(x+handles.xStageRange-(handles.defZeroInXsteps-handles.xHomeInSteps)/handles.steptomm);
XegdePos=-(x-(handles.defZeroInXsteps-handles.xHomeInSteps)/handles.steptomm);
set(handles.x_StE_Neg, 'String', num2str(XegdeNeg));
set(handles.x_StE_Pos, 'String', num2str(XegdePos));

YegdeNeg=-(y+handles.yStageRange-(handles.defZeroInYsteps-handles.yHomeInSteps)/handles.steptomm);
YegdePos=-(y-(handles.defZeroInYsteps-handles.yHomeInSteps)/handles.steptomm);
set(handles.y_StE_Neg, 'String', num2str(YegdeNeg));
set(handles.y_StE_Pos, 'String', num2str(YegdePos));

%get step position and post to GUI
[xposSteps,yposSteps] = Newmark2D_stage_getStepPosition(handles.Newmark2D_stage)
set(handles.tagstepXpos, 'String', num2str(xposSteps));
set(handles.tagstepYpos, 'String', num2str(yposSteps));

function xaxismovedistance_Callback(hObject, eventdata, handles)
% hObject    handle to xaxismovedistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xaxismovedistance as text
%        str2double(get(hObject,'String')) returns contents of xaxismovedistance as a double
moverelativeXPos = str2double(get(hObject, 'String'));
if isnan(moverelativeXPos)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new moverelativeXPos value
handles.metricdata.moverelativeXPos = moverelativeXPos;
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function xaxismovedistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xaxismovedistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function yaxismovedistance_Callback(hObject, eventdata, handles)
% hObject    handle to yaxismovedistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yaxismovedistance as text
%        str2double(get(hObject,'String')) returns contents of yaxismovedistance as a double
moverelativeYPos = str2double(get(hObject, 'String'));
if isnan(moverelativeYPos)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new moverelativeYPos value
handles.metricdata.moverelativeYPos = moverelativeYPos;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function yaxismovedistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yaxismovedistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in zeroaxes.
function zeroaxes_Callback(hObject, eventdata, handles)
% hObject    handle to zeroaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.defZeroInXsteps, handles.defZeroInYsteps]=Newmark2D_stage_zeroAxes(handles.Newmark2D_stage);

% fprintf(['position X set to =\n', handles.defZeroInXsteps]);
% fprintf(['position Y set to =\n', handles.defZeroInYsteps]);

 set(handles.currentXpos, 'String', '0');
 set(handles.currentYpos, 'String', '0');

guidata(hObject,handles)

% --- Executes on button press in movetozerobutton.
function movetozerobutton_Callback(hObject, eventdata, handles)
% hObject    handle to movetozerobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Newmark2D_stage_moveToAbsolute(handles.Newmark2D_stage,handles.metricdata.motorspeed,...
    handles.defZeroInXsteps, handles.defZeroInYsteps, 0, 0);

[x, y] = Newmark2D_stage_getPosition(handles.Newmark2D_stage,handles.defZeroInXsteps,handles.defZeroInYsteps);

set(handles.currentXpos, 'String', num2str(x));
set(handles.currentYpos, 'String', num2str(y));

%show range of stage available
XegdeNeg=-(x+handles.xStageRange-(handles.defZeroInXsteps-handles.xHomeInSteps)/handles.steptomm);
XegdePos=-(x-(handles.defZeroInXsteps-handles.xHomeInSteps)/handles.steptomm);
set(handles.x_StE_Neg, 'String', num2str(XegdeNeg));
set(handles.x_StE_Pos, 'String', num2str(XegdePos));

YegdeNeg=-(y+handles.yStageRange-(handles.defZeroInYsteps-handles.yHomeInSteps)/handles.steptomm);
YegdePos=-(y-(handles.defZeroInYsteps-handles.yHomeInSteps)/handles.steptomm);
set(handles.y_StE_Neg, 'String', num2str(YegdeNeg));
set(handles.y_StE_Pos, 'String', num2str(YegdePos));

%get step position and post to GUI
[xposSteps,yposSteps] = Newmark2D_stage_getStepPosition(handles.Newmark2D_stage);
set(handles.tagstepXpos, 'String', num2str(xposSteps));
set(handles.tagstepYpos, 'String', num2str(yposSteps));
% Update handles structure
guidata(hObject,handles);

function motorspeed_Callback(hObject, eventdata, handles)
% hObject    handle to motorspeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of motorspeed as text
%        str2double(get(hObject,'String')) returns contents of motorspeed as a double
speed = str2double(get(hObject, 'String'));
if isnan(speed)
    set(hObject, 'String', 25);
    errordlg('Input must be a number','Error');
end
if speed>50
    set(hObject, 'String', 50);
    speed=50;
    errordlg('25mm/s is the max speed','Error');
end
if speed<=0
    set(hObject, 'String', 25);
    speed=25;
    errordlg('Speed must be greater than 0mm/s','Error');
end
%set(handles.motorspeed, 'String', speed);
% Save the new motorspeed value
handles.metricdata.motorspeed = speed;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function motorspeed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to motorspeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in moveabsolutebutton.
function moveabsolutebutton_Callback(hObject, eventdata, handles)
% hObject    handle to moveabsolutebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Newmark2D_stage_moveToAbsolute(handles.Newmark2D_stage,handles.metricdata.motorspeed,...
    handles.defZeroInXsteps, handles.defZeroInYsteps,...
    handles.metricdata.moverelativeXPos,handles.metricdata.moverelativeYPos);

[x, y] = Newmark2D_stage_getPosition(handles.Newmark2D_stage,handles.defZeroInXsteps,handles.defZeroInYsteps);;
set(handles.currentXpos, 'String', num2str(x));
set(handles.currentYpos, 'String', num2str(y));

%show range of stage available
XegdeNeg=-(x+handles.xStageRange-(handles.defZeroInXsteps-handles.xHomeInSteps)/handles.steptomm);
XegdePos=-(x-(handles.defZeroInXsteps-handles.xHomeInSteps)/handles.steptomm);
set(handles.x_StE_Neg, 'String', num2str(XegdeNeg));
set(handles.x_StE_Pos, 'String', num2str(XegdePos));

YegdeNeg=-(y+handles.yStageRange-(handles.defZeroInYsteps-handles.yHomeInSteps)/handles.steptomm);
YegdePos=-(y-(handles.defZeroInYsteps-handles.yHomeInSteps)/handles.steptomm);
set(handles.y_StE_Neg, 'String', num2str(YegdeNeg));
set(handles.y_StE_Pos, 'String', num2str(YegdePos));

%get step position and post to GUI
[xposSteps,yposSteps] = Newmark2D_stage_getStepPosition(handles.Newmark2D_stage)
set(handles.tagstepXpos, 'String', num2str(xposSteps));
set(handles.tagstepYpos, 'String', num2str(yposSteps));

function xrangemin_Callback(hObject, eventdata, handles)
% hObject    handle to xrangemin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xrangemin as text
%        str2double(get(hObject,'String')) returns contents of xrangemin as a double
xmin = str2double(get(hObject, 'String'));
if isnan(xmin)
    set(hObject, 'String', 20);
    errordlg('Input must be a number','Error');
end
%set(handles.motorspeed, 'String', speed);
% Save the new motorspeed value
handles.metricdata.xrangemin = xmin;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function xrangemin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xrangemin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xrangemax_Callback(hObject, eventdata, handles)
% hObject    handle to xrangemax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xrangemax as text
%        str2double(get(hObject,'String')) returns contents of xrangemax as a double
xmax = str2double(get(hObject, 'String'));
if isnan(xmax)
    set(hObject, 'String', 20);
    errordlg('Input must be a number','Error');
end
%set(handles.motorspeed, 'String', speed);
% Save the new motorspeed value
handles.metricdata.xrangemax = xmax;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function xrangemax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xrangemax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yrangemin_Callback(hObject, eventdata, handles)
% hObject    handle to yrangemin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yrangemin as text
%        str2double(get(hObject,'String')) returns contents of yrangemin as a double
ymin = str2double(get(hObject, 'String'));
if isnan(ymin)
    set(hObject, 'String', 20);
    errordlg('Input must be a number','Error');
end
%set(handles.motorspeed, 'String', speed);
% Save the new motorspeed value
handles.metricdata.yrangemin = ymin;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function yrangemin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yrangemin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yrangemax_Callback(hObject, eventdata, handles)
% hObject    handle to yrangemax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yrangemax as text
%        str2double(get(hObject,'String')) returns contents of yrangemax as a double
ymax = str2double(get(hObject, 'String'));
if isnan(ymax)
    set(hObject, 'String', 20);
    errordlg('Input must be a number','Error');
end
%set(handles.motorspeed, 'String', speed);
% Save the new motorspeed value
handles.metricdata.yrangemax = ymax;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function yrangemax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yrangemax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function stepsize_Callback(hObject, eventdata, handles)
% hObject    handle to stepsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stepsize as text
%        str2double(get(hObject,'String')) returns contents of stepsize as a double
step = str2double(get(hObject, 'String'));
if isnan(step)
    set(hObject, 'String', 20);
    errordlg('Input must be a number','Error');
end
%set(handles.motorspeed, 'String', speed);
% Save the new motorspeed value
handles.metricdata.stepsize = step;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function stepsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stepsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fsamples_Callback(hObject, eventdata, handles)
% hObject    handle to fsamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fsamples as text
%        str2double(get(hObject,'String')) returns contents of fsamples as a double
fsamples = str2double(get(hObject, 'String'));
if isnan(fsamples)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new moverelativeXPos value
handles.metricdata.fsamples= fsamples;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function fsamples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fsamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fstart_Callback(hObject, eventdata, handles)
% hObject    handle to fstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fstart as text
%        str2double(get(hObject,'String')) returns contents of fstart as a double
fstart = str2double(get(hObject, 'String'));
if isnan(fstart)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new moverelativeXPos value
handles.metricdata.fstart = fstart;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function fstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fstop_Callback(hObject, eventdata, handles)
% hObject    handle to fstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fstop as text
%        str2double(get(hObject,'String')) returns contents of fstop as a double
fstop = str2double(get(hObject, 'String'));
if isnan(fstop)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new moverelativeXPos value
handles.metricdata.fstop = fstop;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function fstop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function switchports_Callback(hObject, eventdata, handles)
% hObject    handle to switchports (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of switchports as text
%        str2double(get(hObject,'String')) returns contents of switchports as a double
sports = str2double(get(hObject, 'String'));
if isnan(sports)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new moverelativeXPos value
handles.metricdata.sports = sports;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function switchports_CreateFcn(hObject, eventdata, handles)
% hObject    handle to switchports (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IFbandwidthTag_Callback(hObject, eventdata, handles)
% hObject    handle to IFbandwidthTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
IFbandwidth=str2double(get(hObject, 'String'));
if isnan(IFbandwidth)
    set(hObject, 'String', 1000);
    errordlg('Input must be a number','Error');
end
handles.metricdata.IFbandInHZ=IFbandwidth
guidata(hObject,handles)

% Hints: get(hObject,'String') returns contents of IFbandwidthTag as text
%        str2double(get(hObject,'String')) returns contents of IFbandwidthTag as a double


% --- Executes during object creation, after setting all properties.
function IFbandwidthTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IFbandwidthTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Meas_Sparam_tag.
function Meas_Sparam_tag_Callback(hObject, eventdata, handles)
% hObject    handle to Meas_Sparam_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Meas_Sparam_tag


% --- Executes during object creation, after setting all properties.
function Meas_Sparam_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Meas_Sparam_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on key press with focus on Meas_Sparam_tag and none of its controls.
function Meas_Sparam_tag_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Meas_Sparam_tag (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function Meas_Sparam_tag_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to Meas_Sparam_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function Calibration_File_Callback(hObject, eventdata, handles)
% hObject    handle to Calibration_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Calibration_File as text
%        str2double(get(hObject,'String')) returns contents of Calibration_File as a double
calfile=get(hObject, 'String');
handles.metricdata.calfile=calfile;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Calibration_File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Calibration_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text16.
function text16_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiobuttonS31.
function radiobuttonS31_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonS31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
buttonState=get(hObject,'Value');
if buttonState==get(hObject,'Max')
    handles.metricdata.buttonState=1;
end
  guidata(hObject,handles)  
% Hint: get(hObject,'Value') returns toggle state of radiobuttonS31


% --- Executes on button press in radiobuttonS23.
function radiobuttonS23_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonS23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
buttonState=get(hObject,'Value');
if buttonState==get(hObject,'Max')
    handles.metricdata.buttonState=2;
end
  guidata(hObject,handles)  
% Hint: get(hObject,'Value') returns toggle state of radiobuttonS23


% --- Executes on button press in radiobuttonS21.
function radiobuttonS21_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonS21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
buttonState=get(hObject,'Value');
if buttonState==get(hObject,'Max')
    handles.metricdata.buttonState=3;
end
  guidata(hObject,handles)  
% Hint: get(hObject,'Value') returns toggle state of radiobuttonS21


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
keyboard

% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function movetozerobutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to movetozerobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in homeStage.
function homeStage_Callback(hObject, eventdata, handles)
% hObject    handle to homeStage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.xHomeInSteps, handles.yHomeInSteps] = Newmark2D_stage_Home(handles.Newmark2D_stage,handles.metricdata.motorspeed);

set(handles.x_StE_Neg, 'Visible','on');
set(handles.x_StE_Pos, 'Visible','on');
XegdeNeg=(handles.xStageRange);
XegdePos=0;
set(handles.x_StE_Neg, 'String', num2str(XegdeNeg));
set(handles.x_StE_Pos, 'String', num2str(XegdePos));

set(handles.y_StE_Neg, 'Visible','on');
set(handles.y_StE_Pos, 'Visible','on');
YegdeNeg=(handles.yStageRange);
YegdePos=0;
set(handles.y_StE_Neg, 'String', num2str(YegdeNeg));
set(handles.y_StE_Pos, 'String', num2str(YegdePos));

%get step position and post to GUI
[xposSteps,yposSteps]=Newmark2D_stage_getStepPosition(handles.Newmark2D_stage);
set(handles.tagstepXpos,'String',num2str(xposSteps));
set(handles.tagstepYpos,'String',num2str(yposSteps));

guidata(hObject,handles)


% --- Executes on button press in stageRange.
function stageRange_Callback(hObject, eventdata, handles)
% hObject    handle to stageRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.xStageRange,handles.yStageRange] = Newmark2D_stage_RangeOfStage(handles.Newmark2D_stage,handles.metricdata.motorspeed)

%show stage range
[xposSteps,yposSteps]=Newmark2D_stage_getStepPosition(handles.Newmark2D_stage)
set(handles.tagstepXpos,'String',num2str(xposSteps))
set(handles.tagstepYpos,'String',num2str(yposSteps))
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function x_StE_Neg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_StE_Neg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function x_StE_Pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_StE_Pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function currentXpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentXpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
