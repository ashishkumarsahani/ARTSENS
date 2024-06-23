function varargout = User_Interface(varargin)
% USER_INTERFACE MATLAB code for User_Interface.fig
%      USER_INTERFACE, by itself, creates a new USER_INTERFACE or raises the existing
%      singleton*.
%
%      H = USER_INTERFACE returns the handle to a new USER_INTERFACE or the handle to
%      the existing singleton*.
%
%      USER_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in USER_INTERFACE.M with the given input arguments.
%
%      USER_INTERFACE('Property','Value',...) creates a new USER_INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before User_Interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to User_Interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help User_Interface

% Last Modified by GUIDE v2.5 29-Jan-2014 12:58:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @User_Interface_OpeningFcn, ...
                   'gui_OutputFcn',  @User_Interface_OutputFcn, ...
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


% --- Executes just before User_Interface is made visible.
function User_Interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to User_Interface (see VARARGIN)

% Choose default command line output for User_Interface
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes User_Interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global US_Sampling_Rate;
global samples_Per_Frame;
global transducer_Frequency;
global ECG_Sampling_Rate;
global live;
live = 0

% --- Outputs from this function are returned to the command line.
function varargout = User_Interface_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in start_But.
function start_But_Callback(hObject, eventdata, handles)
% hObject    handle to start_But (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global stop_All;

last2SecData = zeros(1400,1);
niScope_Link(1,5195,'Dev3');
%ECG_Link(1);
ylim(handles.US_Frame_Plot,[-4 4]);
ylim(handles.ECG_Plot,[-200 200]);
hold(handles.US_Frame_Plot, 'on')
hold(handles.ECG_Plot, 'on')
%hold on;
stop_All=0;
tic
for i = 1:1000
    [frame,current_Time_Stamp] = niScope_Link(2,5195);
    current_Time_Stamp = uint32(current_Time_Stamp);
    time_Axis = [1:5195];
    
%     [A,B,C] = ECG_Link(2);
%     if(A>0)
%         last2SecData = [last2SecData(A+1:end,1);C];
%         cla(handles.ECG_Plot);
%         plot(handles.ECG_Plot,(1:1400),last2SecData);
%     end
    if(mod(i,5)==0)
        cla(handles.US_Frame_Plot);
        plot(handles.US_Frame_Plot,time_Axis,frame);
        drawnow;        
    end
    if(stop_All ==1)
        break; 
    end
end
toc
i
niScope_Link(3,5195);


% --- Executes on button press in File_Sel_But.
function File_Sel_But_Callback(hObject, eventdata, handles)
% hObject    handle to File_Sel_But (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global US_Sampling_Rate;
global samples_Per_Frame;
global transducer_Frequency;
global ECG_Sampling_Rate;

[FileName,PathName] = uigetfile('*.us;*.ust;*.lvsngl','Select the ARTSENS Record');
set(handles.US_File_Addr,'String',strcat(PathName,FileName));

if(strcmp(FileName(end-5:end),'lvsngl')==1)
    US_Sampling_Rate = 100000000;
    samples_Per_Frame = 5195;
    transducer_Frequency = 500000;
    ECG_Sampling_Rate = 0;
end

set(handles.US_Fs_Text,'String',num2str(US_Sampling_Rate));
set(handles.SPF_Text,'String',num2str(samples_Per_Frame));
set(handles.TF_Text,'String',num2str(transducer_Frequency));
set(handles.ECG_Fs_Text,'String',num2str(ECG_Sampling_Rate));


% --- Executes on button press in Analyze_But.
function Analyze_But_Callback(hObject, eventdata, handles)
% hObject    handle to Analyze_But (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global US_Sampling_Rate;
global samples_Per_Frame;
global transducer_Frequency;
global ECG_Sampling_Rate;
global stop_All;
global Avg_Filt_mat;
global shift_Mat;
global last_Track_Frame;
global points_Per_MM;
global current_State;
global measurement_Count;
global dia_Store;
global dist_Store;
global time_Stamp_Store;
global data_Length_Store;
global D_d_Store;
global del_D_Store;
global current_Display_Record_Index;
global live;

points_Per_MM=130;
storage_Mat_Column_Nos = 5;
Avg_Filt_mat = zeros(samples_Per_Frame,storage_Mat_Column_Nos);
shift_Mat= (gen_Shift_Mat(storage_Mat_Column_Nos));

stop_All=0;


total_Frames=0;
US_Fid =0;
measurement_Count = 0;
dia_Store = zeros(100,1);
dist_Store = zeros(100,1);
D_d_Store = 0;
del_D_Store = 0;
if (live == 0)
    current_File_Name = get(handles.US_File_Addr,'String');
    if(strcmp(current_File_Name(end-5:end),'lvsngl')==1 || strcmp(current_File_Name(end-5:end),'us')==1 || strcmp(current_File_Name(end-5:end),'ust')==1)
        [US_Fid,fail_msg] = fopen(current_File_Name);
        if(US_Fid== -1)
            disp(fail_msg)
            return
        end
    else
        disp('invalid file ...');
        return
    end
    
    fseek(US_Fid, 0,'eof');
    filelength = ftell(US_Fid);
    total_Frames = filelength/(4*samples_Per_Frame);
else                                    %live data scope init
%     ECG_Link(1);
  total_Frames = 5000;
    niScope_Link(1,5195,'Dev3');
    current_File_Name = 'live_data';
end
    
    
if(strcmp(current_File_Name(end-5:end),'lvsngl')==1 || live == 1)
    state1 = 0;
    state2 = 0;


    frame_Num=1;
    tic
    set(handles.US_Frame_Plot,'ALimMode','manual');
    ylim(handles.US_Frame_Plot,[-4 4]);
    xlim(handles.US_Frame_Plot,[1 samples_Per_Frame]);
    hold (handles.dist_Axes, 'on');
    hold(handles.US_Frame_Plot, 'on');
    hold (handles.mmt_Store, 'on');
    prev_Time_Stamp=0;
    while(frame_Num < total_Frames)
        if (live== 1)
            [current_Frame,current_Time_Stamp] = niScope_Link(2,5195);
            Elapsed_time = current_Time_Stamp - prev_Time_Stamp
            prev_Time_Stamp = current_Time_Stamp;
            
        else     % file read if live is not equal to 1
            current_Frame = load_LV_Frame(US_Fid, samples_Per_Frame,frame_Num);
        end
        
        if(state1 <2)
            [state1,wall1,wall2] = find_Artery_Walls_SM (current_Frame,frame_Num);
           
        else
            [state1, wall1, wall2, time_Arr, dist_Arr, dia_Arr, full_Cycle_Time_Array, full_Cycle_Dist_Array, full_Cycle_Lumen_Dia_Array, del_D, D_d] = track_Walls_no_ECG_no_TS(current_Frame, wall1, wall2, 50);
            cla(handles.dist_Axes);
            
           
            plot(handles.dist_Axes,time_Arr, dist_Arr,'r');
            plot(handles.dist_Axes,time_Arr, dia_Arr-mean(dia_Arr), 'g');
            if(state1==5)
               [state1, wall1, wall2] = find_Artery_Walls_SM(current_Frame, 1);
            elseif(state1 ==4)
               disp('one cycle extracted');
               measurement_Count = measurement_Count+1;
               data_Length_Store(measurement_Count) = length(full_Cycle_Dist_Array);
               time_Stamp_Store([1:length(full_Cycle_Time_Array)],measurement_Count) = full_Cycle_Time_Array;
               dia_Store([1:length(full_Cycle_Lumen_Dia_Array)],measurement_Count) = full_Cycle_Lumen_Dia_Array;
               dist_Store([1:length(full_Cycle_Dist_Array)],measurement_Count) = full_Cycle_Dist_Array;
               dist_Store(:,measurement_Count)= dist_Store(:,measurement_Count);
               D_d_Store(measurement_Count) = D_d;
               del_D_Store(measurement_Count) = del_D;
               current_Display_Record_Index = measurement_Count;
               
               cla(handles.mmt_Store);
               plot(handles.mmt_Store,dist_Store((1:data_Length_Store(measurement_Count) ),measurement_Count),'r');
               current_State =3;
            end
        end
        if(mod(frame_Num,5)==0)
            cla(handles.US_Frame_Plot)
            



            plot(handles.US_Frame_Plot,(1:samples_Per_Frame),current_Frame);
            plot(handles.US_Frame_Plot,[wall1-points_Per_MM, wall1-points_Per_MM],[4,-4],'r');
            plot(handles.US_Frame_Plot,[wall1+points_Per_MM, wall1+points_Per_MM],[4,-4],'r');
            plot(handles.US_Frame_Plot,[wall2-points_Per_MM, wall2-points_Per_MM],[4,-4],'g');
            plot(handles.US_Frame_Plot,[wall2+points_Per_MM, wall2+points_Per_MM],[4,-4],'g');

            drawnow;
        end   
        if(stop_All ==1) 
            break; 
        end
        last_Track_Frame = current_Frame;
        frame_Num=frame_Num+1;
        
    end
    toc
    frame_Num
end
%hold(handles.US_Frame_Plot);
if (live ==0)
    fclose(US_Fid);
else
    niScope_Link(3,5195);
end
    


% --- Executes on button press in prev_Record.
function prev_Record_Callback(hObject, eventdata, handles)
% hObject    handle to prev_Record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global measurement_Count;
global dia_Store;
global dist_Store;
global time_Stamp_Store;
global data_Length_Store;
global D_d_Store;
global del_D_Store;
global current_Display_Record_Index;

if(current_Display_Record_Index>1)
current_Display_Record_Index = current_Display_Record_Index-1;
end

hold on
cla(handles.mmt_Store);
plot(handles.mmt_Store,dist_Store((1:data_Length_Store(current_Display_Record_Index) ),current_Display_Record_Index),'r');



% --- Executes on button press in next_Record.
function next_Record_Callback(hObject, eventdata, handles)
% hObject    handle to next_Record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global measurement_Count;
global dia_Store;
global dist_Store;
global time_Stamp_Store;
global data_Length_Store;
global D_d_Store;
global del_D_Store;
global current_Display_Record_Index;

if(current_Display_Record_Index<measurement_Count)
current_Display_Record_Index = current_Display_Record_Index+1;
end

hold on
cla(handles.mmt_Store);
plot(handles.mmt_Store,dist_Store((1:data_Length_Store(current_Display_Record_Index) ),current_Display_Record_Index),'r');


% --- Executes on button press in stop_But.
function stop_But_Callback(hObject, eventdata, handles)
% hObject    handle to stop_But (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stop_All;
stop_All=1;


% --- Executes on button press in TCF_But.
function TCF_But_Callback(hObject, eventdata, handles)
% hObject    handle to TCF_But (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global calibration_Data;
if(isempty(calibration_Data))
    [FileName,PathName] = uigetfile('*.tcf','Select the Compensation File');
    calibration_Data = dlmread(strcat(PathName,FileName));
    set(hObject,'String','Disable Transducer Compensation');
else
    calibration_Data=[];
    set(hObject,'String','Enable Transducer Compensation');
end

% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

global US_Sampling_Rate;
global samples_Per_Frame;
global transducer_Frequency;
global ECG_Sampling_Rate;

global live;
if (strcmp(get(eventdata.NewValue, 'Tag'),'rt'))
    live = 1
    US_Sampling_Rate = 100000000;
    samples_Per_Frame = 5195;
    transducer_Frequency = 500000;
    ECG_Sampling_Rate = 0;
else
    live = 0

end




%%%%%%%%%General Functions other than GUI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [state, wall1, wall2, time_Array, distension_Array, dia_Array, full_Cycle_Time_Array, full_Cycle_Dist_Array, full_Cycle_Lumen_Dia_Array, del_D, D_d] = track_Walls_no_ECG_no_TS(frame, Wall1, Wall2, pls_Rate)
global last_Track_Frame;
global current_Dist_Arr;
global PW_Wave;
global DW_Wave;
global pulse_Rate;
global current_State;
global samples_Per_Frame;
global points_Per_MM;
global bad_Signal_Quality_Count;
global dia_Wave;
global dia_Confidence;

distension_Array=0;
dia_Array=0;
full_Cycle_Time_Array = 0;
full_Cycle_Dist_Array =0;
full_Cycle_Lumen_Dia_Array = 0; 
del_D = 0; 
D_d = 0;
wall1=0;
wall2=0;
    
if(current_State ==2)
    pulse_Rate = pls_Rate;
    PW_Wave = [Wall1];
    DW_Wave = [Wall2];
    current_State=3;
    bad_Signal_Quality_Count=0;
    [dia_Wave,dia_Confidence] = find_Dia_L_Fit(frame,PW_Wave(end,1),DW_Wave(end,1),samples_Per_Frame,points_Per_MM);
end
PW_Wave = [PW_Wave; 0];
DW_Wave = [DW_Wave; 0];
dia_Wave = [dia_Wave; 0];
dia_Confidence = [dia_Confidence; 0];
[valid_Position, PW_Wave(end,1), DW_Wave(end,1)] = track_Walls(frame,last_Track_Frame, PW_Wave(end-1,1), DW_Wave(end-1,1),samples_Per_Frame,points_Per_MM);

if(valid_Position == 0)
    current_State = 5;
else
    if(check_Walls_Present(frame,PW_Wave(end),DW_Wave(end))== 0)
        bad_Signal_Quality_Count = bad_Signal_Quality_Count+1;
        if(bad_Signal_Quality_Count >= 10)
            bad_Signal_Quality_Count = 0;
            current_State = 5;
        end
    else 
        bad_Signal_Quality_Count =0;
        wall1= PW_Wave(end);
        wall2= DW_Wave(end);
        [dia_Wave(end,1),dia_Confidence(end,1)] = find_Dia_L_Fit(frame,PW_Wave(end,1),DW_Wave(end,1),samples_Per_Frame,points_Per_MM);
        [state_Change, on_Carotid, change_In_Dia, diastolic_Diameter, distension_Wave] = check_Cycles_ZC_Method(0);
        current_State = current_State+state_Change;
        full_Cycle_Dist_Array = distension_Wave;
        full_Cycle_Time_Array =  (1:length(distension_Wave))/pulse_Rate;
        del_D =  change_In_Dia;
        D_d = diastolic_Diameter;
    end
end
distension_Array = ((DW_Wave - PW_Wave)-mean((DW_Wave - PW_Wave)))/points_Per_MM;
dia_Array = dia_Wave;
time_Array = (1:length(dia_Wave))/pulse_Rate;
state =current_State;



function [valid_Positions, current_PW, current_DW] = track_Walls(current_Frame,previous_Frame, last_PW, last_DW, pts_Per_Frame, pts_Per_MM)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
valid_Positions = 1;
ROI1 = (last_PW - pts_Per_MM : last_PW + pts_Per_MM);  
slot1= previous_Frame(ROI1,1);
slot2= current_Frame(ROI1,1);
shift_Wall1= estimate_Shift(slot1,slot2,floor(pts_Per_MM/2));
current_PW = last_PW- shift_Wall1;

ROI2 = (last_DW - pts_Per_MM : last_DW + pts_Per_MM);
slot1= previous_Frame(ROI2,1);
slot2= current_Frame(ROI2,1);
shift_Wall2= estimate_Shift(slot1,slot2,floor(pts_Per_MM/2));
current_DW = last_DW- shift_Wall2;

if(current_PW < 2*pts_Per_MM || current_DW>pts_Per_Frame-2*pts_Per_MM) 
        valid_Positions =0; % The windows have went out of bounds
end



function [state, Wall1, Wall2] = find_Artery_Walls_SM(frame, reset)
global current_State;
global next_State;
global total_Consecutive_Findings
global Avg_Filt_mat;
global shift_Mat;
global tmp_Wall1;
global tmp_Wall2;
global tmp_Counter;
global samples_Per_Frame;

consecutive_Find_Limit = 3;
second_State_Multiplier = 5;

if(reset ==1)
    storage_Mat_Column_Nos = 5;
    Avg_Filt_mat = zeros(samples_Per_Frame,storage_Mat_Column_Nos);
    shift_Mat= (gen_Shift_Mat(storage_Mat_Column_Nos));
    total_Consecutive_Findings =0;
    tmp_Wall1 = zeros(second_State_Multiplier*consecutive_Find_Limit,1);
    tmp_Wall2 = zeros(second_State_Multiplier*consecutive_Find_Limit,1);
    current_State=0;
    next_State=0;
end

if(current_State == 0)
    [valid, Wall1, Wall2] = find_Artery_Walls(frame);
    total_Consecutive_Findings = total_Consecutive_Findings*valid + valid;
    if(valid ~=0)
        tmp_Wall1(total_Consecutive_Findings) = Wall1;
        tmp_Wall2(total_Consecutive_Findings) = Wall2;
    end
    if(total_Consecutive_Findings >=consecutive_Find_Limit)
        tmp_Counter=total_Consecutive_Findings;
        next_State =1;
    else
        next_State=0;
    end
elseif(current_State == 1)
    [valid,Wall1,Wall2] = find_Artery_Walls(frame);
    tmp_Wall1(total_Consecutive_Findings+1) = Wall1;
    tmp_Wall2(total_Consecutive_Findings+1) = Wall2;
    total_Consecutive_Findings = total_Consecutive_Findings+ valid;  %This is total finds and not consecutiv finds
    tmp_Counter = tmp_Counter +1;
    if(tmp_Counter >= second_State_Multiplier*consecutive_Find_Limit)
        Wall1 = round(median(tmp_Wall1(1:tmp_Counter)));
        Wall2 = round(median(tmp_Wall2(1:tmp_Counter)));
        if(Wall1~=0 && Wall2~=0)       
            next_State = 2;
            tmp_Counter=0;
        else
            total_Consecutive_Findings = 0;
            tmp_Counter=0;
            next_State = 0;
        end
    else
        next_State = 1;
    end
end
current_State = next_State;
state = current_State;


function [ valid, Wall1, Wall2] = find_Artery_Walls(frame)
valid =0;
Wall1=0;
Wall2=0;
global samples_Per_Frame;
global points_Per_MM;
global Avg_Filt_mat;
global shift_Mat;

Avg_Filt_mat = Avg_Filt_mat*shift_Mat;
Avg_Filt_mat(:,1) = frame;
corr_Frame = Corr_Piecewise_Fast(frame,frame);

V = peaks_3mm_Apart(corr_Frame,points_Per_MM);
I=0;
if(~isempty(V))
    I = find_Maximally_Uncorrelated(Avg_Filt_mat,V,samples_Per_Frame,points_Per_MM); 
    if(I>0)
        valid =1;
        Wall1=V(I,1);
        Wall2=V(I,2);
    end
end



function [ corr_Out ] = Corr_Piecewise_Fast(line1,line2)
line_X = zeros(size(line1));
corr_Out = zeros(size(line1));
line_X= line1.*line2;
fNorm = (500000) / (100000000/2);
[b,a] = butter(1, fNorm, 'low');
X1 = filtfilt(b, a, line_X);
corr_Out = X1(:,1);



function [ peak_Array] = peaks_3mm_Apart(inp_Signal,points_Per_MM )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
inp_Signal = inp_Signal/max(abs(inp_Signal));

T = 0.001;
P = [0 0];
L=0;
no_Of_Pairs =0;
peak_Array_Temp = zeros(10,2);
count = 0;
peak_Array=[];

while(T<0.5)
    [P,V] = peakdet(inp_Signal,T);
    L = size(P,1);
    if(L==0)
        peak_Array = [];
        return
    elseif(L >2 && L< 10)
        break;
    elseif(L < 3)
        T = T/5;
    else
        T = T*10;
    end
    count = count+1;
    if(count >10)
        break;
    end

for i = 1:L-1
    if(abs(P(i+1,1) - P(i,1)) > 3*points_Per_MM)
        Val_N = inp_Signal(P(i,1));
        Val_F = inp_Signal(P(i+1,1));
        Val_M = inp_Signal(ceil((P(i,1)+P(i+1,1))/2));
        
        if ((Val_N/Val_M)>10 && (Val_F/Val_M)>10)   % Check SNR for each wall is greater than 20 dB
            no_Of_Pairs = no_Of_Pairs+1;
            peak_Array_Temp(no_Of_Pairs,:) = [P(i,1) P(i+1,1)];
        end
    end
end
peak_Array = peak_Array_Temp([1:no_Of_Pairs],:);
end


function [maxtab, mintab] = peakdet(v, delta)
global maxtab_Mem;
global mintab_Mem;

S = size(v);
if(S(1)*S(2)==0)
    maxtab=[0,0];
    mintab =[0,0];
    return
end
max_Count =0;
min_Count =0;

mn = v(1,1); mx= v(1,1); mxpos=1;mnpos=1;
lookformax = 1;

for i=1:length(v)
  if lookformax
    if v(i) > mx, mx = v(i); mxpos =i;
    elseif v(i) < mx-delta
        if(mxpos ~=1)
          max_Count = max_Count+1;
          maxtab_Mem(max_Count,:) = [mxpos mx];
        end
      mn = v(i); mnpos = i;
      lookformax = 0;
    end  
  else
    if v(i) < mn, mn = v(i); mnpos =i;
    elseif v(i) > mn+delta
        if(mnpos~=S(1))
          min_Count = min_Count+1;
          mintab_Mem(min_Count,:) = [mnpos mn];
        end
      mx = v(i); mxpos = i;
      lookformax = 1;
    end
  end
end
maxtab = maxtab_Mem([1:max_Count],:);
mintab = mintab_Mem([1:min_Count],:);


function [ index ] = find_Maximally_Uncorrelated( block_Of_Frames, probe_Points, pts_Per_Frame, pts_Per_MM)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

L = size(probe_Points,1);
C = ones(L,1)*pts_Per_Frame;
index = 0;

for i = 1: L
[C(i,1),~,~] = distension_Probe(block_Of_Frames, probe_Points(i,1),probe_Points(i,2), pts_Per_Frame, pts_Per_MM);
end

[C1,I] = min(C(:,1));
if(C1<0)
    index = I;
end


function [ corr_Val, Wall1_Distension_Wave, Wall2_Distension_Wave] = distension_Probe(frames, Wall1_Loc, Wall2_Loc, pts_Per_Frame, pts_Per_MM)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

no_Of_Frames = size(frames,2);
corr_Val = 100.00;
Wall2_Distension_Wave = zeros(no_Of_Frames,1);
Wall1_Distension_Wave = zeros(no_Of_Frames,1);

ROI1 = (Wall1_Loc - pts_Per_MM : Wall1_Loc + pts_Per_MM);
ROI2 = (Wall2_Loc - pts_Per_MM: Wall2_Loc + pts_Per_MM);
Wall1_Distension_Wave(1,1) = Wall1_Loc;
Wall2_Distension_Wave(1,1) = Wall2_Loc;
for i = 1: no_Of_Frames-1
    if(ROI1(1) < 1.5*pts_Per_MM)
        return
    end
    slot1= frames(ROI1,i);
    slot2= frames(ROI1,i+1);
    shift_Wall1= estimate_Shift(slot1,slot2,floor(pts_Per_MM/2));
    ROI1= ROI1 + shift_Wall1;
    
    if(ROI2(length(ROI2)) >pts_Per_Frame)
        return
    end
    slot1= frames(ROI2,i);
    slot2= frames(ROI2,i+1);
    shift_Wall2= estimate_Shift(slot1,slot2,floor(pts_Per_MM/2));
    ROI2=ROI2+shift_Wall2;
    Wall1_Distension_Wave(i+1) = Wall1_Distension_Wave(i)+ shift_Wall1;%-(shift_Wall1+shift_Wall2)/2;
    Wall2_Distension_Wave(i+1) = Wall2_Distension_Wave(i)+ shift_Wall2;%-(shift_Wall1+shift_Wall2)/2;
end
W1_Mean_Sub = Wall1_Distension_Wave- mean(Wall1_Distension_Wave);
W2_Mean_Sub = Wall2_Distension_Wave- mean(Wall2_Distension_Wave);
S1_tmp= sum(W1_Mean_Sub.^2);
S2_tmp = sum(W2_Mean_Sub.^2);
S1 = S1_tmp(1,1);
S2 = S2_tmp(1,1);
if (S1>0 && S2 > 0)
    C = find_Corr(W1_Mean_Sub,W2_Mean_Sub,0,1);
else
    C=0;
end
corr_Val = C(1,1);


function [ shift ] = estimate_Shift(sig1,sig2,max_Shift)
% sig1(x+shift) = sig2(x)
%   Detailed explanation goes here
centre = max_Shift;
x = find_Corr(sig1(:,1),sig2(:,1),0,max_Shift);
[m,I] = max(x);
shift = I-centre;


function [ Corr_Sig ] = find_Corr( sig1,sig2,coeff, max_Ret_Len)
l1 = length(sig1);
l2 = length(sig2);
N = max(l1,l2);
if(max_Ret_Len >0)
    ret_Len = max_Ret_Len;
else
    ret_Len = N;
end

nfft = 2^nextpow2(2*N-1);
r_Cross = real(ifft( fft(sig1,nfft) .* conj(fft(sig2,nfft))));
Corr_Sig = [r_Cross(end-ret_Len+2:end) ; r_Cross(1:ret_Len)];

if(coeff ==1)
    r_Auto = real(ifft(fft(sig1,nfft) .* conj(fft(sig1,nfft))));
    Auto_Corr_Max = r_Auto(1);
    Corr_Sig = Corr_Sig/Auto_Corr_Max;
end


function [ frame ] = load_LV_Frame(fid,  no_Of_Pts_Per_Frame, fNo)
global calibration_Data;
frewind(fid);
fseek(fid,no_Of_Pts_Per_Frame*(fNo-1)*4,'bof');
[frame,count]=fread(fid,no_Of_Pts_Per_Frame,'*float','b');
frame= double(frame);

if(count <no_Of_Pts_Per_Frame) 
    disp (strcat('failed to read frame number ',int2str(fNo)));
    return;
end

if(length(calibration_Data)==length(frame)+3)
    hard_Limit_Point = calibration_Data(3);
    transducer_Response_Frame(:,1) = calibration_Data(4:end);
else
    hard_Limit_Point=0;
    transducer_Response_Frame(:,1) = zeros(size(frame));
end
    frame(:,1) = [zeros(hard_Limit_Point,1); frame(hard_Limit_Point+1:end)-transducer_Response_Frame((hard_Limit_Point+1:end))];
    
    
function [ shift_Mat ] = gen_Shift_Mat( dim )
shift_Mat = zeros(dim,dim);
shift_Mat = eye(dim);
shift_Mat = circshift(shift_Mat,[0 -1]);
shift_Mat(1,dim)=0;
shift_Mat = shift_Mat';


function [ status ] = check_Walls_Present( frame,Wall1, Wall2)
global points_Per_MM;
status =0;

if((Wall2-Wall1)/points_Per_MM < 4 || Wall1/points_Per_MM< 3 || (length(frame)-Wall2)/points_Per_MM < 3)
    status =0;
else
    Slot1= frame(Wall1-round(points_Per_MM/2):Wall1+round(points_Per_MM/2));
    Slot2= frame(Wall2-round(points_Per_MM/2):Wall2+round(points_Per_MM/2));
    PW_Slot_Avg = sum(abs(Slot1))/(2*round(points_Per_MM/2)+1);
    DW_Slot_Avg = sum(abs(Slot2))/(2*round(points_Per_MM/2)+1);
    Mid_Slot = frame(((round((Wall1+Wall2)/2))-points_Per_MM):((round((Wall1+Wall2)/2))+points_Per_MM));
    Mid_Slot_Avg = (sum(abs(Mid_Slot)))/(2*points_Per_MM+1);

    if(20*log(PW_Slot_Avg/Mid_Slot_Avg) > 20 && 20*log(DW_Slot_Avg/Mid_Slot_Avg) > 20)
        status =1;
    else
        status =0;
    end
end


function [ diameter_Lumen,confidence] = find_Dia_L_Fit(frame,Wall1,Wall2,No_Of_Samples,points_Per_mm)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global last_Dia;
mid_Pt = round((Wall1+Wall2)/2);
confidence=0;
diameter_Lumen =0;
if(mid_Pt < 2*points_Per_mm || mid_Pt > No_Of_Samples - 2*points_Per_mm)
    return
end

fNorm = (1000000) / (100000000/2);
[b,a] = butter(1, fNorm, 'low');
frame = filtfilt(b, a, frame.*frame);

left_Half = frame ([Wall1-round(1.5*points_Per_mm):mid_Pt],1);
left_Half = flipud(left_Half);
right_Half = frame ([mid_Pt:Wall2+round(1.5*points_Per_mm)], 1);

left_Half_Norm = left_Half/max(abs(left_Half));
right_Half_Norm = right_Half/max(abs(right_Half));

[~,P_R] = max(right_Half_Norm);
[~,P_L] = max(left_Half_Norm);

left_Half_Cut = left_Half_Norm([1: P_L(1,1)],1);
[left_Intima_Pt, left_L, sharpness_Index_Left] = Fit_L(left_Half_Cut);
right_Half_Cut = right_Half_Norm([1:P_R(1,1)],1);
[right_Intima_Pt, right_L, sharpness_Index_Right] = Fit_L(right_Half_Cut);
confidence = min(1/sharpness_Index_Left,1/sharpness_Index_Right);

diameter_Lumen = (left_Intima_Pt+ right_Intima_Pt)/points_Per_mm;



function [ inflection_X, L_Fitted, error] = Fit_L(well_Half)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


Tot_Len = length(well_Half);
sum_Sqr_Prev = Tot_Len;
L = 1*[1:Tot_Len];
L_Fitted = L;
error = 100;
inflection_X = Tot_Len;

if(Tot_Len <=1)
    return
end

well_Half = well_Half-min(well_Half);
well_Half = well_Half/(max(well_Half));
Jump = round(Tot_Len/5);
x= Tot_Len-1;
trial_Count = 0;
while(1)
    x = x- Jump;
    if(x<1  || x>Tot_Len-1)
        break;
    end
    L = L_Gen(well_Half(1),x, well_Half(x),well_Half(Tot_Len),Tot_Len);
    temp_Sqr = sum((L-well_Half).^2);
    if temp_Sqr > sum_Sqr_Prev
        Jump=-round(Jump/5);
        if(Jump == 0)
            inflection_X = x;
            L_Fitted = L;
            error= temp_Sqr;
            break;
        end
    end
    sum_Sqr_Prev = temp_Sqr;
    trial_Count = trial_Count+1;
    if(trial_Count>100)
        break;
    end
end



function [ L_Ret ] = L_Gen(yb,xm,ym,ye,total_Pts )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
xb=1;
xe= total_Pts;
line1 = linspace(yb,ym,xm-xb+1)';
line2 = linspace(ym,ye,xe-xm+1)';
L_Ret = [line1(1:length(line1)-1); line2];


function [ state_Change,  on_Carotid, Del_D, D_d, dist_Wave] = check_Cycles_ZC_Method(carotid_Check_Off0_On1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global PW_Wave ;
global DW_Wave;
global dia_Wave;
global points_Per_MM;
global pulse_Rate;
global dia_Confidence;

on_Carotid = 1;
state_Change=0;
Del_D=0;
D_d=0;
dist_Wave = (DW_Wave(:,1)-PW_Wave(:,1))/points_Per_MM;
dist_Wave = dist_Wave(:,1)-mean(dist_Wave(:,1));
dia_Wave = dia_Wave(:,1);

if (length(dist_Wave(:,1))>pulse_Rate/3)
    P_Wav = PW_Wave(:,1) -  (PW_Wave(:,1)+DW_Wave(:,1))/2;
    D_Wav = DW_Wave -  (PW_Wave(:,1)+DW_Wave(:,1))/2;
    K = find_Corr((P_Wav(:,1)-mean(P_Wav(:,1))),(D_Wav(:,1)-mean(D_Wav(:,1))),1,1);
    if(K(1,1) > -0.5)
        state_Change = 2;
        return;
    end
end

sub_Wave = zeros(1,1);
sub_Wave = zeros(2,1);
positive_Zero_Crossings = zeros(1,1);
positive_Zero_Crossings = zeros(2,1);
negative_Zero_Crossings = zeros(1,1);
negative_Zero_Crossings = zeros(2,1);

neg_Count= 0;
pos_Count = 0;

sign_Derivative_Distension = diff(sign(dist_Wave(:,1)));
last_Val=0;
for i = 1:length(sign_Derivative_Distension(:,1))
    if(sign_Derivative_Distension(i)~=0)
        if(last_Val == sign_Derivative_Distension(i))
            sign_Derivative_Distension(i) = 0;
        else
            if(sign_Derivative_Distension(i)<0)
                neg_Count = neg_Count+1;
                negative_Zero_Crossings(neg_Count,1) =  i;                    
            else
                pos_Count= pos_Count+1;
                positive_Zero_Crossings(pos_Count,1) =  i;
            end
        end
        last_Val = sign_Derivative_Distension(i);
    end
end

all_Crossings = [positive_Zero_Crossings(:,1); -1*negative_Zero_Crossings(:,1)];
[~, all_Crossings_Sort_Locs] = sort(abs(all_Crossings(:,1)));
all_Crossings = all_Crossings(all_Crossings_Sort_Locs(:,1));
all_Crossings_Copy = all_Crossings(:,1);
for i = 1:length(all_Crossings(:,1))-1
    if(abs(all_Crossings(i+1,1))-abs(all_Crossings(i,1))< pulse_Rate/10)
        all_Crossings_Copy(i,1) = 0;
        all_Crossings_Copy(i+1,1)=0;
    end
end

positive_Zero_Crossings = all_Crossings(find(all_Crossings_Copy(:,1)>0));
negative_Zero_Crossings = -1*all_Crossings(find(all_Crossings_Copy(:,1)<0));
pos_Count = length(positive_Zero_Crossings(:,1));
neg_Count = length(negative_Zero_Crossings(:,1));

if (pos_Count>1)
   sub_Wave = dist_Wave(positive_Zero_Crossings(1,1):positive_Zero_Crossings(2,1));
   dia_Sub_Wave = dia_Wave(positive_Zero_Crossings(1,1):positive_Zero_Crossings(2,1));
   confidence_Sub_Wave = dia_Confidence(positive_Zero_Crossings(1,1):positive_Zero_Crossings(2,1));
   [~,max_Loc] = max(sub_Wave(:,1));
   [~,min1_Loc] = min(sub_Wave(1:max_Loc,1));
   [~,min2_Loc] = min(sub_Wave(max_Loc:end,1));
   min2_Loc= min2_Loc + max_Loc-1;
   slope_Rise = abs((sub_Wave(max_Loc,1)-sub_Wave(min1_Loc,1))/(max_Loc-min1_Loc));
   slope_Fall = abs((sub_Wave(max_Loc,1)-sub_Wave(min2_Loc,1))/(min2_Loc-max_Loc));
   
   C_Dist_Dia = find_Corr((sub_Wave(:,1) -mean(sub_Wave(:,1))),(dia_Sub_Wave(:,1)- mean(dia_Sub_Wave(:,1))),1,1);
   if(~(carotid_Check_Off0_On1 && ~(slope_Rise > 3*slope_Fall)) && C_Dist_Dia(1,1) >0.5)
       [distension_Max_Val distension_Max_Point] = max(sub_Wave(:,1));
       [distension_Min_Val, distension_Min_Point] = min(sub_Wave(:,1));
       Del_D = (distension_Max_Val-distension_Min_Val);
       dia_Sub_Wave = corrected_Dia(dia_Sub_Wave,confidence_Sub_Wave, sub_Wave);
       D_d = dia_Sub_Wave(distension_Min_Point,1);
       PW_Wave = PW_Wave(positive_Zero_Crossings(1,1)+1:end);
       DW_Wave = DW_Wave(positive_Zero_Crossings(1,1)+1:end);
       dia_Wave = dia_Wave(positive_Zero_Crossings(1,1)+1:end);
       dia_Confidence = dia_Confidence(positive_Zero_Crossings(1,1)+1:end);
       dist_Wave = sub_Wave(:,1);
       state_Change=1;
   else
       sub_Wave=[0];
   end
elseif(neg_Count>1)
   sub_Wave = dist_Wave(negative_Zero_Crossings(1,1):negative_Zero_Crossings(2,1));
   dia_Sub_Wave = dia_Wave(negative_Zero_Crossings(1,1):negative_Zero_Crossings(2,1));
   confidence_Sub_Wave = dia_Confidence(negative_Zero_Crossings(1,1):negative_Zero_Crossings(2,1));
   [~,min_Loc] = min(sub_Wave(:,1));
   [~,max1_Loc] =  max(sub_Wave(1:min_Loc));
   [~,max2_Loc] =  max(sub_Wave(min_Loc:end));
   slope_Fall = abs((sub_Wave(max1_Loc,1)-sub_Wave(min_Loc,1))/(max1_Loc-min_Loc));
   slope_Rise = abs((sub_Wave(min_Loc,1)-sub_Wave(max2_Loc,1))/(max2_Loc-min_Loc));
   
   C_Dist_Dia = find_Corr((sub_Wave(:,1) -mean(sub_Wave(:,1))),(dia_Sub_Wave(:,1)- mean(dia_Sub_Wave(:,1))),1,1);
   if(~(carotid_Check_Off0_On1 && ~(slope_Rise > 3*slope_Fall)) && C_Dist_Dia(1,1) >0.5)
       [distension_Max_Val distension_Max_Point] = max(sub_Wave);
       [distension_Min_Val, distension_Min_Point] = min(sub_Wave);
       Del_D = (distension_Max_Val-distension_Min_Val);
       dia_Sub_Wave = corrected_Dia(dia_Sub_Wave,confidence_Sub_Wave, sub_Wave);
       D_d = dia_Sub_Wave(distension_Min_Point,1);
       PW_Wave = PW_Wave(negative_Zero_Crossings(1,1)+1:end);
       DW_Wave = DW_Wave(negative_Zero_Crossings(1,1)+1:end);
       dia_Wave = dia_Wave(negative_Zero_Crossings(1,1)+1:end);
       dia_Confidence = dia_Confidence(negative_Zero_Crossings(1,1)+1:end);
       dist_Wave = sub_Wave(:,1);
       state_Change=1;
   else
       sub_Wave=[0];
   end
end

if(length(PW_Wave) > 1.5*pulse_Rate)
    PW_Wave = PW_Wave(end-floor(1.5*pulse_Rate)+1:end);
    DW_Wave = DW_Wave(end-floor(1.5*pulse_Rate)+1:end);
    dia_Wave = dia_Wave(end-floor(1.5*pulse_Rate)+1:end);
    dia_Confidence = dia_Confidence(end-floor(1.5*pulse_Rate)+1:end);
end



function [ correct_Dia_Wav ] = corrected_Dia( dia_Wav, dia_Confidence, dist_Wav )
[~, best_Point] = max(dia_Confidence);

dist_Wav(:,1) = dist_Wav(:,1)-dist_Wav(best_Point);
correct_Dia_Wav = dist_Wav(:,1) + dia_Wav(best_Point);
