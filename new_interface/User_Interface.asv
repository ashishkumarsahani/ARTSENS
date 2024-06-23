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

% Last Modified by GUIDE v2.5 24-Jun-2014 22:52:21

% Begin initialization code - DO NOT EDIT
%                                 

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
%                                 


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
global ah;
global bh;
global al;
global bl;
global b_Corr;
global a_Corr;
global b_L_Fit;
global a_L_Fit;
global b_ecg_h;
global a_ecg_h;
global b_ecg_l;
global a_ecg_l;
global device;
global calibration_Data_loc
global calibration_Data;

fclose ('all');
fNorm = (1000000) / (100000000/2);
[bh,ah] = butter(1, fNorm, 'high');
fNorm = (10000000) / (100000000/2);
[bl,al] = butter(1, fNorm, 'low');
fNorm = (500000) / (100000000/2);
[b_Corr,a_Corr] = butter(1, fNorm, 'low');
fNorm = (1000000) / (100000000/2);
[b_L_Fit,a_L_Fit] = butter(1, fNorm, 'low');
fNorm = (0.5) / (500/2);
[b_ecg_h,a_ecg_h] = butter(1, fNorm, 'high');
fNorm = (30) / (500/2);
[b_ecg_l,a_ecg_l] = butter(4, fNorm, 'low');
global live;
live = 0;
global ecg_on;
ecg_on=0;
global running;
running =0;
set(handles.disp_wfm,'visible','on');
set(handles.carotid,'visible','on');
set(handles.femoral,'visible','on');
set(handles.read_File_Opt,'visible','on');
set(handles.real_Time_Opt,'visible','on');
set(handles.ECG_on,'enable','on');
set(handles.m_track,'visible','off');
set(handles.speed_Test_But,'visible','off');
global c_done;
global f_done;
c_done=0;
f_done=0;
global min_Wall_Dia;
global max_Wall_Dia;
global Wall_Lumen_SNR;

min_Wall_Dia = 4;
max_Wall_Dia = 10;
Wall_Lumen_SNR = 10;

global req_no_cycle;
req_no_cycle = 10;

set(handles.ecg_slider,'Min',1);
set(handles.ecg_slider,'Max',1.01);
set(handles.ecg_slider,'value',1);
global pulse_filt_on;
pulse_filt_on = -1;
if (pulse_filt_on ~= -1)
    set(handles.TCF_But,'String','on pulse filt');
end
set(handles.US_File_Addr,'string',strcat(pwd,'\','temp','.ust'));

global save_conti_RF;
save_conti_RF = 0; % make this value 1 to enable continous saving of rf waveform
if save_conti_RF ==1
    set(handles.Save_continous,'string','saving ON');
    set(handles.Save_continous,'BackgroundColor','g');
else
    set(handles.Save_continous,'string','save');
end

global bg_colour ;
bg_colour =[0.9412	0.9412	0.9412];
prog_bar_position=get (handles.prog_bar_base,'position');
prog_bar_max=prog_bar_position(3);
prog_bar_position(3)=prog_bar_max * 0.05;
set (handles.prog_bar,'position',prog_bar_position);
prog_bar2_position=get (handles.prog_bar2,'position');
prog_bar2_position(3)=prog_bar_max * 0.01;
set (handles.prog_bar2,'position',prog_bar2_position);
global c_HR;
global f_HR;
% c_HR=0;
% f_HR=0;
set(handles.cycle_Info,'String',{['D_d =  '],['Del_D = '], ...
                       ['Beta = '],['PAT = '],['cycle no =','0'],...
                       ['JV cycle no =','0']});
set (handles.cycle_no,'string','0');                   
                   
% accessing system default csv file
fid_sys= fopen('system_default.csv','r+');
if (fid_sys == -1 )
    % no system default file found so creat new system default file
    fid_sys= fopen('system_default.csv','w');
    device = inputdlg('scope device:', 'device input (eg. Dev1)',[1 50]);
    if(isempty(device))
        device = 'Dev1' ;
    end
    s_loc = pwd;
    set(handles.dir_loc,'string', s_loc);
    s_files = '0';
    auto_res = '0';
    ecg_on_str= '0';
    calibration_Data_loc ='';
    system_default = strcat(device,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',')
    fwrite(fid_sys,system_default{1});
    ecg_on= 0;
    
else

    system_default = fgetl(fid_sys);
    if (isempty(system_default) == 1 || length(system_default) < 12)
        fclose (fid_sys);
        fid_sys= fopen('system_default.csv','w');
        device = inputdlg('scope device:', 'device input (eg. Dev1)',[1 50]);
        if(isempty(device))
        device = 'Dev1' ;
        end
        s_loc = pwd;
        s_files = '0';
        auto_res = '0';
        ecg_on_str= '0';
        calibration_Data_loc ='';
        system_default = strcat(device,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',')
        fwrite(fid_sys,system_default{1});
        ecg_on= 0;
    else
        [device_name, s_loc, s_files, auto_res,  ecg_on_str,calibration_Data_loc,~] = strread(system_default, '%s %s %s %s %s %s %s', 'delimiter',','); %system_default = strcat(device,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',')
        device = device_name {1};
        set(handles.dir_loc,'string', s_loc{1});
        set(handles.save_file,'value', str2double (s_files{1}));
        set(handles.auto,'value', str2double (auto_res{1}));
        ecg_on= str2double(ecg_on_str{1});
        set(handles.ECG_on,'value', ecg_on);
        
        if(isempty(calibration_Data_loc)==0 || strcmp(calibration_Data_loc ,'')==0)
            calibration_Data_loc = calibration_Data_loc{1};
            try
                calibration_Data = dlmread(calibration_Data_loc);
            catch
                msgbox('Error in reading calibration file_1','error','error');
                calibration_Data_loc ='';
            end
            set(handles.TCF_But,'String','Disable Transducer Compensation');
        end
    end
end
fclose (fid_sys);            
set(gcf, 'units','normalized','OuterPosition',[0 0 1 1]);



% --- Outputs from this function are returned to the command line.
function varargout = User_Interface_OutputFcn(hObject, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in File_Sel_But.
function File_Sel_But_Callback(hObject, eventdata, handles)
% hObject    handle to File_Sel_But (see GCBO)
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
global PAT_Store;
global ECG_Store;
global ECG_TS_Store;
global current_Display_Record_Index;
global JV_mmt_Count;
global JV_current_Display_Record_Index;
global JV_Dist_Store;
global JV_Data_Length_Store;
global JV_Timestamp_Store;
global JV_ECG_Store;
global JV_ECG_TS_Store;
global consider;

global femoral;
global c_measurement_Count;
global c_consider;
global c_dia;
global c_dist;
global c_time_Stamp;
global c_data_Length;
global c_D_d;
global c_del_D;
global c_PAT;
global c_done;
global c_ECG ;
global c_ECG_TS; 

global c_JV_mmt_Count;
global c_JV_Dist;
global c_JV_Data_Length;
global c_JV_Timestamp;
global c_JV_ECG ;
global c_JV_ECG_TS; 
global c_HR;
c_HR= 0;
global f_measurement_Count;
global f_consider;
global f_dia;
global f_dist;
global f_time_Stamp;
global f_data_Length;
global f_D_d;
global f_del_D;
global f_PAT;
global f_done;
global f_ECG ;
global f_ECG_TS; 
global f_HR;
f_HR= 0;
global start_time;
global live;
global ecg_on;
global ecg_data;
global tstamp_ecg;
global ecg_disp_start;
global tstamp_Rpeak;
global carotid_wave_start;
global JV_wave_start;
global long_dist_store;
global dist_time_store;
global TS_on;
global running;
length_CN = 20;
length_NF = 60;

[FileName,PathName] = uigetfileCD(get(handles.US_File_Addr,'String'),'*.us;*.ust;*.lvsngl;*.mat','Select the ARTSENS Record');
current_File_Name = strcat(PathName,FileName);
set(handles.US_File_Addr,'String',current_File_Name); %
if(strcmp(FileName(end-5:end),'lvsngl')==1)
    US_Sampling_Rate = 100000000;
    samples_Per_Frame = 5195;
    transducer_Frequency = 500000;
    ECG_Sampling_Rate = 0;
    
    set(handles.read_File_Opt,'value',1);
    set(handles.real_Time_Opt,'value',0);
elseif (strcmp(FileName(end-2:end),'ust')==1)  
    US_Sampling_Rate = 100000000;
    samples_Per_Frame = 5195;
    transducer_Frequency = 500000;
    ECG_Sampling_Rate = 500;
    set(handles.read_File_Opt,'value',1);
    set(handles.real_Time_Opt,'value',0);
    [ecg_Fid,fail_msg] = fopen(strcat(current_File_Name(1:end-3),'ecg'));
    if(ecg_Fid== -1)
        
        ecg_on =0;
        set(handles.ECG_on,'value',0);
    else
        ecg_on =1;
        set(handles.ECG_on,'value',1);
        set(handles.ECG_on,'enable','on');
    end
    [mat_Fid,fail_msg] = fopen(strcat(current_File_Name(1:end-3),'mat'));
    if(mat_Fid ~= -1)
        current_File_Name=strcat(current_File_Name(1:end-3),'mat')
        load(current_File_Name,'SBP');
        set(handles.SBP_Edit,'String',num2str(SBP));
        load(current_File_Name,'DBP');
        set(handles.DBP_Edit,'String',num2str(DBP));
        load(current_File_Name,'length_CN');
        set(handles.length_CN_edit,'String',num2str(length_CN));
        load(current_File_Name,'length_NF');
        set(handles.length_NF_edit,'String',num2str(length_NF));
        load(current_File_Name,'name')
        set(handles.file_name,'String',name);
        
        
        try
        load(current_File_Name,'age')
        set(handles.age_yrs,'String',num2str(age));
        load(current_File_Name,'hight')
        set(handles.hight_cm,'String',num2str(hight));
        load(current_File_Name,'weight')
        set(handles.weight_kg,'String',num2str(weight));
        
        load(current_File_Name,'gender')       
        contents_gender = cellstr(get(handles.gender_popup,'String'));
        set(handles.gender_popup,'Value',strmatch(gender,contents_gender));
        load(current_File_Name,'posture')  
        contents_posture = cellstr(get(handles.posture_popup,'String'));
        set(handles.posture_popup,'Value',strmatch(posture,contents_posture));
        
        catch
            set(handles.age_yrs,'String','-');
            set(handles.hight_cm,'String','-');
            set(handles.weight_kg,'String','-');
            set(handles.gender_popup,'Value',1);
            set(handles.posture_popup,'Value',1);
        end
        
        set(handles.pwv_disp,'String','');
        set(handles.pwv_result,'String','');
        
    end
    
elseif (strcmp(FileName(end-2:end),'mat')==1)

    prev_handles = handles;
    load(current_File_Name);
    handles = prev_handles;  
    live = 0;
    running = 0 ;
    set(handles.read_File_Opt,'value',1);
    set(handles.real_Time_Opt,'value',0);
    cla(handles.dist_Axes,'reset');
    cla(handles.mmt_Store, 'reset');
    cla(handles.JV_mmt_Store, 'reset');
    hold(handles.dist_Axes, 'on');
    hold(handles.mmt_Store, 'on');
    hold(handles.JV_mmt_Store, 'on');
%     if (TS_on ==1 || live ==1)
        set(handles.dist_Axes,'ALimMode','manual');
        ylim(handles.dist_Axes,[0 1]);
        xlim(handles.dist_Axes,[0 1500]);
        cla(handles.dist_Axes);
        hold(handles.dist_Axes, 'on');
        cla(handles.mmt_Store, 'reset');
        xlim(handles.mmt_Store,[-100 1000]);
        hold(handles.mmt_Store, 'on');
        cla(handles.JV_mmt_Store, 'reset');
        xlim(handles.JV_mmt_Store,[-100 1000]) ;
        hold(handles.JV_mmt_Store, 'on');
%     end
    set(handles.US_Frame_Plot,'ALimMode','manual');
    ylim(handles.US_Frame_Plot,[-4 4]);
    xlim(handles.US_Frame_Plot,[1 samples_Per_Frame]);
    hold(handles.US_Frame_Plot, 'on');
    cla(handles.US_Frame_Plot);
    set(handles.ECG_Plot,'ALimMode','manual');
    ylim(handles.ECG_Plot,[0 2]);
    xlim(handles.ECG_Plot,[1 (14*500 * 2)]);
    hold(handles.ECG_Plot,'on');
    
    set(handles.ecg_slider,'Min',1);
    set(handles.ecg_slider,'Max',1.01);
    set(handles.ecg_slider,'value',1);
    if (measurement_Count >= 1)
       current_Display_Record_Index = 1;
       cla(handles.mmt_Store);
       plot(handles.mmt_Store,time_Stamp_Store(1:data_Length_Store(current_Display_Record_Index),current_Display_Record_Index), dist_Store((1:data_Length_Store(current_Display_Record_Index)),current_Display_Record_Index),'r');
       temp_length= length(find(ECG_TS_Store(:,current_Display_Record_Index)~=0));
       plot(handles.mmt_Store,ECG_TS_Store(1:temp_length,current_Display_Record_Index),ECG_Store(1:temp_length,current_Display_Record_Index),'b');
       [beta_Arr, ~, ~, best_beta, avg_Dia,avg_del_D,avg_PAT,SD_PAT, beta_Confidence] = calculate_Stiffness(del_D_Store(1:measurement_Count), D_d_Store(1:measurement_Count),PAT_Store(1:measurement_Count), SBP, DBP);
       set(handles.cycle_Info,'String',{['D_d =  ',num2str(D_d_Store(current_Display_Record_Index))],['Del_D = ', num2str(del_D_Store(current_Display_Record_Index))], ...
                       ['Beta = ', num2str(beta_Arr(current_Display_Record_Index))],['PAT = ', num2str(PAT_Store(current_Display_Record_Index))],['cycle no =',num2str(current_Display_Record_Index),' /',num2str(measurement_Count) ]...
                       ,['JV cycle no =',num2str(JV_current_Display_Record_Index),' /',num2str(JV_mmt_Count)]});
       set(handles.final_Results_Text,'String',{['Best Beta = ',num2str(best_beta)], ['Beta Confidence = ',num2str(beta_Confidence)],...
                       ['Average D_d = ',num2str(avg_Dia)],['Average Del_d = ',num2str(avg_del_D)], ['Average PAT = ',num2str(avg_PAT)],...
                       ['S.D. PAT = ',num2str(SD_PAT)]});
       set(handles.beta_disp,'String',num2str(best_beta));
       set(handles.cycle_no,'String',strcat(num2str(current_Display_Record_Index),' /',num2str(measurement_Count))); %
    end
    if(JV_mmt_Count >= 1)
        JV_current_Display_Record_Index=1;
       cla(handles.JV_mmt_Store);
       plot(handles.JV_mmt_Store,JV_Timestamp_Store(1:JV_Data_Length_Store(JV_current_Display_Record_Index),JV_current_Display_Record_Index), JV_Dist_Store((1:JV_Data_Length_Store(JV_current_Display_Record_Index)),JV_current_Display_Record_Index),'r');
       temp_length= length(find(JV_ECG_TS_Store(:,JV_current_Display_Record_Index)~=0));
       plot(handles.JV_mmt_Store,JV_ECG_TS_Store(1:temp_length, JV_current_Display_Record_Index), JV_ECG_Store(1:temp_length,JV_current_Display_Record_Index),'b');
       cycle_info_cell = get(handles.cycle_Info,'String');
       cycle_info_cell{end,1}= strcat('JV cycle no =',num2str(JV_current_Display_Record_Index),' /',num2str(JV_mmt_Count));
       set(handles.cycle_Info,'String',cycle_info_cell);
    end
     if (ecg_on ==1)  
        cla(handles.ECG_Plot);
        
        ECG_Data_Ch = (ecg_data(ecg_disp_start:end) - min(ecg_data(ecg_disp_start:end))) /(max(ecg_data(ecg_disp_start:end))- min(ecg_data(ecg_disp_start:end)));
        
        plot(handles.ECG_Plot,(tstamp_ecg(ecg_disp_start:end)-tstamp_ecg(ecg_disp_start)),ECG_Data_Ch);
%                 tempdata = ( long_dist_store -min(long_dist_store)) ;
        dist_disp_start = 1 + length (find(dist_time_store < tstamp_ecg(ecg_disp_start)));
        tempdata = ((long_dist_store(dist_disp_start:end))- min(long_dist_store(dist_disp_start:end)))/(max(abs(long_dist_store(dist_disp_start:end)))- min(long_dist_store(dist_disp_start:end))) +1;% ))-min(tempdata);

        plot(handles.ECG_Plot,(dist_time_store(dist_disp_start:end) - tstamp_ecg(ecg_disp_start)),tempdata,'r');

        for (i= 1:length(tstamp_Rpeak))
            if ((tstamp_Rpeak(i)-tstamp_ecg(ecg_disp_start) )> 0)
              plot(handles.ECG_Plot,[(tstamp_Rpeak(i)-tstamp_ecg(ecg_disp_start)),( tstamp_Rpeak(i)-tstamp_ecg(ecg_disp_start))],[2,0],'c');
            end
        end
        for (i= 1:length(carotid_wave_start))
            if ((carotid_wave_start(i)-tstamp_ecg(ecg_disp_start)) > 0)
              plot(handles.ECG_Plot,[(carotid_wave_start(i)-tstamp_ecg(ecg_disp_start)),( carotid_wave_start(i)-tstamp_ecg(ecg_disp_start))],[2,0],'k');
            end
        end
        for (i= 1:length(JV_wave_start))
            if ((JV_wave_start(i)-tstamp_ecg(ecg_disp_start)) > 0)
              plot(handles.ECG_Plot,[(JV_wave_start(i)-tstamp_ecg(ecg_disp_start)),( JV_wave_start(i)-tstamp_ecg(ecg_disp_start))],[2,0],'m');
            end
        end
        set(handles.ecg_slider,'Min',1);
        set(handles.ecg_slider,'Max',ecg_disp_start );
        set(handles.ecg_slider,'value',ecg_disp_start );

     end
        if femoral==1
            set(handles.disp_femoral,'value',1); 
           set(handles.disp_carotid,'value',0); 
        else
            set(handles.disp_femoral,'value',0); 
           set(handles.disp_carotid,'value',1); 
        end
        
    drawnow ;  
    set(handles.SBP_Edit,'String',num2str(SBP));
    set(handles.DBP_Edit,'String',num2str(DBP));
    set(handles.length_CN_edit,'String',num2str(length_CN));
    set(handles.length_NF_edit,'String',num2str(length_NF));
%     set(handles.file_name,'String',name);
    try
        set(handles.file_name,'String',name);
        set(handles.age_yrs,'String',num2str(age));
        set(handles.hight_cm,'String',num2str(hight));
        set(handles.weight_kg,'String',num2str(weight));
        contents_gender = cellstr(get(handles.gender_popup,'String'));
        set(handles.gender_popup,'Value',strmatch(gender,contents_gender));
        contents_posture = cellstr(get(handles.posture_popup,'String'));
        set(handles.posture_popup,'Value',strmatch(posture,contents_posture));
        set(handles.hr_disp,'string',num2str(Heart_Rate));
    catch
        set(handles.age_yrs,'String','-');
        set(handles.hight_cm,'String','-');
        set(handles.weight_kg,'String','-');
        set(handles.gender_popup,'Value',1);
        set(handles.posture_popup,'Value',1);
        set(handles.hr_disp,'string','');
    end
    
    set(handles.pwv_disp,'String','');
    set(handles.pwv_result,'String','');
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
global running;
running =1;
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
global PAT_Store;
global ECG_Store;
global ECG_TS_Store;
global current_Display_Record_Index;
global JV_mmt_Count;
global JV_current_Display_Record_Index;
global JV_Dist_Store;
global JV_Data_Length_Store;
global JV_Timestamp_Store;
global JV_ECG_Store;
global JV_ECG_TS_Store;
global consider;

global femoral;
global c_measurement_Count;
global c_consider;
global c_dia;
global c_dist;
global c_time_Stamp;
global c_data_Length;
global c_D_d;
global c_del_D;
global c_PAT;
global c_done;
global c_ECG ;
global c_ECG_TS; 

global c_JV_mmt_Count;
global c_JV_Dist;
global c_JV_Data_Length;
global c_JV_Timestamp;
global c_JV_ECG ;
global c_JV_ECG_TS; 
global c_HR;
        
global f_measurement_Count;
global f_consider;
global f_dia;
global f_dist;
global f_time_Stamp;
global f_data_Length;
global f_D_d;
global f_del_D;
global f_PAT;
global f_done;
global f_ECG ;
global f_ECG_TS; 
global f_HR;
      
global start_time;
global prev_Time_Stamp;
global ecg_data;
global tstamp_ecg;
global ecg_disp_start;
global tstamp_Rpeak;
global carotid_wave_start;
global JV_wave_start;
global long_dist_store;
global dist_time_store;

tstamp_Rpeak = zeros( 2, 1);
global prev_Rpeak_ind;
prev_Rpeak_ind=1;
global ah;
global bh;
global al;
global bl;
global b_ecg_l;
global a_ecg_l;
global device;
global calibration_Data;
global min_Wall_Dia;
global max_Wall_Dia;
global Wall_Lumen_SNR;
global TS_on;
global save_conti_RF;
global bg_colour;
cla(handles.mmt_Store);
cla(handles.JV_mmt_Store);
global live;
global ecg_on;

ecg_on =get(handles.ECG_on,'value');
ecg_lead_status =0;
ecg_data = zeros( 14*500, 1);
tstamp_ecg = ecg_data;
femoral = get(handles.femoral,'value');
global pulse_filt_on;
if (pulse_filt_on ~= -1)
load('pulse.mat','pulse');
pulse_Template = double(pulse(33:150));
pulse_Template = double(pulse_Template - mean(pulse_Template));
[~,shift_Dist] = max(abs(pulse_Template));
end
if (femoral==1)    
    min_Wall_Dia = 4;
    max_Wall_Dia = 12;
    Wall_Lumen_SNR = 5;
    set(handles.carotid,'visible','off');
    set(handles.femoral,'visible','on'); 
else
    set(handles.carotid,'visible','on');
    set(handles.femoral,'visible','off');
    min_Wall_Dia = 4;
    max_Wall_Dia = 10;
    Wall_Lumen_SNR = 10;
end
      
set(handles.disp_wfm,'visible','off');

live = get(handles.real_Time_Opt,'Value');
if (live==1)
    US_Sampling_Rate = 100000000;
    samples_Per_Frame = 5195;       %5195*3
    transducer_Frequency = 500000;
    ECG_Sampling_Rate = 500;
    set(handles.m_track,'visible','off');
    set(handles.read_File_Opt,'visible','off');
    set(handles.real_Time_Opt,'visible','on');
    if (ecg_on==1)
        ECG_Link(1,'test.ecg');
    end
else
    set(handles.m_track,'visible','on');
    set(handles.read_File_Opt,'visible','on');
    set(handles.real_Time_Opt,'visible','off');
end

points_Per_MM=130;
storage_Mat_Column_Nos = 15;
Avg_Filt_mat = zeros(samples_Per_Frame,storage_Mat_Column_Nos);
shift_Mat= (gen_Shift_Mat(storage_Mat_Column_Nos));

stop_All=0;

total_Frames=0;
US_Fid =0;
measurement_Count = 0;
dia_Store = zeros(100,1);
dist_Store = zeros(100,1);
JV_mmt_Count =0;
JV_current_Display_Record_Index = 0;
JV_Dist_Store = zeros(100,1);
JV_Timestamp_Store = zeros(100,1);
JV_Data_Length_Store=0;
D_d_Store = 0;
del_D_Store = 0;
PAT_Store=0;
ECG_Store =[];
ECG_TS_Store=[] ;
JV_ECG_Store=[];
JV_ECG_TS_Store=[];


frame_Rate = 0;
SBP = str2double(get(handles.SBP_Edit,'String'));
DBP = str2double(get(handles.DBP_Edit,'String'));
length_CN=str2double(get(handles.length_CN_edit,'String'));
length_NF=str2double(get(handles.length_NF_edit,'String'));

age=str2double(get(handles.age_yrs,'String'));
hight=str2double(get(handles.hight_cm,'String'));
weight=str2double(get(handles.weight_kg,'String'));

contents_gender = cellstr(get(handles.gender_popup,'String'));
gender=contents_gender{get(handles.gender_popup,'Value')};
contents_posture = cellstr(get(handles.posture_popup,'String'));
posture=contents_posture{get(handles.posture_popup,'Value')};

if  (isnan(SBP) ==1 || isnan(DBP) ==1)
    msgbox('please enter BP values ( SBP, DBP)','error','error');
    set(handles.disp_wfm,'visible','on');
    set(handles.carotid,'visible','on');
    set(handles.femoral,'visible','on');
    set(handles.read_File_Opt,'visible','on');
    set(handles.real_Time_Opt,'visible','on');
    return
end
set(handles.file_name,'Enable' ,'inactive');

if (live == 0)
    current_File_Name = get(handles.US_File_Addr,'String');% 
    if(strcmp(current_File_Name(end-5:end),'lvsngl')==1 || strcmp(current_File_Name(end-1:end),'us')==1 )
        [US_Fid,fail_msg] = fopen(current_File_Name);
        if(US_Fid== -1)
            disp(fail_msg)
            return
        end
        TS_on =0;
    elseif (strcmp(current_File_Name(end-2:end),'ust')==1)
            [US_Fid,fail_msg] = fopen(current_File_Name);
        if(US_Fid== -1)
            disp(fail_msg)
            return
        end
        TS_on =1;
        if (strcmp(current_File_Name(end-7:end-4),'_fem')==1) % file name - *****_fem.ust
            femoral= 1;
            min_Wall_Dia = 4;
            max_Wall_Dia = 12;
            Wall_Lumen_SNR = 5;
            set(handles.carotid,'visible','off','value',0);
            set(handles.femoral,'visible','on','value',1); 

        else
            femoral= 0;
            set(handles.carotid,'visible','on','value',1);
            set(handles.femoral,'visible','off','value',0);

            min_Wall_Dia = 4;
            max_Wall_Dia = 10;
            Wall_Lumen_SNR = 10;
        end
        
        [ecg_Fid,fail_msg] = fopen(strcat(current_File_Name(1:end-3),'ecg'));
        if(ecg_Fid== -1)
            disp(fail_msg)
            ecg_on =0;
            set(handles.ECG_on,'value',0);
            set(handles.ECG_on,'enable','inactive');
        else
            ecg_on = get(handles.ECG_on,'value');
%             set(handles.ECG_on,'value',1);
            set(handles.ECG_on,'enable','inactive');
            fseek(ecg_Fid, 0,'eof');
            ecg_length = (ftell(ecg_Fid))/4;
            frewind(ecg_Fid);
            ecg_info=fread(ecg_Fid,28,'*float','b');
            ECG_Sampling_Rate = ecg_info(3) ;
            ecg_info (27)= ecg_length;
            ecg_info (28)= 28;
        end
    elseif(strcmp(current_File_Name(end-2:end),'mat')==1)
%         load(current_File_Name)
        msgbox('invalid file name,  .mat not allowed');
        set(handles.disp_wfm,'visible','on');
        set(handles.carotid,'visible','on');
        set(handles.femoral,'visible','on');
        set(handles.read_File_Opt,'visible','on');
        set(handles.real_Time_Opt,'visible','on');
        set(handles.m_track,'visible','off');
        return

%         goto  display_data;
    else
        msgbox('invalid file name ...');
        set(handles.disp_wfm,'visible','on');
        set(handles.carotid,'visible','on');
        set(handles.femoral,'visible','on');
        set(handles.read_File_Opt,'visible','on');
        set(handles.real_Time_Opt,'visible','on');
        set(handles.m_track,'visible','off');
        return
    end
    fseek(US_Fid, 0,'eof');
    filelength = ftell(US_Fid);
    total_Frames = filelength/(4*samples_Per_Frame);
    prev_Time_Stamp = 0;
    TS =0;
    dist_time_store =[0];
    long_dist_store =[0];
    cla(handles.US_Frame_Plot,'reset');
    cla(handles.ECG_Plot);
    cla(handles.dist_Axes,'reset');
    cla(handles.mmt_Store, 'reset');
    cla(handles.JV_mmt_Store, 'reset');
    hold(handles.dist_Axes, 'on');
    hold(handles.mmt_Store, 'on');
    hold(handles.JV_mmt_Store, 'on');
    if (TS_on ==1)
        set(handles.dist_Axes,'ALimMode','manual');
        ylim(handles.dist_Axes,[0 1]);
        xlim(handles.dist_Axes,[0 1500]);
        cla(handles.dist_Axes);
        hold(handles.dist_Axes, 'on');
        cla(handles.mmt_Store, 'reset');
        xlim(handles.mmt_Store,[-100 1000]);
        hold(handles.mmt_Store, 'on');
        cla(handles.JV_mmt_Store, 'reset');
        xlim(handles.JV_mmt_Store,[-100 1000]) ;
        hold(handles.JV_mmt_Store, 'on');
    end
    
else         %(live == 1)                %for real time data signal
    
    total_Frames = 1000;
    save_file = get(handles.save_file,'value');
    if (save_file == 0)
        file_name = 'temp';
        if  (exist(strcat(file_name,'.ust'))==2)
            choice = questdlg({'Would you like to repeace temp file?',' No => give new file name & start again',' Auto rename => add date & time to temp file name'}, ...
                                'Replace saved file', ...
                                'Yes','No','Auto rename new file','Yes')
            if(strcmp(choice,'No')==1)
%                 file_name = inputdlg('New file Name:', 'File name input',[1 50])
%                 if(isempty(file_name))
                    return
%                 end
            elseif (strcmp(choice,'Auto rename new file')==1)
                formatOut = '_mmddyy_HH_MM';
                file_name = strcat(file_name,datestr(now,formatOut))
            end
        end
    else
        name = get(handles.file_name,'string');
        if (isempty(get(handles.dir_loc,'string')))
            set(handles.dir_loc,'string',pwd)
        end
        if (isempty(name))
            set(handles.file_name,'string','untitled')
            name = 'untitled';
        end
        [stat_dir,mess_dir,messid_dir] = mkdir(get(handles.dir_loc,'string'),name);
        if (stat_dir ==0)
            msgbox('folder can not created.. try again whith other folder name/loc');
            return
        elseif ((stat_dir ==1)&& (isempty(mess_dir)))
            [cleared_all]=clear_all_prev_data();
            disp('cleared_all.. new started');
        end
        file_name = strcat(get(handles.dir_loc,'string'),'\',name,'\',name,'_',get(handles.reading_no,'string'));
        formatOut = '_mmddyy_HH_MM';
        file_time =datestr(now,formatOut);
        file_name = strcat(file_name,file_time);
    end
    if (femoral==1)
        file_name = strcat(file_name,'_fem');
    else
        file_name = strcat(file_name,'_car');
    end
    if (save_file == 1)
        save(file_name);
    end
    [start_time] = niScope_Link(1,samples_Per_Frame,device,strcat(file_name,'.ust'));
    current_File_Name = strcat(file_name,'.ust');
    set(handles.US_File_Addr,'String',current_File_Name);
    TS_on =1;
    if (ecg_on ==1)
        niScope_Link(4,start_time);
        ECG_Link(5,(start_time +150),strcat(file_name,'.ecg'));
        set(handles.ECG_on,'enable','inactive');
%         save_ECG_Store = zeros( 28, 1);
%         save_ECG_Store(1) = 28;
%         save_ECG_Store(2) = 1;
%         save_ECG_Store(3) = ECG_Sampling_Rate;
%         save_ECG_Store(4) = 1;
%         save_ECG_Store(5) = 1;
%         save_ECG_Store(6) = start_time;
%         fs = fopen(strcat(file_name,'.ecg'), 'w+');
    end
    
    total_Elapsed_time=0;
    [current_Frame,TS] = niScope_Link(2,samples_Per_Frame,1);
    
    
    set(handles.dist_Axes,'ALimMode','manual');
    ylim(handles.dist_Axes,[0 1]);
    xlim(handles.dist_Axes,[0 1500]);
    cla(handles.dist_Axes);
    hold(handles.dist_Axes, 'on');
    cla(handles.mmt_Store, 'reset');
    xlim(handles.mmt_Store,[-100 1000]);
    hold(handles.mmt_Store, 'on');
    cla(handles.JV_mmt_Store, 'reset');
    xlim(handles.JV_mmt_Store,[-100 1000]) ;
    hold(handles.JV_mmt_Store, 'on');
    
    
end

set(handles.US_Frame_Plot,'ALimMode','manual');
ylim(handles.US_Frame_Plot,[-4 4]);
xlim(handles.US_Frame_Plot,[1 samples_Per_Frame]);
hold(handles.US_Frame_Plot, 'on');
set(handles.ECG_Plot,'ALimMode','manual');
ylim(handles.ECG_Plot,[0 2]);
xlim(handles.ECG_Plot,[1 (length(ecg_data) * 2)]);
hold(handles.ECG_Plot,'on');
set(handles.ecg_slider,'Min',0.9);
set(handles.ecg_slider,'Max',1.01);
set(handles.ecg_slider,'value',1);
if(strcmp(current_File_Name(end-5:end),'lvsngl')==1 || (strcmp(current_File_Name(end-2:end),'ust')==1) || live == 1)
    state1 = 0;
    state2 = 0;
    proc_time=0;
    ECG_Data_Ch = zeros( 1, 1);
    tstamp_ecg = ECG_Data_Ch;
    ecg_data = ECG_Data_Ch;
    frame_Num=1;
    ecg_sub_start =0;
    dist_time_store =[TS];
    long_dist_store =[0.5];
    JV_wave_start = [0];
    carotid_wave_start =[0];
    save_RF =save_conti_RF;
    set(handles.cycle_Info,'String',{['D_d =  '],['Del_D = '], ...
                       ['Beta = '],['PAT = '],['cycle no =','0'],...
                       ['JV cycle no =','0']});
    set(handles.final_Results_Text,'String', '');
    set(handles.pwv_disp,'String','');
    set(handles.beta_disp,'String','');
    set(handles.cycle_no,'String','0');
    set(handles.pwv_result,'String','');
    set(handles.hr_disp,'string','');
    dist_disp_start = 1;
    prev_manu_track= 0;% get(handles.track,'value'); 
    track_start=0;
    manu_track=0;
% profile on
    tic
    set(handles.state_ind,'BackgroundColor',bg_colour);
    set(handles.state_ind,'string','Finding Artery');
    
    prog_bar_position=get (handles.prog_bar_base,'position');
    prog_bar_max=prog_bar_position(3);
    prog_bar_position(3)=prog_bar_max * 0.01;
    set (handles.prog_bar,'position',prog_bar_position);
    set (handles.prog_bar,'string','0');
    prog_bar2_position=get (handles.prog_bar2,'position');
    prog_bar2_position(3)=prog_bar_max * 0.01;
    set (handles.prog_bar2,'position',prog_bar2_position);
    set (handles.prog_bar2,'string','');
    save_conti_RF =1;
    
    set(handles.Save_continous,'string','saving ON');
    set(handles.Save_continous,'BackgroundColor','g');
    while(frame_Num < total_Frames || (live== 1))
        frame_Rate = 1/toc;
        proc_time = proc_time + (1/frame_Rate);
        tic;
        
        
        
        if (live== 1)
            prev_Time_Stamp = single (TS);
            [Frame,TS] = niScope_Link(2,samples_Per_Frame,save_RF);
            TS = single(TS);
            Elapsed_time = TS - prev_Time_Stamp;
%             prev_Time_Stamp = TS;
            if((tstamp_ecg(end)-dist_time_store(dist_disp_start)) > 14000)
                i =dist_disp_start;
                while((tstamp_ecg(end)-dist_time_store(i))>14000)&& (i <length(dist_time_store))
                    i= i+1;
                end
                dist_disp_start = i;
%                 dist_time_store =[dist_time_store(i:end) ];
%                 long_dist_store =[long_dist_store(i:end) ];           
            end
            total_Elapsed_time = total_Elapsed_time + Elapsed_time  ;
            Frame(:,1) = Frame - mean(Frame);
            
            if((length(calibration_Data)==length(Frame)+3)&& (pulse_filt_on== -1) )
                hard_Limit_Point = calibration_Data(3);
                transducer_Response_Frame(:,1) = calibration_Data(4:end);
                Frame(:,1) = [zeros(hard_Limit_Point,1); Frame(hard_Limit_Point+1:end)-transducer_Response_Frame((hard_Limit_Point+1:end))];
               
            else
                hard_Limit_Point=0;
                transducer_Response_Frame(:,1) = zeros(size(Frame));
            end
                
               current_Frame = double(Frame);
%             current_Frame = double(current_Frame/max(abs(current_Frame)));%current_Frame = double(current_Frame);

%             current_Frame = filtfilt(bh, ah,current_Frame);
            
        elseif (live== 0)     % file read if live is not equal to 1
            prev_Time_Stamp = TS;
%             if track_start==1
%                 frame_Num = frame_Num-15;
%                   [~,prev_Time_Stamp]= load_LV_Frame(US_Fid, samples_Per_Frame,frame_Num,TS_on);
%                   frame_Num = frame_Num+1;
%             end
            [Frame,TS] = load_LV_Frame(US_Fid, samples_Per_Frame,frame_Num,TS_on);
            Elapsed_time = TS - prev_Time_Stamp;
            current_Frame = double (Frame);
            manu_track=get(handles.track,'value'); %prev_manu_track=(get(handles.track,'value');
            if (manu_track==1 && (state1 <2 || prev_manu_track==0))
                W1_tmp= str2double(get(handles.W1,'string'));
                W2_tmp= str2double(get(handles.W2,'string'));
                if ((W2_tmp > W1_tmp) && (1<W1_tmp && W1_tmp<5000)&& (1<W2_tmp && W2_tmp<5000) )
                    wall1= W1_tmp;
                    wall2= W2_tmp;
                    state1 = 2;
                    current_State=2;
                    prev_manu_track=1;
                end
            elseif (manu_track==0 && prev_manu_track==1)
                prev_manu_track=0;
%             elseif (manu_track==1 && prev_manu_track==1 && (state1 < 2 || state1 ==5 ))
%                 state1 = 3;
%                 %current_State=2;
            end
            prog_bar2_position(3)=frame_Num * prog_bar_max / total_Frames ;
            set (handles.prog_bar2,'position',prog_bar2_position);
            
    
        end
        if (pulse_filt_on==1)
            Corr_Frame = xcorr(current_Frame,pulse_Template);        
            current_Frame = Corr_Frame(samples_Per_Frame-shift_Dist:end-shift_Dist);
            current_Frame = double(4 * current_Frame/max(abs(current_Frame)));
        elseif (pulse_filt_on==0)
            current_Frame = double(4 * current_Frame/max(abs(current_Frame)));
        end
        if (ecg_on ==1)
            if (live==1) 
            [a,~, ch2, tstamp,ecg_lead_status] = ECG_Link(2,1);
%                 if (a ~=0)
%                     for i =1 :(a) ;% (i =0 ; i < (n/samples_per_frame) ; i++)
%                        save_ECG_Store = [save_ECG_Store ;ch2(i )] ;
%                        save_ECG_Store = [save_ECG_Store ;tstamp(i)] ;
%                     end
%                   end
                elseif (live==0)
                    [a,ch2,tstamp] = load_ecg_data(ecg_Fid, ecg_info,TS);
                    ecg_info(28)= ecg_info(28)+ (a*2);
                    if((tstamp_ecg(end)-dist_time_store(dist_disp_start)) > 14000)
                        i =dist_disp_start;
                        while((tstamp_ecg(end)-dist_time_store(i))>14000)&& (i <length(dist_time_store))
                            i= i+1;
                        end
                        dist_disp_start = i;
%                         dist_time_store =[dist_time_store(i:end) ];
%                         long_dist_store =[long_dist_store(i:end) ];
    
                    end
            end
            if  (a ~=0)
                if (length(tstamp_ecg)<7000)
                    ecg_disp_start = 1;
                    ecg_data = [ecg_data; ch2()];
                    tstamp_ecg = [tstamp_ecg; tstamp()];
                else
                    ecg_disp_start = length(ecg_data)-6999+a;
                    ecg_data = [ecg_data; ch2()];           % [ecg_data((end-6999+a):end); ch2()];
                    tstamp_ecg = [tstamp_ecg; tstamp()];    %[tstamp_ecg((end-6999+a):end); tstamp()];
                end
                ecg_sub_start = ecg_sub_start+a-1;
                if max (abs(ecg_data(ecg_disp_start:end)) ~=0)
                [ecg_sub_start] = ecg_processing(ecg_data(end-ecg_sub_start: end),tstamp_ecg(end-ecg_sub_start: end),ecg_sub_start);
                
                ECG_Data_Ch = (ecg_data(ecg_disp_start:end) - min(ecg_data(ecg_disp_start:end))) /(max(ecg_data(ecg_disp_start:end))- min(ecg_data(ecg_disp_start:end)));
                ECG_Data_Ch  = filtfilt(b_ecg_l, a_ecg_l, double( ECG_Data_Ch));
                else
                    ECG_Data_Ch = ecg_data(ecg_disp_start:end);
                end
            end
        end
        
        if(state1 <2)
            [state1,wall1,wall2] = find_Artery_Walls_SM (current_Frame,frame_Num);

            if ((live== 1)&&(mod(frame_Num,50)==0))
                    niScope_Link(4,start_time);
                    set(handles.state_ind,'BackgroundColor',bg_colour);
                    set(handles.state_ind,'string','Finding Artery');
            end
            save_RF = save_conti_RF;
            if state1 >=2
                track_start=1;
                save_RF = 1;
            end
%         
        else
            [state1, wall1, wall2, on_Carotid, time_Arr, dist_Arr, dia_Arr, full_Cycle_Time_Array, full_Cycle_Dist_Array, full_Cycle_Lumen_Dia_Array, del_D, D_d,full_Cycle_ECG,full_Cycle_ECG_Time,PAT] = track_Walls_wth_ECG(current_Frame, wall1, wall2, 50,TS);
%             long_dist_store(end) = dist_Arr(end);
            dist_time_store =[dist_time_store ;TS];
            long_dist_store =[long_dist_store ;dist_Arr(end)];
            save_RF =1;
            track_start=0;
            if (length(time_Arr)== length(dia_Arr) && time_Arr(1)~=0 && (mod(frame_Num,5)==0))
                if (ecg_on ==1 )
                    ecg_time_Arr = (tstamp_ecg(ecg_disp_start:end) - time_Arr(1));
                    time_Arr = (time_Arr - time_Arr(1));
                    
                    dist_Arr = (dist_Arr - min(dist_Arr)) /(max(dist_Arr)- min(dist_Arr));
                    cla(handles.dist_Axes);
                    plot(handles.dist_Axes,time_Arr, dist_Arr,'r');
                    plot(handles.dist_Axes,ecg_time_Arr, ECG_Data_Ch, 'b');
%                     drawnow;
                else %if (mod(frame_Num,5)==0)
                    dist_Arr = (dist_Arr - min(dist_Arr)) /(max(dist_Arr)- min(dist_Arr));
                    dia_Arr = (dia_Arr - min(dia_Arr)) /(max(dia_Arr)- min(dia_Arr));
                    time_Arr = (time_Arr - time_Arr(1));
                    cla(handles.dist_Axes);
                    plot(handles.dist_Axes,time_Arr, dist_Arr,'r');
                    plot(handles.dist_Axes,time_Arr, dia_Arr, 'g');
%                     drawnow;
                end
            elseif (mod(frame_Num,5)==0)
                time_Arr = (1:length(dia_Arr))/50;
                dist_Arr = (dist_Arr - min(dist_Arr)) /(max(dist_Arr)- min(dist_Arr));
                cla(handles.dist_Axes);
                plot(handles.dist_Axes,time_Arr, dist_Arr,'r');
                plot(handles.dist_Axes,time_Arr, dia_Arr-mean(dia_Arr), 'g');
            end
            if(state1==5)
                if (live== 1)
                    niScope_Link(4,start_time);
%                     save_RF =save_conti_RF;                
%                     if (ecg_on== 1)
%                         count = fwrite(fs,save_ECG_Store(1:end-29),'*float','b');
%                         save_ECG_Store= save_ECG_Store(end-28 :end);
%                     end
                end
               [state1, wall1, wall2] = find_Artery_Walls_SM(current_Frame, 1);
               set(handles.state_ind,'BackgroundColor',bg_colour);
               set(handles.state_ind,'string','Finding Artery');
            elseif(state1 ==4)
%                disp('one cycle extracted');
               if(on_Carotid==1)
                   measurement_Count = measurement_Count+1;
                   data_Length_Store(measurement_Count) = length(full_Cycle_Dist_Array);
                   consider(measurement_Count)=1;
                   time_Stamp_Store([1:length(full_Cycle_Time_Array)],measurement_Count) = full_Cycle_Time_Array;
                   dia_Store([1:length(full_Cycle_Lumen_Dia_Array)],measurement_Count) = full_Cycle_Lumen_Dia_Array;
                   dist_Store([1:length(full_Cycle_Dist_Array)],measurement_Count) = full_Cycle_Dist_Array;
                   dist_Store(:,measurement_Count)= dist_Store(:,measurement_Count);
                   D_d_Store(measurement_Count) = D_d;
                   del_D_Store(measurement_Count) = del_D;
                   PAT_Store(measurement_Count) = PAT;
                   current_Display_Record_Index = measurement_Count; 
                   carotid_wave_start(measurement_Count)=tstamp_Rpeak(prev_Rpeak_ind);
                   ECG_Store([1:length(full_Cycle_ECG)],measurement_Count) = full_Cycle_ECG;
                   ECG_TS_Store([1:length(full_Cycle_ECG)],measurement_Count)=full_Cycle_ECG_Time ;
                   

                   [beta_Arr, ~, ~, best_beta, avg_Dia,avg_del_D,avg_PAT,SD_PAT, beta_Confidence] = calculate_Stiffness(del_D_Store(1:measurement_Count), D_d_Store(1:measurement_Count),PAT_Store(1:measurement_Count), SBP, DBP);
%                    cla(H1,(handles.mmt_Store));
%                    cla(H2,(handles.mmt_Store));%[handles.AX h1 h2]=
                   cla(handles.mmt_Store);
                   plot(handles.mmt_Store,full_Cycle_Time_Array,full_Cycle_Dist_Array,'r');
                   plot(handles.mmt_Store,full_Cycle_ECG_Time,full_Cycle_ECG,'b');
                   set(handles.cycle_Info,'String',{['D_d =  ',num2str(D_d_Store(measurement_Count))],['Del_D = ', num2str(del_D_Store(measurement_Count))], ...
                       ['Beta = ', num2str(beta_Arr(measurement_Count))],['PAT = ', num2str(PAT_Store(measurement_Count))],...
                       ['cycle no =',num2str(measurement_Count)],['JV cycle no =',num2str(JV_mmt_Count)]});
                   set(handles.state_ind,'BackgroundColor','g');
                   set(handles.state_ind,'string','Artery found');
                   set(handles.final_Results_Text,'String',{['Best Beta = ',num2str(best_beta)], ['Beta Confidence = ',num2str(beta_Confidence)],...
                       ['Average D_d = ',num2str(avg_Dia)],['Average Del_d = ',num2str(avg_del_D)], ['Average PAT = ',num2str(avg_PAT)],...
                       ['S.D. PAT = ',num2str(SD_PAT)]});
                   set (handles.beta_disp,'string',num2str(best_beta));
                   set (handles.cycle_no,'string',num2str(measurement_Count));
                   
                   prog_bar_position(3)=beta_Confidence*prog_bar_max/100;
                   if (prog_bar_position(3) >=prog_bar_max) && (beta_Confidence ~=inf)
                       prog_bar_position(3) =prog_bar_max;
                       beta_Confidence = 100;
                       set (handles.prog_bar,'position',prog_bar_position);
                       set (handles.prog_bar,'string',num2str(beta_Confidence));
                   elseif ((prog_bar_position(3) <prog_bar_max)&& (beta_Confidence ~=inf))
                       set (handles.prog_bar,'position',prog_bar_position);
                       set (handles.prog_bar,'string',num2str(beta_Confidence));
                   end
                   
               elseif(on_Carotid==-1)
                   JV_mmt_Count = JV_mmt_Count+1;
                   JV_Data_Length_Store(JV_mmt_Count)=  length(full_Cycle_Dist_Array);
                   JV_current_Display_Record_Index = JV_mmt_Count;
                   JV_Dist_Store([1:JV_Data_Length_Store(JV_mmt_Count)],JV_mmt_Count) = full_Cycle_Dist_Array;
                   JV_Timestamp_Store([1:JV_Data_Length_Store(JV_mmt_Count)],JV_mmt_Count) = full_Cycle_Time_Array;
                   JV_ECG_Store([1:length(full_Cycle_ECG)],JV_mmt_Count) = full_Cycle_ECG;
                   JV_ECG_TS_Store([1:length(full_Cycle_ECG)],JV_mmt_Count)=full_Cycle_ECG_Time ;
                   JV_wave_start(JV_mmt_Count)=tstamp_Rpeak(prev_Rpeak_ind);
%                    cla(BX(1),(handles.JV_mmt_Store));
%                    cla(BX(2),(handles.JV_mmt_Store));
%                   
%                    [BX hb1 hb2]=plotyy(handles.JV_mmt_Store,full_Cycle_Time_Array, full_Cycle_Dist_Array,full_Cycle_ECG_Time,full_Cycle_ECG);
                   cla(handles.JV_mmt_Store);
                   plot(handles.JV_mmt_Store,JV_Timestamp_Store(1:JV_Data_Length_Store(JV_mmt_Count),JV_mmt_Count), JV_Dist_Store((1:JV_Data_Length_Store(JV_mmt_Count)),JV_mmt_Count),'r');
                   plot(handles.JV_mmt_Store,JV_ECG_TS_Store(1:length(full_Cycle_ECG), JV_mmt_Count), JV_ECG_Store(1:length(full_Cycle_ECG),JV_mmt_Count),'b');
                   cycle_info_cell = get(handles.cycle_Info,'String');
                   cycle_info_cell{end,1}= strcat('JV cycle no =',num2str(JV_mmt_Count));
                   set(handles.cycle_Info,'String',cycle_info_cell);
                   set(handles.state_ind,'BackgroundColor','m');
                   set(handles.state_ind,'string','Vein found');
               elseif(on_Carotid==0) 
                   set(handles.state_ind,'BackgroundColor',bg_colour);
                   set(handles.state_ind,'string','Not Artery!');
               end
               current_State =3;
            end,
        end
        disp_no=5;
        if(mod(frame_Num,disp_no)==0)
            set(handles.speed_Indicator_Text,'String',{['Frame Rate : ', num2str(frame_Rate)],['No : ',num2str(frame_Num),'/',num2str(total_Frames),' lead:',num2str(ecg_lead_status)]})
            cla(handles.US_Frame_Plot);
            
            current_Frame = filtfilt(bh, ah, current_Frame);
            plot(handles.US_Frame_Plot,(1:samples_Per_Frame),current_Frame,'b');
            plot(handles.US_Frame_Plot,[wall1-points_Per_MM, wall1-points_Per_MM],[4,-4],'r');
            plot(handles.US_Frame_Plot,[wall1+points_Per_MM, wall1+points_Per_MM],[4,-4],'r');
            plot(handles.US_Frame_Plot,[wall2-points_Per_MM, wall2-points_Per_MM],[4,-4],'g');
            plot(handles.US_Frame_Plot,[wall2+points_Per_MM, wall2+points_Per_MM],[4,-4],'g');
            
            
            if (ecg_on ==1 && (mod(frame_Num,disp_no*16)==0  ))   %|| live==0
                if (manu_track==0 && state1>=2 && wall1 ~=0 && get(handles.disp_tracking,'value')==1)
                set(handles.W1,'string',num2str(wall1));
                set(handles.W2,'string',num2str(wall2));%wall1, wall2 manu_track==1
                end
                cla(handles.ECG_Plot);
                set(handles.ecg_slider,'Max',ecg_disp_start );
                set(handles.ecg_slider,'value',ecg_disp_start );
                plot(handles.ECG_Plot,(tstamp_ecg(ecg_disp_start:end)-tstamp_ecg(ecg_disp_start)),ECG_Data_Ch);
%                 tempdata = ( long_dist_store -min(long_dist_store)) ;
                tempdata = ((long_dist_store(dist_disp_start:end))- min(long_dist_store(dist_disp_start:end)))/(max(abs(long_dist_store(dist_disp_start:end)))- min(long_dist_store(dist_disp_start:end))) +1;% ))-min(tempdata);

                plot(handles.ECG_Plot,(dist_time_store(dist_disp_start:end) - tstamp_ecg(ecg_disp_start)),tempdata,'r');
           
                for (i= 1:length(tstamp_Rpeak))
                    if ((tstamp_Rpeak(i)-tstamp_ecg(ecg_disp_start) )> 0)
                      plot(handles.ECG_Plot,[(tstamp_Rpeak(i)-tstamp_ecg(ecg_disp_start)),( tstamp_Rpeak(i)-tstamp_ecg(ecg_disp_start))],[2,0],'c');
                    end
                end
                if length(tstamp_Rpeak) >3
                Rpeak_time_diff=diff(tstamp_Rpeak);
                best_HR= round((1000/median(Rpeak_time_diff))*60);
                set(handles.hr_disp,'string',num2str(best_HR));
                end
                for (i= 1:length(carotid_wave_start))
                    if ((carotid_wave_start(i)-tstamp_ecg(ecg_disp_start)) > 0)
                      plot(handles.ECG_Plot,[(carotid_wave_start(i)-tstamp_ecg(ecg_disp_start)),( carotid_wave_start(i)-tstamp_ecg(ecg_disp_start))],[2,0],'k');
                    end
                end
                for (i= 1:length(JV_wave_start))
                    if ((JV_wave_start(i)-tstamp_ecg(ecg_disp_start)) > 0)
                      plot(handles.ECG_Plot,[(JV_wave_start(i)-tstamp_ecg(ecg_disp_start)),( JV_wave_start(i)-tstamp_ecg(ecg_disp_start))],[2,0],'m');
                    end
                end
                
            end
        
        
            drawnow ;   
        end
        if(stop_All ==1) 
            break; 
        end
        last_Track_Frame = current_Frame;
        frame_Num=frame_Num+1;
    end
    toc;
%     profile viewer
    frame_Num = frame_Num
    if (ecg_on ==1)
    Rpeak_time_diff=diff(tstamp_Rpeak);
    Rpeak_time_diff=Rpeak_time_diff(find(Rpeak_time_diff > 400));
    Rpeak_time_diff=Rpeak_time_diff(find(Rpeak_time_diff <= 2100));
    Heart_Rate= round((1000/median(Rpeak_time_diff))*60);
    set(handles.hr_disp,'string',num2str(Heart_Rate));
    else
        Heart_Rate = '--';
    end
    if (femoral== 0)
       c_measurement_Count=measurement_Count ;
       c_consider = ones(measurement_Count,1);
       c_data_Length = data_Length_Store;
       c_time_Stamp= time_Stamp_Store;
       c_dia = dia_Store;
       c_dist=dist_Store;
       c_ECG =  ECG_Store;
       c_ECG_TS =  ECG_TS_Store;
       
       c_D_d = D_d_Store;
       c_del_D = del_D_Store;
       c_PAT = PAT_Store;
       c_done=1; 
       
       c_JV_mmt_Count=JV_mmt_Count;
       c_JV_Data_Length=  JV_Data_Length_Store;
       c_JV_Dist = JV_Dist_Store;
       c_JV_Timestamp = JV_Timestamp_Store;
       c_JV_ECG =  JV_ECG_Store;
       c_JV_ECG_TS =  JV_ECG_TS_Store;
       c_HR= Heart_Rate;
       set(handles.disp_carotid,'value',1);
       set(handles.disp_femoral,'value',0); 
    elseif(femoral== 1)
       f_measurement_Count=measurement_Count ;
       f_consider = ones(measurement_Count,1);
       f_data_Length = data_Length_Store;
       f_time_Stamp= time_Stamp_Store;
       f_dia = dia_Store;
       f_dist=dist_Store;
       f_ECG =  ECG_Store;
       f_ECG_TS =  ECG_TS_Store;
       
       f_D_d = D_d_Store;
       f_del_D = del_D_Store;
       f_PAT = PAT_Store;
       f_done=1;
       f_HR= Heart_Rate;
       set(handles.disp_femoral,'value',1); 
       set(handles.disp_carotid,'value',0); 
    end
    set(handles.disp_wfm,'visible','on');
    set(handles.carotid,'visible','on');
    set(handles.femoral,'visible','on');
    set(handles.read_File_Opt,'visible','on');
    set(handles.real_Time_Opt,'visible','on');
    SBP = str2double(get(handles.SBP_Edit,'String'));
   DBP = str2double(get(handles.DBP_Edit,'String'));
   length_CN=str2double(get(handles.length_CN_edit,'String'));
   length_NF=str2double(get(handles.length_NF_edit,'String'));
end
if (live ==0)
    fclose(US_Fid);
    if ecg_on==1
        fclose(ecg_Fid);
    end
%     set(handles.ECG_on,'value',0);
    set(handles.ECG_on,'enable','on');
    save_file = get(handles.save_file,'value');
    if save_file ==1
        save(strcat(current_File_Name(1:end-4),'_simu'));
        set(handles.US_File_Addr,'String',strcat(current_File_Name(1:end-4),'_simu.mat'));
    elseif save_file ==0
        choice = questdlg({'Would you like to save result?'}, ...
                                'save result file', ...
                                'Yes','No','Yes');
        if(strcmp(choice,'Yes')==1) 
            save(strcat(current_File_Name(1:end-4),'_simu'));
            set(handles.US_File_Addr,'String',strcat(current_File_Name(1:end-4),'_simu.mat'));
        end
    end
    set(handles.m_track,'visible','off');
    
else %(live ==1)
   
    total_Elapsed_time = total_Elapsed_time
   
    disp((total_Elapsed_time/1000)/frame_Num)
     proc_time =proc_time
    disp((proc_time)/frame_Num)
    
    
   if ecg_on==1
       ECG_Link (6,1);
       ECG_Link(3,0);
       set(handles.ECG_on,'enable','on');
%        set(handles.ECG_on,'value',1);
%     count = fwrite(fs,save_ECG_Store,'*float','b');
%     fclose(fs);
%      save_ECG_Store =[0]
   end
   niScope_Link(3,samples_Per_Frame,1);
   
   
   if save_file ==1
%        save(current_File_Name(1:end-4));
       current_File_Name = strcat (get(handles.dir_loc,'string'),'\',name,'\',name,'_',get(handles.reading_no,'string'),'_',file_name(end-15:end));%
       save(current_File_Name);
%        set(handles.US_File_Addr,'String',strcat(current_File_Name,'.mat'));
   elseif save_file ==0
        choice = questdlg({'Would you like to save result?'}, ...
                                'save result file', ...
                                'Yes','No','Yes');
        if(strcmp(choice,'Yes')==1) 
            save(current_File_Name(1:end-4));
            
            set(handles.US_File_Addr,'String',strcat(current_File_Name(1:end-4),'.mat'));
        end
   end
end
running =0;




% --- Executes when selected object is changed in disp_wfm.
function disp_wfm_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in disp_wfm 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global disp_femoral;
global measurement_Count;
global consider;
global dia_Store;
global dist_Store;
global time_Stamp_Store;
global data_Length_Store;
global D_d_Store;
global del_D_Store;
global PAT_Store;
global current_Display_Record_Index;
global JV_mmt_Count;
global JV_current_Display_Record_Index;
global JV_Dist_Store;
global JV_Data_Length_Store;
global JV_Timestamp_Store;
global ECG_TS_Store;
global ECG_Store;
global JV_ECG_Store ;
global JV_ECG_TS_Store; 



global c_measurement_Count;
global c_consider;
global c_dia;
global c_dist;
global c_time_Stamp;
global c_data_Length;
global c_D_d;
global c_del_D;
global c_PAT;
global c_done;
global c_ECG ;
global c_ECG_TS ;

global c_JV_mmt_Count;
global c_JV_Dist;
global c_JV_Data_Length;
global c_JV_Timestamp;
global c_JV_ECG  ;
global c_JV_ECG_TS ;
global c_HR;
        
global f_measurement_Count;
global f_consider;
global f_dia;
global f_dist;
global f_time_Stamp;
global f_data_Length;
global f_D_d;
global f_del_D;
global f_PAT;
global f_done;
global f_ECG ;
global f_ECG_TS ;
global f_HR;

disp_femoral = get(handles.disp_femoral,'value');
disp_carotid = get(handles.disp_carotid,'value');
if (disp_femoral== 0)
    %c_done=1; 
    measurement_Count=c_measurement_Count ;
    consider = c_consider;
    data_Length_Store=c_data_Length ;
    time_Stamp_Store=c_time_Stamp;
    dia_Store=c_dia ;
    dist_Store=c_dist;
    ECG_Store =  c_ECG;
    ECG_TS_Store =c_ECG_TS ;

    D_d_Store=c_D_d;
    del_D_Store=c_del_D ;
    PAT_Store= c_PAT ;
    
    JV_mmt_Count = c_JV_mmt_Count;
    JV_Data_Length_Store= c_JV_Data_Length ;
    JV_Dist_Store = c_JV_Dist  ;
    JV_Timestamp_Store= c_JV_Timestamp  ;
    JV_ECG_Store = c_JV_ECG  ;
    JV_ECG_TS_Store=   c_JV_ECG_TS ;
    Heart_Rate=c_HR;

elseif(disp_femoral== 1)
    %f_done=1;
    measurement_Count=f_measurement_Count ;
    consider = f_consider;
    data_Length_Store=f_data_Length  ;
    time_Stamp_Store= f_time_Stamp ;
    dia_Store=f_dia  ;
    dist_Store=f_dist;
    ECG_Store =  f_ECG;
    ECG_TS_Store =f_ECG_TS ;

    D_d_Store=f_D_d  ;
    del_D_Store=f_del_D  ;
    PAT_Store=f_PAT  ;
    JV_mmt_Count =0;
    Heart_Rate=f_HR;
end
set(handles.hr_disp,'string',num2str(Heart_Rate));
%%%% displaying 1st waveform with total avg vlaue%%%
if(measurement_Count>0)
    current_Display_Record_Index =1;
    axes(handles.mmt_Store);
    hold on
    cla(handles.mmt_Store);
    plot(handles.mmt_Store,time_Stamp_Store(1:data_Length_Store(current_Display_Record_Index),current_Display_Record_Index), dist_Store((1:data_Length_Store(current_Display_Record_Index)),current_Display_Record_Index),'r');
    temp_length= length(find(ECG_TS_Store(:,current_Display_Record_Index)~=0));
    plot(handles.mmt_Store,ECG_TS_Store(1:temp_length,current_Display_Record_Index),ECG_Store(1:temp_length,current_Display_Record_Index),'b');
    SBP = str2double(get(handles.SBP_Edit,'String'));
    DBP = str2double(get(handles.DBP_Edit,'String'));
    [beta_Arr, ~, ~, best_beta, avg_Dia,avg_del_D,avg_PAT,SD_PAT, beta_Confidence] = calculate_Stiffness(del_D_Store(1:measurement_Count), D_d_Store(1:measurement_Count),PAT_Store(1:measurement_Count), SBP, DBP);
    
    set(handles.final_Results_Text,'String',{['Best Beta = ',num2str(best_beta)], ['Beta Confidence = ',num2str(beta_Confidence)], ['Average D_d = ',num2str(avg_Dia)],...
        ['Average Del_d = ',num2str(avg_del_D)], ['Average PAT = ',num2str(avg_PAT)], ['S.D. PAT = ',num2str(SD_PAT)]}); 
    set (handles.beta_disp,'string',num2str(best_beta));
    median_del_D=median(del_D_Store(1:measurement_Count));
    set(handles.consider,'value',consider(current_Display_Record_Index));
    set(handles.cycle_Info,'String',{['D_d =  ',num2str(D_d_Store(current_Display_Record_Index))],['Del_D = ', num2str(del_D_Store(current_Display_Record_Index))], ...
                       ['Beta = ', num2str(beta_Arr(current_Display_Record_Index))],['PAT = ', num2str(PAT_Store(current_Display_Record_Index))],['cycle no =',num2str(current_Display_Record_Index),' /',num2str(measurement_Count) ]...
                       ,['JV cycle no =',num2str(JV_current_Display_Record_Index),' /',num2str(JV_mmt_Count)]});
    set (handles.cycle_no,'string',strcat(num2str(current_Display_Record_Index),' /',num2str(measurement_Count)));
else
    cla(handles.mmt_Store);
end

if(JV_mmt_Count>0)
    
    JV_current_Display_Record_Index=1;
    
    axes(handles.JV_mmt_Store);
    hold on
    cla(handles.JV_mmt_Store);
    plot(handles.JV_mmt_Store,JV_Timestamp_Store(1:JV_Data_Length_Store(JV_current_Display_Record_Index),JV_current_Display_Record_Index), JV_Dist_Store((1:JV_Data_Length_Store(JV_current_Display_Record_Index)),JV_current_Display_Record_Index),'r');
    temp_length= length(find(JV_ECG_TS_Store(:,JV_current_Display_Record_Index)~=0));
    plot(handles.JV_mmt_Store,JV_ECG_TS_Store(1:temp_length,current_Display_Record_Index),JV_ECG_Store(1:temp_length,current_Display_Record_Index),'b');
    cycle_info_cell = get(handles.cycle_Info,'String');
    cycle_info_cell{end,1}= strcat('JV cycle no =',num2str(JV_mmt_Count));
    set(handles.cycle_Info,'String',cycle_info_cell);
else
    cla(handles.JV_mmt_Store);

end






% --- Executes on button press in prev_Record.
function prev_Record_Callback(hObject, eventdata, handles)
% hObject    handle to prev_Record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global measurement_Count;
global consider;
global dia_Store;
global dist_Store;
global PAT_Store;
global time_Stamp_Store;
global data_Length_Store;
global D_d_Store;
global del_D_Store;
global current_Display_Record_Index;
global JV_mmt_Count;
global JV_current_Display_Record_Index;
global JV_Dist_Store;
global JV_Data_Length_Store;
global JV_Timestamp_Store;
global ECG_TS_Store;
global ECG_Store;
global JV_ECG_Store ;
global JV_ECG_TS_Store; 

if(measurement_Count>0)
    if(current_Display_Record_Index>1)
        current_Display_Record_Index = current_Display_Record_Index-1;
    end
    axes(handles.mmt_Store);
    hold on
    cla(handles.mmt_Store);
    plot(handles.mmt_Store,time_Stamp_Store(1:data_Length_Store(current_Display_Record_Index),current_Display_Record_Index), dist_Store((1:data_Length_Store(current_Display_Record_Index)),current_Display_Record_Index),'r');
    temp_length= length(find(ECG_TS_Store(:,current_Display_Record_Index)~=0));
    plot(handles.mmt_Store,ECG_TS_Store(1:temp_length,current_Display_Record_Index),ECG_Store(1:temp_length,current_Display_Record_Index),'b');
    SBP = str2double(get(handles.SBP_Edit,'String'));
    DBP = str2double(get(handles.DBP_Edit,'String'));
    [beta_Arr, compliance_Arr, distensibility_Arr, best_beta, avg_Dia,avg_del_D,avg_PAT,SD_PAT, beta_Confidence] = calculate_Stiffness(del_D_Store(1:measurement_Count), D_d_Store(1:measurement_Count),PAT_Store(1:measurement_Count), SBP, DBP);
    set(handles.cycle_Info,'String',{['D_d =  ',num2str(D_d_Store(current_Display_Record_Index))],['Del_D = ', num2str(del_D_Store(current_Display_Record_Index))], ...
                       ['Beta = ', num2str(beta_Arr(current_Display_Record_Index))],['PAT = ', num2str(PAT_Store(current_Display_Record_Index))],['cycle no =',num2str(current_Display_Record_Index),' /',num2str(measurement_Count)]...
                       ,['JV cycle no =',num2str(JV_current_Display_Record_Index),' /',num2str(JV_mmt_Count)]});

    set(handles.consider,'value',consider(current_Display_Record_Index));
    set (handles.cycle_no,'string',strcat(num2str(current_Display_Record_Index),' /',num2str(measurement_Count)));
end

if(JV_mmt_Count>0)
    if(JV_current_Display_Record_Index>1)
    JV_current_Display_Record_Index=JV_current_Display_Record_Index-1;
end
    axes(handles.JV_mmt_Store);
    hold on
    cla(handles.JV_mmt_Store);
    plot(handles.JV_mmt_Store,JV_Timestamp_Store(1:JV_Data_Length_Store(JV_current_Display_Record_Index),JV_current_Display_Record_Index), JV_Dist_Store((1:JV_Data_Length_Store(JV_current_Display_Record_Index)),JV_current_Display_Record_Index),'r');
    temp_length= length(find(JV_ECG_TS_Store(:,JV_current_Display_Record_Index)~=0));
    plot(handles.JV_mmt_Store,JV_ECG_TS_Store(1:temp_length,JV_current_Display_Record_Index),JV_ECG_Store(1:temp_length,JV_current_Display_Record_Index),'b');
    cycle_info_cell = get(handles.cycle_Info,'String');
    cycle_info_cell{end,1}= strcat('JV cycle no =',num2str(JV_current_Display_Record_Index),' /',num2str(JV_mmt_Count));
    set(handles.cycle_Info,'String',cycle_info_cell);
end




% --- Executes on button press in next_Record.
function next_Record_Callback(hObject, eventdata, handles)
% hObject    handle to next_Record (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global measurement_Count;
global consider;
global dia_Store;
global dist_Store;
global time_Stamp_Store;
global data_Length_Store;
global D_d_Store;
global del_D_Store;
global PAT_Store;
global current_Display_Record_Index;
global JV_mmt_Count;
global JV_current_Display_Record_Index;
global JV_Dist_Store;
global JV_Data_Length_Store;
global JV_Timestamp_Store;
global ECG_TS_Store;
global ECG_Store;
global JV_ECG_Store ;
global JV_ECG_TS_Store; 

if(measurement_Count>0)
    if(current_Display_Record_Index<measurement_Count)
        current_Display_Record_Index = current_Display_Record_Index+1;
    end
    axes(handles.mmt_Store);
    hold on
    cla(handles.mmt_Store);
    plot(handles.mmt_Store,time_Stamp_Store(1:data_Length_Store(current_Display_Record_Index),current_Display_Record_Index), dist_Store((1:data_Length_Store(current_Display_Record_Index)),current_Display_Record_Index),'r');
    temp_length= length(find(ECG_TS_Store(:,current_Display_Record_Index)~=0));
    plot(handles.mmt_Store,ECG_TS_Store(1:temp_length,current_Display_Record_Index),ECG_Store(1:temp_length,current_Display_Record_Index),'b');
    SBP = str2double(get(handles.SBP_Edit,'String'));
    DBP = str2double(get(handles.DBP_Edit,'String'));
    [beta_Arr, compliance_Arr, distensibility_Arr, best_beta, avg_Dia,avg_del_D,avg_PAT,SD_PAT, beta_Confidence] = calculate_Stiffness(del_D_Store(1:measurement_Count), D_d_Store(1:measurement_Count),PAT_Store(1:measurement_Count), SBP, DBP);
    set(handles.cycle_Info,'String',{['D_d =  ',num2str(D_d_Store(current_Display_Record_Index))],['Del_D = ', num2str(del_D_Store(current_Display_Record_Index))], ... 
                       ['Beta = ', num2str(beta_Arr(current_Display_Record_Index))],['PAT = ', num2str(PAT_Store(current_Display_Record_Index))],['cycle no =',num2str(current_Display_Record_Index),' /',num2str(measurement_Count)]...
                       ,['JV cycle no =',num2str(JV_current_Display_Record_Index),' /',num2str(JV_mmt_Count)]});
               
    set(handles.consider,'value',consider(current_Display_Record_Index));
    set (handles.cycle_no,'string',strcat(num2str(current_Display_Record_Index),' /',num2str(measurement_Count)));
end

if(JV_mmt_Count>0)
    if(JV_current_Display_Record_Index<JV_mmt_Count)
        JV_current_Display_Record_Index=JV_current_Display_Record_Index+1;
    end
    axes(handles.JV_mmt_Store);
    hold on
    cla(handles.JV_mmt_Store);
    plot(handles.JV_mmt_Store,JV_Timestamp_Store(1:JV_Data_Length_Store(JV_current_Display_Record_Index),JV_current_Display_Record_Index), JV_Dist_Store((1:JV_Data_Length_Store(JV_current_Display_Record_Index)),JV_current_Display_Record_Index),'r');
    temp_length= length(find(JV_ECG_TS_Store(:,JV_current_Display_Record_Index)~=0));
    plot(handles.JV_mmt_Store,JV_ECG_TS_Store(1:temp_length,JV_current_Display_Record_Index),JV_ECG_Store(1:temp_length,JV_current_Display_Record_Index),'b');
    cycle_info_cell = get(handles.cycle_Info,'String');
    cycle_info_cell{end,1}= strcat('JV cycle no =',num2str(JV_current_Display_Record_Index),' /',num2str(JV_mmt_Count));
    set(handles.cycle_Info,'String',cycle_info_cell);
end



% --- Executes on button press in consider.
function consider_Callback(hObject, eventdata, handles)
% hObject    handle to consider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of consider
global consider;
global c_consider;
global f_consider;
global disp_femoral;
global current_Display_Record_Index;


consider(current_Display_Record_Index)= get(handles.consider,'value');
if (disp_femoral ==1)
    f_consider(current_Display_Record_Index) = consider(current_Display_Record_Index);
elseif (disp_femoral ==0)
    c_consider(current_Display_Record_Index) = consider(current_Display_Record_Index);
end


% --- Executes on button press in result.
function result_Callback(hObject, eventdata, handles)
% hObject    handle to result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global c_PAT;
global f_PAT;
global PAT_Store;
global c_done;
global f_done;
global consider;
global c_consider;
global f_consider;
global c_HR;
global f_HR;


global c_measurement_Count;
global c_dia;
global c_dist;
global c_D_d;
global c_del_D;
global f_measurement_Count;
global f_dia;
global f_dist;
global f_D_d;
global f_del_D;
global live;
live = get(handles.real_Time_Opt,'Value');



length_CN = str2double(get(handles.length_CN_edit,'String'));
length_NF = str2double(get(handles.length_NF_edit,'String'));

if  (isnan(length_CN) ==1 || isnan(length_NF) ==1)
    msgbox('please enter distance values (length to calculate PWV)','error','error');
    
    return
end

length = length_NF - length_CN ;
SBP = str2double(get(handles.SBP_Edit,'String'));
DBP = str2double(get(handles.DBP_Edit,'String'));

auto = get(handles.auto,'value');

if (auto ==1)
    c_consider = ones(c_measurement_Count,1);
    f_consider = ones(f_measurement_Count,1);
end
if(c_measurement_Count >= 1)
    [beta_Arr, compliance_Arr, distensibility_Arr, best_beta, avg_Dia,avg_del_D,avg_PAT,SD_PAT, beta_Confidence] = calculate_Stiffness(c_del_D(1:c_measurement_Count), c_D_d(1:c_measurement_Count),c_PAT(1:c_measurement_Count), SBP, DBP);
    PAT_low = avg_PAT - SD_PAT;
    PAT_high = avg_PAT + SD_PAT;
    count =0;
    for (i=1 : c_measurement_Count)
        if (((PAT_low <= c_PAT(i)) && (c_PAT(i) <= PAT_high) &&(auto ==1)) || ((c_consider(i) == 1)&&(auto ==0)))
            c_consider(i) = 1;
            count = count +1;
            PAT_carotid(count) = c_PAT(i);
            
        else
            c_consider(i) = 0;
        end
    end
    [beta_Arr, compliance_Arr, distensibility_Arr, best_beta, avg_Dia,avg_del_D,c_avg_PAT,c_SD_PAT, beta_Confidence] = calculate_Stiffness(c_del_D(1:c_measurement_Count), c_D_d(1:c_measurement_Count),PAT_carotid(1:count), SBP, DBP);
else
    c_measurement_Count =0;
    c_avg_PAT = 0;
    c_SD_PAT =0;
    best_beta = 0;
    beta_Confidence =0;
    avg_Dia = 0;
end
if(f_measurement_Count >= 1)
    avg_PAT = mean(f_PAT(1:f_measurement_Count));
    SD_PAT = std(f_PAT(1:f_measurement_Count));
    PAT_low = avg_PAT - SD_PAT;
    PAT_high = avg_PAT + SD_PAT;
    f_count =0;
    for (i=1 : f_measurement_Count)
        if (((PAT_low <= f_PAT(i)) && (f_PAT(i) <= PAT_high)&&(auto ==1)) || ((f_consider(i) == 1)&&(auto ==0)))
            f_consider(i) = 1;
            f_count = f_count +1;
            PAT_femoral(f_count) = f_PAT(i);
        else
            f_consider(i) = 0;
        end
    end
    f_avg_PAT = mean(PAT_femoral(1:f_count));
    f_SD_PAT = std(PAT_femoral(1:f_count));
else
    f_measurement_Count =0;
    f_avg_PAT = 0;
    f_SD_PAT =0;
end
% set(handles.final_Results_Text,'String',{['Carotid'],['Best Beta = ',num2str(best_beta)], ['Beta Confidence = ',num2str(beta_Confidence)], ['Average D_d = ',num2str(avg_Dia)], ['PAT car = ',num2str(c_avg_PAT),' | ' ,num2str(c_SD_PAT)], ['PAT fem = ',num2str(f_avg_PAT),' | ' ,num2str(f_SD_PAT)]}); 
    

if ((f_measurement_Count ~= 0) && (c_measurement_Count ~= 0))    
    time_diff= f_avg_PAT - c_avg_PAT ;
    pwv = (length/100)/(time_diff/1000);

    set(handles.pwv_result,'String',{['Carotid Beta = ',num2str(best_beta)],['carotid PAT=',num2str(c_avg_PAT),' ms | ',num2str(c_SD_PAT)],['femoral PAT=',num2str(f_avg_PAT),' ms | ',num2str(f_SD_PAT)],['time delay=  ',num2str(time_diff),' ms'],['PWV =  ',num2str(pwv),' m/s']});
%     best_beta
    median_del_D=median(c_del_D(1:c_measurement_Count));
    delta_p= ((pwv * pwv)* 1050 * (((2* median_del_D * avg_Dia)+ (median_del_D * median_del_D))/(avg_Dia * avg_Dia))  /(1000 * 0.1333322));
    [beta_Arr, compliance_Arr, distensibility_Arr, best_beta_calcu, avg_Dia,avg_del_D,c_avg_PAT,c_SD_PAT, beta_Confidence] = calculate_Stiffness(c_del_D(1:c_measurement_Count), c_D_d(1:c_measurement_Count),PAT_carotid(1:count), (DBP + delta_p) , DBP);
    set(handles.pwv_result,'String',{['Carotid Beta = ',num2str(best_beta)],['carotid PAT=',num2str(c_avg_PAT),' ms | ',num2str(c_SD_PAT)],['femoral PAT=',num2str(f_avg_PAT),' ms | ',num2str(f_SD_PAT)],['time delay=  ',num2str(time_diff),' ms'],['PWV =  ',num2str(pwv),' m/s']...
        ,['del_P calcu = ',num2str(delta_p),' mm/Hg'], ['Car Beta calcu = ',num2str(best_beta_calcu)], ['car HR=',num2str(c_HR),'   fem HR=',num2str(f_HR)] });
    if (abs(c_HR - f_HR)>=10)
        msgbox('heart rate has varied (>10) between two site measurement','warning !!','warn');
    end


else
    set(handles.pwv_result,'String',{['Carotid Beta = ',num2str(best_beta)],['carotid PAT=',num2str(c_avg_PAT),' ms | ',num2str(c_SD_PAT)],['femoral PAT=',num2str(f_avg_PAT),' ms | ',num2str(f_SD_PAT)],['car HR=',num2str(c_HR),'fem HR=',num2str(f_HR)],[ 'measurement not yet completed']});
    time_diff ='measurement not yet completed';
    pwv ='--';
    delta_p = 'measurement not yet completed';
    best_beta_calcu = 'measurement not yet completed';
end
% set(handles.cycle_Info,'String',{['D_d =  ',num2str(D_d_Store(current_Display_Record_Index))],['Del_D = ', num2str(del_D_Store(current_Display_Record_Index))], ...
%                        ['Beta = ', num2str(beta_Arr(current_Display_Record_Index))],['PAT = ', num2str(PAT_Store(current_Display_Record_Index))],['cycle no =',num2str(current_Display_Record_Index),' /',num2str(measurement_Count) ]});
set (handles.beta_disp,'string',num2str(best_beta));
set (handles.pwv_disp,'string',strcat( num2str(pwv),' m/s'));
current_File_Name = get(handles.US_File_Addr,'String');
if live==1
    patient_name =get(handles.file_name,'String');    
    
else
    [pathstr, name_str, file_ext]=fileparts(current_File_Name);
    patient_name = name_str(1:end-17);
end
age=str2double(get(handles.age_yrs,'String'));
hight=str2double(get(handles.hight_cm,'String'));
weight=str2double(get(handles.weight_kg,'String'));

contents_gender = cellstr(get(handles.gender_popup,'String'));
gender=contents_gender{get(handles.gender_popup,'Value')};
contents_posture = cellstr(get(handles.posture_popup,'String'));
posture=contents_posture{get(handles.posture_popup,'Value')};
result_table= {'name:',patient_name;'age',age;'hight',hight; 'weight',weight; 'gender',gender;'posture',posture;...
                'SBP' ,SBP; 'DBP' , DBP; 'length carotid-Snotch (cm)',length_CN;...
                'length Snotch-femoral(cm)', length_NF; 'Carotid Beta', best_beta; 'carotid PAT(ms)',c_avg_PAT;...
                'femoral PAT(ms)',f_avg_PAT; 'time delay(ms)', time_diff; 'PWV (m/s)' ,  pwv; 'del_P calcu (mm/Hg)',delta_p;...
                'Car Beta calcu',best_beta_calcu; 'carotid HR=',c_HR;'fememoral HR=',f_HR }%  'SD carotid PAT(ms)',c_SD_PAT; 'SD femoral PAT(ms)',f_SD_PAT }

save_file = get(handles.save_file,'value');
current_File_Name = strcat(current_File_Name(1:end-8),'.xls');
if save_file == 0
    choice = questdlg({'Would you like to save result?'}, ...
                            'save result file', ...
                            'Yes','No','Yes');
else
    choice = 'Yes';
end
if((strcmp(choice,'Yes')==1) || save_file ==1)

    xlswrite(current_File_Name,result_table);
%     current_File_Name=get(handles.US_File_Addr,'String'); %
%     if((strcmp(current_File_Name(end-2:end),'mat')==1)|| (save_file ==1))
%         f_consider_backup = f_consider;
%         c_consider_backup = c_consider;
%         prev_handles = handles;
%         load(current_File_Name);
%         saved_handles = handles;
%         handles = prev_handles;
%         SBP = str2double(get(handles.SBP_Edit,'String'));
%         DBP = str2double(get(handles.DBP_Edit,'String'));
%         length_CN=str2double(get(handles.length_CN_edit,'String'));
%         length_NF=str2double(get(handles.length_NF_edit,'String'));
%         f_consider= f_consider_backup  ;
%         c_consider = c_consider_backup  ;
%         current_File_Name=get(handles.US_File_Addr,'String')
%         handles = saved_handles;
%         save (current_File_Name)
%         handles = prev_handles;
%     end

end
live = 0;


      


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
global pulse_filt_on;
global calibration_Data;
global calibration_Data_loc;
if (pulse_filt_on == -1)
    if(isempty(calibration_Data_loc) || strcmp(calibration_Data_loc ,'')==1)
        [FileName,PathName] = uigetfile('*.tcf','Select the Compensation File');
        try
            calibration_Data = dlmread(strcat(PathName,FileName));
        catch
            msgbox('Error in reading calibration file_2','error','error');
            calibration_Data_loc ='';
        end
        
        calibration_Data_loc =strcat(PathName,FileName);
        set(hObject,'String','Disable Transducer Compensation');
        
        fid_sys= fopen('system_default.csv', 'r+');
        system_default = fgetl(fid_sys);
        [device_name, s_loc, s_files, auto_res,  ecg_on_str,~,~] = strread(system_default, '%s %s %s %s %s %s %s', 'delimiter',','); %system_default = strcat(device,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',')

        system_default = strcat(device_name,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',');
        fclose (fid_sys);
        fid_sys = fopen('system_default.csv','w');
        fwrite(fid_sys,system_default{1});
        fclose (fid_sys);
        
    elseif (isempty(calibration_Data))
        try
            calibration_Data = dlmread(calibration_Data_loc);
        catch
            msgbox('Error in reading calibration file_3','error','error');    
        end
        set(hObject,'String','Disable Transducer Compensation');
    else
        calibration_Data=[];
        set(hObject,'String','Enable Transducer Compensation');
    end

elseif (pulse_filt_on == 0)
    pulse_filt_on = 1;
    set(handles.TCF_But,'String','off pulse filt');
elseif (pulse_filt_on == 1)
    pulse_filt_on = 0;
    set(handles.TCF_But,'String','on pulse filt');
end


function SBP_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to SBP_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SBP_Edit as text
%        str2double(get(hObject,'String')) returns contents of SBP_Edit as a double


% --- Executes during object creation, after setting all properties.
function SBP_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SBP_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DBP_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to DBP_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DBP_Edit as text
%        str2double(get(hObject,'String')) returns contents of DBP_Edit as a double


% --- Executes during object creation, after setting all properties.
function DBP_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DBP_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in read_File_Opt.
function read_File_Opt_Callback(hObject, eventdata, handles)
% hObject    handle to read_File_Opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of read_File_Opt


% --- Executes on button press in real_Time_Opt.
function real_Time_Opt_Callback(hObject, eventdata, handles)
% hObject    handle to real_Time_Opt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of real_Time_Opt


% --- Executes during object creation, after setting all properties.
function Analyze_But_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Analyze_But (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes on button press in speed_Test_But.
function speed_Test_But_Callback(hObject, eventdata, handles)
% hObject    handle to speed_Test_But (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global samples_Per_Frame;
global live;
frame_Rate=0;
read_File_Opt = get(handles.read_File_Opt,'Value');
live = get(handles.real_Time_Opt,'Value');

if(read_File_Opt ==1 && live ==0)
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
    
    set(handles.speed_Indicator_Text,'String','Calculating Frame Rate. Please Wait ...')
    drawnow
    frame_Num=1;
    tic;
    while(frame_Num < min(1000,total_Frames))
        [current_Frame,TS] = load_LV_Frame(US_Fid, samples_Per_Frame,frame_Num,1);
        frame_Num = frame_Num+1;
    end
    frame_Rate = min(1000,total_Frames)/toc;
    set(handles.speed_Indicator_Text,'String',['Frame Rate : ', num2str(frame_Rate)])
    fclose(US_Fid);

end


function length_CN_edit_Callback(hObject, eventdata, handles)
% hObject    handle to length_CN_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of length_CN_edit as text
%        str2double(get(hObject,'String')) returns contents of length_CN_edit as a double


% --- Executes during object creation, after setting all properties.
function length_CN_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to length_CN_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gen_calib.
function gen_calib_Callback(hObject, eventdata, handles)
% hObject    handle to gen_calib (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global live;
live = get(handles.real_Time_Opt,'Value');
if (live ==1)
global US_Sampling_Rate;
global samples_Per_Frame;
global transducer_Frequency;
global calibration_Data;
global device;
global stop_All;
global calibration_Data_loc;
stop_All=0;
US_Sampling_Rate = 100000000;
samples_Per_Frame = 5195;
transducer_Frequency = 500000;
cla(handles.US_Frame_Plot,'reset');
% ylim(handles.US_Frame_Plot,[-4 4]);
% xlim(handles.US_Frame_Plot,[1 samples_Per_Frame]);
hold(handles.US_Frame_Plot,'on')

total_Frames = 1000;
    [start_time] = niScope_Link(1,samples_Per_Frame, device,'a');
    tic
    [current_Frame,~] = niScope_Link(2,samples_Per_Frame,0);
    current_Frame = current_Frame - mean(current_Frame);
    plot(handles.US_Frame_Plot,(1:samples_Per_Frame),current_Frame);
    [ calibration_Done, cut_Point, calibration_Frame, Array_for_File] = gen_Calibration_Data(current_Frame, US_Sampling_Rate, samples_Per_Frame, total_Frames, 1);
    while(calibration_Done ~=1)
        [current_Frame,~] = niScope_Link(2,samples_Per_Frame,0);
        current_Frame = current_Frame - mean(current_Frame);
        [ calibration_Done, cut_Point, calibration_Frame, Array_for_File] = gen_Calibration_Data(current_Frame, US_Sampling_Rate, samples_Per_Frame, total_Frames, 0);
        cla(handles.US_Frame_Plot);
        plot(handles.US_Frame_Plot,(1:samples_Per_Frame),current_Frame);
        drawnow;
        if (stop_All==1)
        break
        end
    end
    toc
    niScope_Link(3,samples_Per_Frame,0);
    if (stop_All==0)
        
%     cla(handles.US_Frame_Plot);
    csvwrite('dev_newprobe.tcf',Array_for_File);
    calibration_Data_loc = 'dev_newprobe.tcf';
    
    fid_sys= fopen('system_default.csv', 'r+');
    system_default = fgetl(fid_sys);
    [device_name, s_loc, s_files, auto_res,  ecg_on_str,~,~] = strread(system_default, '%s %s %s %s %s %s %s', 'delimiter',','); %system_default = strcat(device,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',')

    system_default = strcat(device_name,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',');
    fclose (fid_sys);
    fid_sys = fopen('system_default.csv','w');
    fwrite(fid_sys,system_default{1});
    fclose (fid_sys);
    
    calibration_Data = Array_for_File;
    set(handles.TCF_But,'String','Disable Transducer Compensation');
    end
end


function file_name_Callback(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of file_name as text
%        str2double(get(hObject,'String')) returns contents of file_name as a double


% --- Executes during object creation, after setting all properties.
function file_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function length_NF_edit_Callback(hObject, eventdata, handles)
% hObject    handle to length_NF_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of length_NF_edit as text
%        str2double(get(hObject,'String')) returns contents of length_NF_edit as a double


% --- Executes during object creation, after setting all properties.
function length_NF_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to length_NF_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_file.
function save_file_Callback(hObject, eventdata, handles)
% hObject    handle to save_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_file
fid_sys= fopen('system_default.csv', 'r+');
system_default = fgetl(fid_sys);
[device_name, s_loc, s_files, auto_res,  ecg_on_str,calibration_Data_loc,~] = strread(system_default, '%s %s %s %s %s %s %s', 'delimiter',','); %system_default = strcat(device,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',')

s_files = num2str( get(hObject,'Value'));
system_default = strcat(device_name,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',');
fclose (fid_sys);
fid_sys = fopen('system_default.csv','w');
fwrite(fid_sys,system_default{1});
fclose (fid_sys);


% --- Executes on button press in ECG_on.
function ECG_on_Callback(hObject, eventdata, handles)
% hObject    handle to ECG_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ECG_on
global ecg_on;
if (get(handles.real_Time_Opt,'value')==1)
    start = get(handles.ECG_on,'value');
    if (start==1)
%         ECG_Link(1,'test.ecg');
        ecg_on=1;
    else 
        ecg_on=0;
%         ECG_Link(3,0);
    end
    fid_sys= fopen('system_default.csv', 'r+');
    system_default = fgetl(fid_sys);
    [device_name, s_loc, s_files, auto_res,  ecg_on_str,calibration_Data_loc,~] = strread(system_default, '%s %s %s %s %s %s %s', 'delimiter',','); %system_default = strcat(device,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',')

    ecg_on_str = num2str( ecg_on);
    system_default = strcat(device_name,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',');
    fclose (fid_sys);
    fid_sys = fopen('system_default.csv','w');
    fwrite(fid_sys,system_default{1});
    fclose (fid_sys);
else 
    ecg_on = get(handles.ECG_on,'value'); 
end



% --- Executes on slider movement.
function ecg_slider_Callback(hObject, eventdata, handles)
% hObject    handle to ecg_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global prev_Time_Stamp;
global ecg_data;
global tstamp_ecg;
global ecg_disp_start;
global tstamp_Rpeak;
global carotid_wave_start;
global JV_wave_start;
global long_dist_store;
global dist_time_store;
global b_ecg_l;
global a_ecg_l;

ecg_disp_start = round(get(handles.ecg_slider,'value'));
ecg_disp_end = ecg_disp_start + 6999;
ECG_Data_Ch = double(( ecg_data - min(ecg_data(ecg_disp_start:ecg_disp_end))) /(max(ecg_data(ecg_disp_start:ecg_disp_end))- min(ecg_data(ecg_disp_start:ecg_disp_end))));
ECG_Data_Ch  = filtfilt(b_ecg_l, a_ecg_l,  ECG_Data_Ch);               
cla(handles.ECG_Plot);
xlim(handles.ECG_Plot,[ecg_disp_start*2 (ecg_disp_end*2)]);
plot(handles.ECG_Plot,(tstamp_ecg),ECG_Data_Ch);

dist_disp_start = 1 + length (find(dist_time_store < tstamp_ecg(ecg_disp_start)));
dist_disp_end = length (find(dist_time_store < tstamp_ecg(ecg_disp_end)));
if (dist_disp_end > dist_disp_start)
tempdata = ((long_dist_store)- min(long_dist_store(dist_disp_start:dist_disp_end)))/(max(abs(long_dist_store(dist_disp_start:dist_disp_end)))- min(long_dist_store(dist_disp_start:dist_disp_end))) +1;% ))-min(tempdata);

plot(handles.ECG_Plot,(dist_time_store ),tempdata,'r');
end
plot(handles.ECG_Plot,[(tstamp_Rpeak),( tstamp_Rpeak )],[2,0],'c');

plot(handles.ECG_Plot,[(carotid_wave_start)',( carotid_wave_start)'],[2,0],'k');

plot(handles.ECG_Plot,[(JV_wave_start)',( JV_wave_start)'],[2,0],'m');

ecg_disp_start;



% --- Executes during object creation, after setting all properties.
function ecg_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ecg_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in auto.
function auto_Callback(hObject, eventdata, handles)
% hObject    handle to auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto
fid_sys= fopen('system_default.csv', 'r+');
system_default = fgetl(fid_sys);
[device_name, s_loc, s_files, auto_res,  ecg_on_str,calibration_Data_loc,~] = strread(system_default, '%s %s %s %s %s %s %s', 'delimiter',','); %system_default = strcat(device,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',')

auto_res = num2str( get(hObject,'Value'));
system_default = strcat(device_name,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',');
fclose (fid_sys);
fid_sys = fopen('system_default.csv','w');
fwrite(fid_sys,system_default{1});
fclose (fid_sys);


% --- Executes on button press in save_loc.
function save_loc_Callback(hObject, eventdata, handles)
% hObject    handle to save_loc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder_name = uigetdir(get(handles.dir_loc,'string'),'select Data save location');
if (folder_name ~=0)
    set(handles.dir_loc,'string',folder_name);
    
else
    folder_name =pwd;
    set(handles.dir_loc,'string',pwd);
end
fid_sys= fopen('system_default.csv', 'r+');
system_default = fgetl(fid_sys);
[device_name, s_loc, s_files, auto_res,  ecg_on_str,calibration_Data_loc,~] = strread(system_default, '%s %s %s %s %s %s %s', 'delimiter',','); %system_default = strcat(device,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',')

s_loc = folder_name;
system_default = strcat(device_name,',',s_loc,',',s_files,',',auto_res,',',ecg_on_str,',',calibration_Data_loc,',');
fclose (fid_sys);
fid_sys = fopen('system_default.csv','w');
fwrite(fid_sys,system_default{1});
fclose (fid_sys);




% --- Executes during object creation, after setting all properties.
function dir_loc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dir_loc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function reading_no_Callback(hObject, eventdata, handles)
% hObject    handle to reading_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reading_no as text
%        str2double(get(hObject,'String')) returns contents of reading_no as a double


% --- Executes during object creation, after setting all properties.
function reading_no_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reading_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Save_continous.
function Save_continous_Callback(hObject, eventdata, handles)
% hObject    handle to Save_continous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global bg_colour;
global save_conti_RF;
global running ;
global c_consider;
global f_consider;
if running==1
if save_conti_RF ==1
    save_conti_RF =0;
    set(handles.Save_continous,'string','save');
    set(handles.Save_continous,'BackgroundColor',bg_colour);
elseif save_conti_RF ==0
    save_conti_RF =1;
    
    set(handles.Save_continous,'string','saving ON');
    set(handles.Save_continous,'BackgroundColor','g');
end
elseif running==0
current_File_Name=get(handles.US_File_Addr,'String'); %
if(strcmp(current_File_Name(end-2:end),'mat')==1)
    f_consider_backup = f_consider;
    c_consider_backup = c_consider;
    prev_handles = handles;
    load(current_File_Name);
    old_handles = handles;
    handles = prev_handles;
    SBP = str2double(get(handles.SBP_Edit,'String'));
    DBP = str2double(get(handles.DBP_Edit,'String'));
    length_CN=str2double(get(handles.length_CN_edit,'String'));
    length_NF=str2double(get(handles.length_NF_edit,'String'));
    age=str2double(get(handles.age_yrs,'String'));
    hight=str2double(get(handles.hight_cm,'String'));
    weight=str2double(get(handles.weight_kg,'String'));

    contents_gender = cellstr(get(handles.gender_popup,'String'));
    gender=contents_gender{get(handles.gender_popup,'Value')};
    contents_posture = cellstr(get(handles.posture_popup,'String'));
    posture=contents_posture{get(handles.posture_popup,'Value')};
    f_consider= f_consider_backup  ;
    c_consider = c_consider_backup  ;
    current_File_Name=get(handles.US_File_Addr,'String')
    
    handles = old_handles ;
    save (current_File_Name)
    handles = prev_handles;
end
end

% --- Executes on selection change in posture_popup.
function posture_popup_Callback(hObject, eventdata, handles)
% hObject    handle to posture_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns posture_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from posture_popup
contents = cellstr(get(hObject,'String'));
posture=contents{get(hObject,'Value')};


% --- Executes during object creation, after setting all properties.
function posture_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posture_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in gender_popup.
function gender_popup_Callback(hObject, eventdata, handles)
% hObject    handle to gender_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns gender_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gender_popup
contents = cellstr(get(hObject,'String'))
get(hObject,'Value')
gender=contents{get(hObject,'Value')}
try
    
catch
end

% --- Executes during object creation, after setting all properties.
function gender_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gender_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hight_cm_Callback(hObject, eventdata, handles)
% hObject    handle to hight_cm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hight_cm as text
%        str2double(get(hObject,'String')) returns contents of hight_cm as a double


% --- Executes during object creation, after setting all properties.
function hight_cm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hight_cm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function weight_kg_Callback(hObject, eventdata, handles)
% hObject    handle to weight_kg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of weight_kg as text
%        str2double(get(hObject,'String')) returns contents of weight_kg as a double


% --- Executes during object creation, after setting all properties.
function weight_kg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to weight_kg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function age_yrs_Callback(hObject, eventdata, handles)
% hObject    handle to age_yrs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of age_yrs as text
%        str2double(get(hObject,'String')) returns contents of age_yrs as a double


% --- Executes during object creation, after setting all properties.
function age_yrs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to age_yrs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cycle_no_Callback(hObject, eventdata, handles)
% hObject    handle to cycle_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cycle_no as text
%        str2double(get(hObject,'String')) returns contents of cycle_no as a double


% --- Executes during object creation, after setting all properties.
function cycle_no_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cycle_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on cycle_no and none of its controls.
function cycle_no_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to cycle_no (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in new_patient.
function new_patient_Callback(hObject, eventdata, handles)
% hObject    handle to new_patient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.file_name,'String',''); 
set(handles.file_name,'Enable' ,'on');
set(handles.reading_no,'String','1');
set(handles.SBP_Edit,'String','-');
set(handles.DBP_Edit,'String','-');
set(handles.length_CN_edit,'String','-');
set(handles.length_NF_edit,'String','-');

set(handles.age_yrs,'String','-');
set(handles.hight_cm,'String','-');
set(handles.weight_kg,'String','-');

set(handles.gender_popup,'Value',1);
set(handles.posture_popup,'Value',1);
[cleared_all]=clear_all_prev_data();


%%%%%%%%%General Functions other than GUI%%%%%%%%%%%%%%%%%%%
function [state, wall1, wall2, on_Carotid, time_Array, distension_Array, dia_Array, full_Cycle_Time_Array, full_Cycle_Dist_Array, full_Cycle_Lumen_Dia_Array, del_D, D_d,full_Cycle_ECG,full_Cycle_ECG_Time,PAT] = track_Walls_wth_ECG(frame, Wall1, Wall2, pls_Rate,TS)
global last_Track_Frame;
global PW_Wave;
global DW_Wave;
global TS_Wave;
global pulse_Rate;
global current_State;
global samples_Per_Frame;
global points_Per_MM;
global bad_Signal_Quality_Count;
global dia_Wave;
global dia_Confidence;
global prev_Time_Stamp;
global ecg_on;

distension_Array=0;
dia_Array=0;
full_Cycle_Time_Array = 0;
full_Cycle_Dist_Array =0;
full_Cycle_Lumen_Dia_Array = 0; 
full_Cycle_ECG =0;
full_Cycle_ECG_Time=0;
sub_TS_Wave =0;
del_D = 0; 
D_d = 0;
PAT =0;
wall1=0;
wall2=0;
on_Carotid=1;
ecg_sub_Wave=0;
sub_tstamp_ecg=0;
    
if(current_State ==2)
    pulse_Rate = pls_Rate;
    PW_Wave = [Wall1];
    DW_Wave = [Wall2];
    TS_Wave = [prev_Time_Stamp];
    current_State=3;
    bad_Signal_Quality_Count=0;
    dia_Wave =[0];
    dia_Confidence =[0];
     [dia_Wave,dia_Confidence] = find_Dia_L_Fit(frame,PW_Wave(end,1),DW_Wave(end,1),samples_Per_Frame,points_Per_MM);
end
PW_Wave = [PW_Wave; 0];
DW_Wave = [DW_Wave; 0];
dia_Wave = [dia_Wave; 0];
TS_Wave = [TS_Wave; TS]; % add time stamp here
dia_Confidence = [dia_Confidence; 0];
[valid_Position, PW_Wave(end,1), DW_Wave(end,1)] = track_Walls(frame,last_Track_Frame, PW_Wave(end-1,1), DW_Wave(end-1,1),samples_Per_Frame,points_Per_MM);

if(valid_Position == 0)
    current_State = 5;
else
    if(check_Walls_Present(frame,PW_Wave(end),DW_Wave(end))== 0)
        bad_Signal_Quality_Count = bad_Signal_Quality_Count+1;
        if(bad_Signal_Quality_Count >= 20)
            bad_Signal_Quality_Count = 0;
            current_State = 5;
        end
    else 
        bad_Signal_Quality_Count =0;
        wall1= PW_Wave(end);
        wall2= DW_Wave(end);
         [dia_Wave(end,1),dia_Confidence(end,1)] = find_Dia_L_Fit(frame,PW_Wave(end,1),DW_Wave(end,1),samples_Per_Frame,points_Per_MM);
        if (ecg_on==1)
            [state_Change, on_Carotid, change_In_Dia, diastolic_Diameter, distension_Wave,sub_TS_Wave,ecg_sub_Wave,sub_tstamp_ecg,PAT] = check_Cycles_ECG_Method();
        else
            [state_Change, on_Carotid, change_In_Dia, diastolic_Diameter, distension_Wave,sub_TS_Wave] = check_Cycles_ZC_Method();
        end
        current_State = current_State+state_Change;
        if(current_State==4)
            full_Cycle_Dist_Array = distension_Wave;
            if (sub_TS_Wave(1) == 0 && sub_TS_Wave(end)==0)
                full_Cycle_Time_Array = (1:length(distension_Wave))/pulse_Rate;
            else
                full_Cycle_Time_Array = sub_TS_Wave - sub_TS_Wave(1);% (1:length(distension_Wave))/pulse_Rate;
            end
            del_D =  change_In_Dia;
            D_d = diastolic_Diameter;
            full_Cycle_ECG = ecg_sub_Wave;
            full_Cycle_ECG_Time = sub_tstamp_ecg -sub_TS_Wave(1);
        end
    end
end
distension_Array = ((DW_Wave - PW_Wave))/points_Per_MM;%-mean(DW_Wave - PW_Wave)
dia_Array = dia_Wave;
time_Array = TS_Wave;
%time_Array = (1:length(dia_Wave))/pulse_Rate;
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

consecutive_Find_Limit = 2;
second_State_Multiplier = 3;

if(reset ==1)
    storage_Mat_Column_Nos = 15;
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
global b_Corr;
global a_Corr;
line_X = zeros(size(line1));
corr_Out = zeros(size(line1));
line_X= line1.*line2;
X1 = filtfilt(b_Corr, a_Corr, line_X);
corr_Out = X1(:,1);



function [ peak_Array] = peaks_3mm_Apart(inp_Signal,points_Per_MM )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
inp_Signal = inp_Signal/max(abs(inp_Signal));
global min_Wall_Dia;
global max_Wall_Dia;
global Wall_Lumen_SNR;
T = 0.01;
P = [0 0];
peak_No=0;
no_Of_Pairs =0;
peak_Array_Temp = zeros(10,2);
count = 0;
peak_Array=[];

while(T<0.5)
   [peak_Locs,peak_Vals, valley_Locs, valley_Vals, peak_No, valley_No] = peak_detect(inp_Signal,T,length(inp_Signal));
   P = [peak_Locs(1:peak_No) , peak_Vals(1:peak_No)];
   
%     peak_No = size(P,1);
    if(peak_No==0)
        peak_Array = [];
        return
    elseif(peak_No >2 && peak_No< 10)
        break;
    elseif(peak_No < 3)
        T = T/5;
    else
        T = T*2;
    end
    count = count+1;
    if(count >10)
        break;
    end
end

for i = 1:peak_No-1
    if ((abs(P(i+1,1) - P(i,1)) > min_Wall_Dia*points_Per_MM) &&  (abs(P(i+1,1) - P(i,1)) <max_Wall_Dia*points_Per_MM))
        Val_N = inp_Signal(P(i,1));
        Val_F = inp_Signal(P(i+1,1));
        Val_M = inp_Signal(ceil((P(i,1)+P(i+1,1))/2));
        
        if (20*log(Val_N/Val_M)>Wall_Lumen_SNR && 20*log(Val_F/Val_M)>Wall_Lumen_SNR)   % Check SNR for each wall is greater than 20 dB
            no_Of_Pairs = no_Of_Pairs+1;
            peak_Array_Temp(no_Of_Pairs,:) = [P(i,1) P(i+1,1)];
        end
    end
end
peak_Array = peak_Array_Temp([1:no_Of_Pairs],:);


% function [maxtab, mintab] = peakdet(v, delta)
% global maxtab_Mem;
% global mintab_Mem;
% 
% S = size(v);
% if(S(1)*S(2)==0)
%     maxtab=[0,0];
%     mintab =[0,0];
%     return
% end
% max_Count =0;
% min_Count =0;
% 
% mn = v(1,1); mx= v(1,1); mxpos=1;mnpos=1;
% lookformax = 1;
% 
% for i=1:length(v)
%   if lookformax
%     if v(i) > mx, mx = v(i); mxpos =i;
%     elseif v(i) < mx-delta
%         if(mxpos ~=1)
%           max_Count = max_Count+1;
%           maxtab_Mem(max_Count,:) = [mxpos mx];
%         end
%       mn = v(i); mnpos = i;
%       lookformax = 0;
%     end  
%   else
%     if v(i) < mn, mn = v(i); mnpos =i;
%     elseif v(i) > mn+delta
%         if(mnpos~=S(1))
%           min_Count = min_Count+1;
%           mintab_Mem(min_Count,:) = [mnpos mn];
%         end
%       mx = v(i); mxpos = i;
%       lookformax = 1;
%     end
%   end
% end
% maxtab = maxtab_Mem([1:max_Count],:);
% mintab = mintab_Mem([1:min_Count],:);


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
if(C1<-0.5)
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
    C = corr(W1_Mean_Sub,W2_Mean_Sub,'type','spearman');
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
    r_Auto1 = real(ifft(fft(sig1,nfft) .* conj(fft(sig1,nfft))));
    r_Auto2 = real(ifft(fft(sig2,nfft) .* conj(fft(sig2,nfft))));
    Auto_Corr_Max = max(r_Auto1(1), r_Auto2(1));
    Corr_Sig = Corr_Sig/Auto_Corr_Max;
end


function [ frame,TS ] = load_LV_Frame(fid,  no_Of_Pts_Per_Frame, fNo,TS_on)
global calibration_Data;
global ah;
global bh;
global al;
global bl;
global pulse_filt_on;

frewind(fid);
fseek(fid,no_Of_Pts_Per_Frame*(fNo-1)*4,'bof');
[frame(:,1),count]=fread(fid,no_Of_Pts_Per_Frame,'*float','b');


if(count <no_Of_Pts_Per_Frame) 
    disp (strcat('failed to read frame number ',int2str(fNo)));
    return;
end
if(TS_on ==1)
    TS = frame(1);
    frame(1) =0;
else
    TS =0;
end
frame(:,1)= double(frame - mean(frame)) ;
if((length(calibration_Data)==length(frame)+3)&& pulse_filt_on==-1)
    hard_Limit_Point = calibration_Data(3);
    transducer_Response_Frame(:,1) = calibration_Data(4:end);
    frame(:,1) = [zeros(hard_Limit_Point,1); frame(hard_Limit_Point+1:end)-transducer_Response_Frame((hard_Limit_Point+1:end))];
else
    hard_Limit_Point=0;
    transducer_Response_Frame(:,1) = zeros(size(frame));
end
    

% frame = filtfilt(bh, ah, frame);

function [a,data,tstamp] = load_ecg_data(ecg_Fid, ecg_info,TS)

global ecg_data;
global tstamp_ecg;
a =0;
data = zeros(1,1);
tstamp = zeros(1,1);
if (tstamp_ecg(end) > TS)
    return
end
if ((TS - tstamp_ecg(end)) > 1000)
    no_of_data = 1000;
elseif ((TS - tstamp_ecg(end) ) > 500)
    no_of_data = 500;
else
    no_of_data = 100;
end
frewind(ecg_Fid);
fseek(ecg_Fid,ecg_info(28)*4,'bof');
[read_data,count]=fread(ecg_Fid,no_of_data,'*float','b');
a = count/2;
for (i = 0:((count/2)-1))
data(i+1,1) = read_data(i*2 +1);
tstamp (i+1,1) = read_data((i+1)*2);
end


    
function [ shift_Mat ] = gen_Shift_Mat( dim )
shift_Mat = zeros(dim,dim);
shift_Mat = eye(dim);
shift_Mat = circshift(shift_Mat,[0 -1]);
shift_Mat(1,dim)=0;
shift_Mat = shift_Mat';


function [ status ] = check_Walls_Present( frame,Wall1, Wall2)
global points_Per_MM;
global min_Wall_Dia;
global max_Wall_Dia;
global Wall_Lumen_SNR;
status =0;

if((Wall2-Wall1)/points_Per_MM < min_Wall_Dia || (Wall2-Wall1)/points_Per_MM > max_Wall_Dia || Wall1/points_Per_MM< 3 || (length(frame)-Wall2)/points_Per_MM < 3)
    status =0;
else
    Slot1= frame(Wall1-round(points_Per_MM/2):Wall1+round(points_Per_MM/2));
    Slot2= frame(Wall2-round(points_Per_MM/2):Wall2+round(points_Per_MM/2));
    PW_Slot_Avg = sum(abs(Slot1))/(2*round(points_Per_MM/2)+1);
    DW_Slot_Avg = sum(abs(Slot2))/(2*round(points_Per_MM/2)+1);
    Mid_Slot = frame(((round((Wall1+Wall2)/2))-points_Per_MM):((round((Wall1+Wall2)/2))+points_Per_MM));
    Mid_Slot_Avg = (sum(abs(Mid_Slot)))/(2*points_Per_MM+1);

    if(20*log(PW_Slot_Avg/Mid_Slot_Avg) > Wall_Lumen_SNR && 20*log(DW_Slot_Avg/Mid_Slot_Avg) > Wall_Lumen_SNR)
        status =1;
    else
        status =0;
    end
end


function [ diameter_Lumen,confidence] = find_Dia_L_Fit(frame,Wall1,Wall2,No_Of_Samples,points_Per_mm)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global last_Dia;
global b_L_Fit;
global a_L_Fit;
mid_Pt = round((Wall1+Wall2)/2);
confidence=0;
diameter_Lumen =0;
if(Wall1 < 2*points_Per_mm || mid_Pt < 2*points_Per_mm || Wall2 > No_Of_Samples - 2*points_Per_mm)
    return
end

frame = filtfilt(b_L_Fit, a_L_Fit, frame.*frame);

left_Half = frame ([Wall1-round(1.5*points_Per_mm):mid_Pt],1);
left_Half = flipud(left_Half);
right_Half = frame ([mid_Pt:Wall2+round(1.5*points_Per_mm)], 1);

left_Half_Norm = left_Half/max(abs(left_Half));
right_Half_Norm = right_Half/max(abs(right_Half));

[~,P_R] = max(right_Half_Norm);
[~,P_L] = max(left_Half_Norm);

left_Half_Cut = left_Half_Norm([1: P_L(1,1)],1);
[left_Intima_Pt, sharpness_Index_Left ,left_L,L_trial_Count ] = fitL(left_Half_Cut); %Inflection point, Error, Curve, trial_Count
right_Half_Cut = right_Half_Norm([1:P_R(1,1)],1);
[right_Intima_Pt,sharpness_Index_Right, right_L,R_trial_Count ] = fitL(right_Half_Cut);
confidence = min(1/sharpness_Index_Left,1/sharpness_Index_Right);

diameter_Lumen = (left_Intima_Pt+ right_Intima_Pt)/points_Per_mm;



% function [ inflection_X, L_Fitted, error] = Fit_L(well_Half)
% %UNTITLED2 Summary of this function goes here
% %   Detailed explanation goes here
% 
% Tot_Len = length(well_Half);
% L = 1*[1:Tot_Len];
% L_Fitted = L;
% error = 100;
% inflection_X = Tot_Len;
% 
% if(Tot_Len <=1)
%     return
% end
% 
% well_Half = well_Half-min(well_Half);
% well_Half = well_Half/(max(well_Half));
% Jump = 1;
% x= round(2*Tot_Len/3);
% trial_Count = 0;
% 
% while(1)
%     trial_Count = trial_Count+1;
%     x = x- Jump;
%     if(x<1  || x>Tot_Len-1)
%         break;
%     end
%     L1 = L_Gen(well_Half(1),x, well_Half(x),well_Half(Tot_Len),Tot_Len);
%     Err1 = sum((L1-well_Half).^2);
%     if Err1 < error
%         inflection_X = x;
%         L_Fitted = L1;
%         error= Err1;
%         if(error <0.5)
%             break;
%         end
%     end
%     L2 = L_Gen(well_Half(1),x+1, well_Half(x),well_Half(Tot_Len),Tot_Len);
%     Err2 = sum((L2-well_Half).^2);
%     Jump = round(Err1/(Err2-Err1));
%     if(trial_Count>=20)
%         break;
%     end
% end
% 
% 
% 
% 
% 
% function [ L_Ret ] = L_Gen(yb,xm,ym,ye,total_Pts )
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% xb=1;
% xe= total_Pts;
% line1 = linspace(yb,ym,xm-xb+1)';
% line2 = linspace(ym,ye,xe-xm+1)';
% L_Ret = [line1(1:length(line1)-1); line2];


function [ state_Change,  on_Carotid, Del_D, D_d, dist_Wave,sub_TS_Wave] = check_Cycles_ZC_Method()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global PW_Wave ;
global DW_Wave;
global dia_Wave;
global TS_Wave;
global points_Per_MM;
global pulse_Rate;
global dia_Confidence;
global TS_on;

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
sub_TS_Wave = zeros(1,1);
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
   sub_TS_Wave =TS_Wave((positive_Zero_Crossings(1,1):positive_Zero_Crossings(2,1)));
   dia_Sub_Wave = dia_Wave(positive_Zero_Crossings(1,1):positive_Zero_Crossings(2,1));
   confidence_Sub_Wave = dia_Confidence(positive_Zero_Crossings(1,1):positive_Zero_Crossings(2,1));
   [~,max_Loc] = max(sub_Wave(:,1));
   [~,min1_Loc] = min(sub_Wave(1:max_Loc,1));
   [~,min2_Loc] = min(sub_Wave(max_Loc:end,1));
   min2_Loc= min2_Loc + max_Loc-1;
   slope_Rise = abs((sub_Wave(max_Loc,1)-sub_Wave(min1_Loc,1))/abs((max_Loc-min1_Loc)));
   slope_Fall = abs((sub_Wave(max_Loc,1)-sub_Wave(min2_Loc,1))/abs((min2_Loc-max_Loc)));
   
   C_Dist_Dia = find_Corr((sub_Wave(:,1) -mean(sub_Wave(:,1))),(dia_Sub_Wave(:,1)- mean(dia_Sub_Wave(:,1))),1,1);
   if( C_Dist_Dia(1,1) >0)
       [distension_Max_Val distension_Max_Point] = max(sub_Wave(:,1));
       [distension_Min_Val, distension_Min_Point] = min(sub_Wave(:,1));
       Del_D = (distension_Max_Val-distension_Min_Val);
       dia_Sub_Wave = corrected_Dia(dia_Sub_Wave,confidence_Sub_Wave, sub_Wave);
       D_d = dia_Sub_Wave(distension_Min_Point,1);
       PW_Wave = PW_Wave(positive_Zero_Crossings(1,1)+1:end);
       DW_Wave = DW_Wave(positive_Zero_Crossings(1,1)+1:end);
       dia_Wave = dia_Wave(positive_Zero_Crossings(1,1)+1:end);
       TS_Wave = TS_Wave(positive_Zero_Crossings(1,1)+1:end);
       dia_Confidence = dia_Confidence(positive_Zero_Crossings(1,1)+1:end);
       dist_Wave = sub_Wave(:,1);
       state_Change=1;
       if(slope_Rise > 2*slope_Fall)
           on_Carotid =1;
       elseif(slope_Fall > slope_Rise)
           on_Carotid =-1;
       else
           on_Carotid =0;
       end 
   else
       sub_Wave=[0];
   end
elseif(neg_Count>1)
   sub_Wave = dist_Wave(negative_Zero_Crossings(1,1):negative_Zero_Crossings(2,1));
   sub_TS_Wave =TS_Wave(negative_Zero_Crossings(1,1):negative_Zero_Crossings(2,1));
   dia_Sub_Wave = dia_Wave(negative_Zero_Crossings(1,1):negative_Zero_Crossings(2,1));
   confidence_Sub_Wave = dia_Confidence(negative_Zero_Crossings(1,1):negative_Zero_Crossings(2,1));
   [~,min_Loc] = min(sub_Wave(:,1));
   [~,max1_Loc] =  max(sub_Wave(1:min_Loc));
   [~,max2_Loc] =  max(sub_Wave(min_Loc:end));
   max2_Loc = max2_Loc+min_Loc-1;
   slope_Fall = abs((sub_Wave(max1_Loc,1)-sub_Wave(min_Loc,1))/abs((max1_Loc-min_Loc)));
   slope_Rise = abs((sub_Wave(min_Loc,1)-sub_Wave(max2_Loc,1))/abs((max2_Loc-min_Loc)));
   
   C_Dist_Dia = find_Corr((sub_Wave(:,1) -mean(sub_Wave(:,1))),(dia_Sub_Wave(:,1)- mean(dia_Sub_Wave(:,1))),1,1);
   if( C_Dist_Dia(1,1) >0)
       [distension_Max_Val distension_Max_Point] = max(sub_Wave);
       [distension_Min_Val, distension_Min_Point] = min(sub_Wave);
       Del_D = (distension_Max_Val-distension_Min_Val);
       dia_Sub_Wave = corrected_Dia(dia_Sub_Wave,confidence_Sub_Wave, sub_Wave);
       D_d = dia_Sub_Wave(distension_Min_Point,1);
       PW_Wave = PW_Wave(negative_Zero_Crossings(1,1)+1:end);
       DW_Wave = DW_Wave(negative_Zero_Crossings(1,1)+1:end);
       dia_Wave = dia_Wave(negative_Zero_Crossings(1,1)+1:end);
       TS_Wave = TS_Wave(negative_Zero_Crossings(1,1)+1:end); 
       dia_Confidence = dia_Confidence(negative_Zero_Crossings(1,1)+1:end);
       dist_Wave = sub_Wave(:,1);
       state_Change=1;
       if(slope_Rise > 2*slope_Fall)
           on_Carotid =1;
       elseif(slope_Fall > slope_Rise)
           on_Carotid =-1;
       else
           on_Carotid =0;
       end
   else
       sub_Wave=[0];
   end
end

if(((TS_Wave(end)-TS_Wave(1)) > 1500)&& TS_on==1)
    i =2;
    while((TS_Wave(end)-TS_Wave(i))>1500)
        i= i+1;
    end
    PW_Wave = PW_Wave(i:end);
    DW_Wave = DW_Wave(i:end);
    dia_Wave = dia_Wave(i:end);
    TS_Wave = TS_Wave (i:end);
    dia_Confidence = dia_Confidence(i:end);
elseif((length(PW_Wave) > 1.5*pulse_Rate)&& (TS_on==0))
    PW_Wave = PW_Wave(end-floor(1.5*pulse_Rate)+1:end);
    DW_Wave = DW_Wave(end-floor(1.5*pulse_Rate)+1:end);
    dia_Wave = dia_Wave(end-floor(1.5*pulse_Rate)+1:end);
    dia_Confidence = dia_Confidence(end-floor(1.5*pulse_Rate)+1:end);
end

function [ state_Change,  on_Carotid, Del_D, D_d, dist_Wave,sub_TS_Wave,ecg_sub_Wave,sub_tstamp_ecg,PAT] = check_Cycles_ECG_Method()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global PW_Wave ;
global DW_Wave;
global dia_Wave;
global TS_Wave;
global points_Per_MM;
global pulse_Rate;
global dia_Confidence;
global ecg_data;
global tstamp_ecg;
global Rpeak_detected;
global tstamp_Rpeak;
global ecg_disp_start;
global prev_Rpeak_ind;
global femoral;

global b_ecg_l;
global a_ecg_l;

on_Carotid = 1;
state_Change=0;
Del_D=0;
D_d=0;
PAT=0;
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

sub_Wave = zeros(2,1);
sub_TS_Wave = zeros(1,1);

ecg_sub_Wave = zeros(1,1);
sub_tstamp_ecg = zeros(1,1);

for (i=prev_Rpeak_ind:(length(tstamp_Rpeak)-1))
    if ((TS_Wave(1)<= tstamp_Rpeak(i) )&& (TS_Wave(end)>= tstamp_Rpeak(i+1)) )
        prev_Rpeak_ind =i;
        j=0;
        while (TS_Wave(end-j)> tstamp_Rpeak(i+1))
            j=j+1;
        end
        k=1;
       while (TS_Wave(k)< tstamp_Rpeak(prev_Rpeak_ind))
            k=k+1;
       end
       sub_Wave = dist_Wave(k:end-j);
       sub_TS_Wave =TS_Wave(k:end-j);
       dia_Sub_Wave = dia_Wave(k:end-j);
       confidence_Sub_Wave = dia_Confidence(k:end-j);
       if (length(sub_TS_Wave)<= 5)
           break
       end
       start_ecg_subW= length(find(tstamp_ecg <(sub_TS_Wave(1)-100)));% max((floor((sub_TS_Wave(1)-tstamp_ecg(1))/2)-50) ,1)  ; % no of previous samples 1000(ms)/ecg sampling rate =2
       end_ecg_subW=  length(find(tstamp_ecg <(sub_TS_Wave(end)+100))); %(floor((sub_TS_Wave(length(sub_TS_Wave))-tstamp_ecg(1))/2) +100) ;
       ecg_sub_Wave = filtfilt(b_ecg_l, a_ecg_l,double( ecg_data));
       if (length(ecg_data)> end_ecg_subW)
            ecg_sub_Wave = ecg_sub_Wave(start_ecg_subW:end_ecg_subW);
            sub_tstamp_ecg = tstamp_ecg(start_ecg_subW:end_ecg_subW);
       else
            ecg_sub_Wave = ecg_sub_Wave(start_ecg_subW:end);
            sub_tstamp_ecg = tstamp_ecg(start_ecg_subW:end);
       end
       
       
       ecg_sub_Wave =ecg_sub_Wave * (max(abs(sub_Wave))/max(abs(ecg_sub_Wave))) ;
       
       [~,max_Loc] = max(sub_Wave(:,1));
       [~,min1_Loc] = min(sub_Wave(1:max_Loc,1));
       [~,min2_Loc] = min(sub_Wave(max_Loc:end,1));
       
       min2_Loc= min2_Loc + max_Loc-1;
       slope_Rise = abs((sub_Wave(max_Loc,1)-sub_Wave(min1_Loc,1))/abs((max_Loc-min1_Loc)));
       slope_Fall = abs((sub_Wave(max_Loc,1)-sub_Wave(min2_Loc,1))/abs((min2_Loc-max_Loc)));
        j=j+1;
        
       dist_Wave = sub_Wave(:,1);
       sub_Wave(:,1)=(sub_Wave(:,1)-(sub_Wave(1,1)));%          
       
       above_ZC = length (find(sub_Wave(:,1)>0)) / length(sub_Wave(:,1));
       length_dist=length(sub_Wave(:,1)); %
       if((slope_Rise > 2*slope_Fall)&& (above_ZC >= 0.5)) %((0.75 *length_dist)<=min2_Loc) &&(max_Loc <=(0.5 * length_dist))&&
           on_Carotid =1;
       elseif (above_ZC <= 0.5) %(((0.75 *length_dist)>min2_Loc) || ((0.5 * length_dist)<=max_Loc)) && 
%            on_Carotid =-1;
           [~,min_Loc] = min(sub_Wave(:,1));
           [~,max1_Loc] =  max(sub_Wave(1:min_Loc));
           [~,max2_Loc] =  max(sub_Wave(min_Loc:end));
           max2_Loc = max2_Loc+min_Loc-1;
           slope_Fall = abs((sub_Wave(max1_Loc,1)-sub_Wave(min_Loc,1))/abs((max1_Loc-min_Loc)));
           slope_Rise = abs((sub_Wave(min_Loc,1)-sub_Wave(max2_Loc,1))/abs((max2_Loc-min_Loc)));
           if(slope_Fall > slope_Rise)
               on_Carotid =-1;
           else
               on_Carotid =0;
           end
       else
           on_Carotid =0;
       end 
        
       C_Dist_Dia = find_Corr((sub_Wave(:,1) -mean(sub_Wave(:,1))),(dia_Sub_Wave(:,1)- mean(dia_Sub_Wave(:,1))),1,1);
       if( (C_Dist_Dia(1,1) >0 ||  femoral == 1)&& on_Carotid ==1)
           [distension_Max_Val distension_Max_Point] = max(sub_Wave(:,1));
           [distension_Min_Val, distension_Min_Point] = min(sub_Wave(:,1));
           Del_D = (distension_Max_Val-distension_Min_Val);
           dia_Sub_Wave = corrected_Dia(dia_Sub_Wave,confidence_Sub_Wave, sub_Wave);
           D_d = dia_Sub_Wave(distension_Min_Point,1);
           [sub_TS_Wave_inter,sub_Wave_inter,N_inter] = fftinterpol(sub_TS_Wave,sub_Wave,1);
           [~,max_Loc_int] = max(sub_Wave_inter(:,1));
           [~,min1_Loc_int] = min(sub_Wave_inter(1:max_Loc_int,1));
           diff_dist = diff(sub_Wave_inter,2);
           diff_dist = [0;0;(diff_dist / max(diff_dist))];
           diff_dist_rise = diff_dist((min1_Loc_int):(max_Loc_int));
%            diff_dist(3:end+2) = diff_dist(1:end);
%            diff_dist(1:2)= [0 0];
           [peak_Locs,peak_Vals, valley_Locs, valley_Vals, peak_No, valley_No]=peak_detect(double(diff_dist),0.7,length(diff_dist));
           max_2nd_diff= find(peak_Locs(1:peak_No)>= (min1_Loc_int),1,'first');
           if isempty(max_2nd_diff) 
               max_2nd_diff=1
           end
           if (peak_No>=1 && (peak_Locs(max_2nd_diff) < (max_Loc_int)) )
           PAT = sub_TS_Wave_inter(peak_Locs(max_2nd_diff)) - tstamp_Rpeak(prev_Rpeak_ind) ;
           elseif (peak_Locs(max_2nd_diff) >= (max_Loc_int))
               PAT_grTmax = sub_TS_Wave_inter(peak_Locs(1)) - tstamp_Rpeak(prev_Rpeak_ind) 
               PAT = PAT_grTmax;
           else

            [~, ind] =  max(diff_dist_rise);
             PAT_pk0 = sub_TS_Wave_inter(ind+(min1_Loc_int)) - tstamp_Rpeak(prev_Rpeak_ind) 
             PAT = PAT_pk0;
           end
           PW_Wave = PW_Wave(end-j:end);
           DW_Wave = DW_Wave(end-j:end);
           dia_Wave = dia_Wave(end-j:end);
           TS_Wave = TS_Wave(end-j:end);
           dia_Confidence = dia_Confidence(end-j:end);
           state_Change=1;
       elseif (on_Carotid ==-1)
           PW_Wave = PW_Wave(end-j:end);
           DW_Wave = DW_Wave(end-j:end);
           dia_Wave = dia_Wave(end-j:end);
           TS_Wave = TS_Wave(end-j:end);
           dia_Confidence = dia_Confidence(end-j:end);
           state_Change=1;
       elseif (on_Carotid ==0)
           state_Change=1;
       else
           sub_Wave=[0];
       end
        
        break
    end
end
if((TS_Wave(end)-TS_Wave(1)) > 1500)
    i =2;
    while((TS_Wave(end)-TS_Wave(i))>1500)
        i= i+1;
    end
    PW_Wave = PW_Wave(i:end);
    DW_Wave = DW_Wave(i:end);
    dia_Wave = dia_Wave(i:end);
    TS_Wave = TS_Wave (i:end);
    dia_Confidence = dia_Confidence(i:end);
end

   



function [ correct_Dia_Wav ] = corrected_Dia( dia_Wav, dia_Confidence, dist_Wav )
[~, best_Point] = max(dia_Confidence);

dist_Wav(:,1) = dist_Wav(:,1)-dist_Wav(best_Point);
correct_Dia_Wav = dist_Wav(:,1) + dia_Wav(best_Point);


function [beta_Arr, compliance_Arr, distensibility_Arr, best_beta, avg_Dia,avg_del_D,avg_PAT,SD_PAT, beta_Confidence] = calculate_Stiffness(del_D_Arr, D_d_Arr,PAT_Arr, SBP, DBP)
global stop_All;
global running;
global live;
global req_no_cycle ;
no_cycle= length(del_D_Arr);
for i = 1:no_cycle
    compliance_Arr(i,1) = (del_D_Arr(i))/(SBP-DBP);
    distensibility_Arr(i,1) = compliance_Arr(i,1)/D_d_Arr(i);
    beta_Arr(i,1) = log(SBP/DBP)/((del_D_Arr(i))/D_d_Arr(i));
end

best_beta =  median(beta_Arr);
avg_Dia = median(D_d_Arr);
avg_del_D= median(del_D_Arr);
avg_PAT = median(PAT_Arr);
SD_PAT = std(PAT_Arr);

cycle_ratio = no_cycle / (2*req_no_cycle);
beta_Confidence = (cycle_ratio * 100)/std(beta_Arr);
if (no_cycle < (req_no_cycle))
    beta_Confidence = (cycle_ratio)*100;
elseif ((no_cycle >= (2*req_no_cycle)) && running ==1 && live ==1 )
    stop_All=1 ;
    msgbox(strcat('target no of cycle achieved but SD beta is ',num2str(std(beta_Arr))),'Done');
elseif (((beta_Confidence >= 100 )&& (beta_Confidence ~= inf)) && running ==1 && live ==1)
    stop_All=1 ;
    msgbox('target no of cycle achieved with good beta confidence','Done');
elseif (beta_Confidence == inf)
    beta_Confidence = (cycle_ratio);
end




function [ecg_sub_start] = ecg_processing(sub_ecg_data, sub_ecg_ts,~)
global Rpeak_detected;
global tstamp_Rpeak;

Rpeak_detected = 0;

if (length(sub_ecg_data)>=1400)
            sub_ecg_data   = sub_ecg_data(end-1399:end);
            sub_ecg_ts =    sub_ecg_ts(end-1399:end);
            
            
end
if ((length(sub_ecg_data) >=200) && max (abs(sub_ecg_data)) ~=0)
    
       
        
    sub_ecg_data  = sub_ecg_data /(max (abs(sub_ecg_data)));
    [Rpeak]=ecg_Rpeakdet(sub_ecg_data, 0.9);
    pospeak = Rpeak(:, 1);

    if (length(pospeak)>=2)


        n = floor(pospeak(2))-100;
        if ((n)>100)
        Rpeak_detected = 1;
%         sub_ecg_ts(floor(pospeak(2)));
        tstamp_Rpeak(1:end+1) = [tstamp_Rpeak(1:end); sub_ecg_ts(floor(pospeak(2)))];
        sub_ecg_data  = sub_ecg_data(n+1:end);

        end
    end

end
ecg_sub_start = length( sub_ecg_data);


function [Rpeak]=ecg_Rpeakdet(data, delta)
% Detect R peak (QRS) in ECG waveform
%
%
global b_ecg_l;
global a_ecg_l;
global b_ecg_h;
global a_ecg_h;
% zeros_pad=8;
% temp_4zero=zeros(zeros_pad,1);
% temp_data= [temp_4zero ; data; temp_4zero];
if data==0
    Rpeak_mem = [0, 0];
    Rpeak = Rpeak_mem;
    return
end

data  = filtfilt(b_ecg_l, a_ecg_l,double( data));
% data= temp_data(zeros_pad+1:end-zeros_pad);
d = diff(data);
abs_d= abs(d)- mean(d);
norm_d = d /(max (abs_d));

[peak_Locs,~, valley_Locs, ~, peak_No, valley_No]= peak_detect(double(norm_d), delta,length(norm_d));

    peak_count =0; %p = floor p[]
    
for (i = 1:peak_No)
    for (j = 0:1)
        if (valley_No < (i+j))
            break
        end
        if (((valley_Locs(i+j) - peak_Locs(i))<=20) && ((valley_Locs(i+j) - peak_Locs(i))>=0) )
            Rpeakind = floor((valley_Locs(i+j) + peak_Locs(i))/2);
            peak_count = peak_count + 1;
            Rpeak_mem(peak_count, :) = [Rpeakind data(Rpeakind)];
            break
        end
        
    end
end
if (peak_count>7 || peak_count ==0 )
    Rpeak_mem = [0, 0];
    Rpeak = Rpeak_mem;
    return
else
    Rpeak = Rpeak_mem([1:peak_count],:);
end


function [ calibration_Done, cut_Point, calibration_Frame, Array_for_File] = gen_Calibration_Data(current_Frame, sampling_Rate, samples_Per_Frame, total_No_Of_Frames, reset)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global average_Frame;
global total_Frames;
global current_Frame_Num;
global cut_Location;
global ZCR_Frame;
cut_Point =0;
calibration_Frame = zeros(samples_Per_Frame,1);
Array_for_File = zeros(samples_Per_Frame+2,1);
calibration_Done = 0;
ZCR_Win_Size = 50;

if (reset ==1)
   average_Frame = zeros(samples_Per_Frame,1);
   total_Frames = total_No_Of_Frames;
   cut_Location = 0;
   current_Frame_Num = 0;
   ZCR_Frame = zeros(floor(samples_Per_Frame/ZCR_Win_Size),total_Frames);
end

current_Frame_Num = current_Frame_Num+1;

if(current_Frame_Num < total_Frames)
    average_Frame = average_Frame + (current_Frame-mean(current_Frame))/total_Frames;
else
    cut_Point = find((current_Frame>=(max(average_Frame)/2)),1,'last');
    calibration_Frame = average_Frame;
    calibration_Done = 1;
    Array_for_File = [sampling_Rate;samples_Per_Frame;cut_Point;calibration_Frame];
end


% --- Executes on mouse press over axes background.
function JV_mmt_Store_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to JV_mmt_Store (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% a=1


% --- Executes during object creation, after setting all properties.
function File_Sel_But_CreateFcn(hObject, eventdata, handles)
% hObject    handle to File_Sel_But (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

function [FileName,FilePath] = uigetfileCD(FileSpec, extension, message) 
  %[FileName,PathName] = uigetfileCD(PathName,'*.us;*.ust;*.lvsngl;*.mat','Select the ARTSENS Record');
  [fPath, fName, fExt] = fileparts(FileSpec);
  backDir = cd;
  cd(fPath);
  [FileName, FilePath] = uigetfile(extension, message);
  if ischar(FileName)
    File = fullfile(FilePath, FileName);
    
  else
    File = 0;
    FileName = strcat(fName,fExt);
    FilePath =strcat(fPath,'\');
  end
  cd(backDir);
return;



function W1_Callback(hObject, eventdata, handles)
% hObject    handle to W1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of W1 as text
%        str2double(get(hObject,'String')) returns contents of W1 as a double


% --- Executes during object creation, after setting all properties.
function W1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to W1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function W2_Callback(hObject, eventdata, handles)
% hObject    handle to W2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of W2 as text
%        str2double(get(hObject,'String')) returns contents of W2 as a double


% --- Executes during object creation, after setting all properties.
function W2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to W2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in track.
function track_Callback(hObject, eventdata, handles)
% hObject    handle to track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of track


% --- Executes on button press in disp_tracking.
function disp_tracking_Callback(hObject, eventdata, handles)
% hObject    handle to disp_tracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of disp_tracking



function [cleared_all]=clear_all_prev_data()

global c_measurement_Count;
global c_consider;
global c_dia;
global c_dist;
global c_time_Stamp;
global c_data_Length;
global c_D_d;
global c_del_D;
global c_PAT;
global c_done;
global c_ECG ;
global c_ECG_TS ;

global c_JV_mmt_Count;
global c_JV_Dist;
global c_JV_Data_Length;
global c_JV_Timestamp;
global c_JV_ECG  ;
global c_JV_ECG_TS ;
        
global f_measurement_Count;
global f_consider;
global f_dia;
global f_dist;
global f_time_Stamp;
global f_data_Length;
global f_D_d;
global f_del_D;
global f_PAT;
global f_done;
global f_ECG ;
global f_ECG_TS ;

global c_HR;
global f_HR;
%%--

c_measurement_Count =0;
c_consider =[];
c_dia=[];
c_dist=[];
c_time_Stamp=[];
c_data_Length=[];
c_D_d=[];
c_del_D=[];
c_PAT=[];
c_done=[];
c_ECG =[];
c_ECG_TS =[];

c_JV_mmt_Count=0;
c_JV_Dist=[];
c_JV_Data_Length=[];
c_JV_Timestamp=[];
c_JV_ECG =[] ;
c_JV_ECG_TS=[] ;
c_HR= [];

f_measurement_Count=0;
f_consider=[];
f_dia=[];
f_dist=[];
f_time_Stamp=[];
f_data_Length=[];
f_D_d=[];
f_del_D=[];
f_PAT=[];
f_done=[];
f_ECG =[];
f_ECG_TS =[];
f_HR=[];
cleared_all =1;



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on edit17 and none of its controls.
function edit17_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on file_name and none of its controls.
function file_name_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to file_name (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
global stop_All;
global samples_Per_Frame;
global live;
global running;
stop_All=1;
pause(1);
if (live == 1) 
    if(running ==1)
    ECG_Link(3,0);
    niScope_Link(3,samples_Per_Frame,1);
    end
end
delete(hObject);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on key press with focus on carotid and none of its controls.
function carotid_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to carotid (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
