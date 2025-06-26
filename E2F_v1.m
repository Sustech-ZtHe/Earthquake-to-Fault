function varargout = E2F_v1(varargin)
% E2F_V1 MATLAB code for E2F_v1.fig
%      E2F_V1, by itself, creates a new E2F_V1 or raises the existing
%      singleton*.
%
%      H = E2F_V1 returns the handle to a new E2F_V1 or the handle to
%      the existing singleton*.
%
%      E2F_V1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in E2F_V1.M with the given input arguments.
%
%      E2F_V1('Property','Value',...) creates a new E2F_V1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before E2F_v1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to E2F_v1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help E2F_v1

% Last Modified by GUIDE v2.5 29-Apr-2025 20:05:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @E2F_v1_OpeningFcn, ...
                   'gui_OutputFcn',  @E2F_v1_OutputFcn, ...
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


% --- Executes just before E2F_v1 is made visible.
function E2F_v1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to E2F_v1 (see VARARGIN)

% Choose default command line output for E2F_v1
handles.output = hObject;

% 初始化状态变量
handles.status = 1;
disp('E2F has been launched');

% 动态绑定 CloseRequestFcn
set(hObject, 'CloseRequestFcn', @myGUI_CloseRequestFcn);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes E2F_v1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = E2F_v1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when user attempts to close the GUI.
function myGUI_CloseRequestFcn(hObject, eventdata, handles)
% 在关闭 GUI 时执行清理操作

% 显示关闭状态
disp('E2F is about to close...');

% 删除 GUI
delete(hObject);

% 清除所有变量和环境
clear all;

% 提示清除完成
disp('All variables have been cleared and E2F has been closed');


%% STEP 0
function input_Callback(hObject, eventdata, handles)
[filename,path]=uigetfile('*','Earthquake Catalog')
global  minPoi hypo evla evlo_km evla_km depth mag la0 lo0 colortemplate mi evlo evdp
hypo =load([path filename]);
if size(hypo,2)>11
    error('Please Check Input File Format!')
end
minPoi=10;
colortemplate=rand(10000,10000);mi=4;
lon=hypo(:,9);
lat=hypo(:,8);
depth=hypo(:,10);
mag=hypo(:,11);
la0=min(lat);lo0=min(lon);
evla = hypo(:,8);evlo = hypo(:,9);evdp = hypo(:,10);
evlo_km = deg2km(distance([la0*ones(length(evla),1),evlo],[la0*ones(length(evla),1),lo0*ones(length(evla),1)]));
evla_km = deg2km(distance([evla,lo0*ones(length(evla),1)],[la0*ones(length(evla),1),lo0*ones(length(evla),1)]));
% Spatial Location
figure('Position', [50, 20, 880, 860]);
ax1=axes('Position', [0.1, 0.4, 0.6, 0.5]);
scatter3(evlo_km,evla_km,depth,power(exp(mag),1/5),'filled','MarkerEdgeColor',[153/255 0 0],'MarkerFaceColor',[255/255 51/255 51/255]);
xlim(ax1,[min(evlo_km) max(evlo_km)]);xlabel('Easting (km)')
ylim(ax1,[min(evla_km) max(evla_km)]);ylabel('Northing (km)')
zlim(ax1,[min(depth) max(depth)]);zlabel(ax1,'Depth (km)')
set(ax1,'ZDir','reverse','FontSize',15,'Color',[0.9 0.9 0.9],'LineWidth',1.5,'box','on','XGrid','off','YGrid','off','ZGrid','off');
view(ax1,0,90)
% Magnitude 
ax2=axes('Position', [0.1, 0.1, 0.6, 0.15]);
histogram(mag,'BinWidth', (max(mag)-min(mag))/40,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0 0 0],'LineWidth',1.5);
xlabel('Magnitude');ylabel('Counts')
set(ax2,'fontsize',15,'LineWidth',1.5,'Color',[0.9 0.9 0.9],'box','on');
% Depth Location
ax3=axes('Position', [0.8, 0.4, 0.15, 0.5]);
histogram(depth,'BinWidth', (max(depth)-min(depth))/40,'LineWidth',1.5);hold on
set(ax3, 'color', [0.9 0.9 0.9], 'fontsize', 15, 'linewidth', 1.5,'YDir', 'reverse','XDir','reverse'); box on;
view(ax3,-90,90)
xlabel('Depth');ylabel('Counts')
% myVar
vars = whos;  
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

%% STEP 1
function run_MODE1_Callback(hObject, eventdata, handles)
global fai_MODE
if (get(hObject,'Value') == get (hObject,'Max'))
    fai_MODE=1;
else
    fai_MODE=2;
end
% myVar
vars = whos;  
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function PBAD_factor_Callback(hObject, eventdata, handles)
global PBAD_Multiple main_fnum
inputStr = get(hObject, 'String');
PBAD_Multiple = str2double(inputStr)
main_fnum=1;
% myVar
vars = whos;  
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function PBAD_factor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% myVar
vars = whos;  
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function run_MODE2_Callback(hObject, eventdata, handles)
global fai_MODE
if (get(hObject,'Value') == get (hObject,'Max'))
    fai_MODE=2;
else
    fai_MODE=1;
end
% myVar
vars = whos;  
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function classify_Callback(hObject, eventdata, handles)
global fai PBAD_Multiple fai_MODE main_fnum A B line_nums_1 f_sort la0 lo0 line_nums_2 f_sort_E2F main_fault angle_tolerance hypo evlo_km evla_km  mag hypo_out_points sita currentPath evdp
lim=[];faiup=[];faidown=[];uplim=[];downlim=[];faia=[];faiu=[];
MEDIAN0=[];MEDIAN1=[];f_sort0=[];f_sort1=[];Values=[];rangedown=[];rangeup=[];
if isempty(fai_MODE)
        error('No MODE Selected')
end
hypo_out_points=[evlo_km,evla_km,evdp,mag,hypo(:,2),hypo(:,3),hypo(:,4),hypo(:,5),hypo(:,6),hypo(:,7)];
% 3D Hough Transform - Write Hough Transform Result to 'fault_try.dat'

currentPath = pwd;
% 判断系统平台
if ispc
    sep = '\';
    exe = '.\hough-3d-lines-master\hough3dlines.exe';
    sp_path = [currentPath, '\hough-3d-lines-master\sp.output.txt'];
    ab_path = [currentPath, '\hough-3d-lines-master\ab.output.txt'];
    dat_path = [currentPath, '\hough-3d-lines-master\data\fault_try.dat'];
else
    sep = '/';
    exe = './hough-3d-lines-master/hough3dlines';
    sp_path = [currentPath, '/hough-3d-lines-master/sp.output.txt'];
    ab_path = [currentPath, '/hough-3d-lines-master/ab.output.txt'];
    dat_path = [currentPath, '/hough-3d-lines-master/data/fault_try.dat'];
end
% 写入 fault_try.dat
fopen(dat_path, 'w');
dlmwrite(dat_path, hypo_out_points, 'precision', 8);
% 调用 hough3dlines
command = [exe, ' ', dat_path];
system(command);
% 读取输出结果
f_sort = load(sp_path);
for i = 1:f_sort(end)
    indices = find(f_sort(:,end) == i);
    line_nums_1(i,1) = indices(1);
    line_nums_1(i,2) = indices(end);
end
disp(['lines: ', num2str(line_nums_1(:,1)')])

ab_output = load(ab_path);

A=ab_output(:,1:3);
B=ab_output(:,4:6);
sita=real(acos(sqrt(B(:,1).^2+B(:,2).^2))*180/pi);
fai=rad2deg(atan2(B(:,1),B(:,2)));

for i=1:main_fnum
    fu=[];fd=[];
    fai_rad_ori=atan2(B(:,1),B(:,2));
    fai_deg_ori=rad2deg(fai_rad_ori);
    fai=fai_deg_ori;
    fai(find(fai<0))=fai(find(fai<0))+180;
    [bin_counts,bin_edges,bin_numbers]=histcounts(fai,0:10:180);
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    if fai_MODE==1
        if isempty(PBAD_Multiple)
            error('No Factor ')
        end
        idd=find(bin_counts==max(bin_counts),1);
        PBAD=median(abs(fai-bin_centers(idd)));
        faiup=bin_centers(idd)+PBAD_Multiple*PBAD;
        faidown=bin_centers(idd)-PBAD_Multiple*PBAD;
        nbu=find(bin_edges<=bin_centers(idd),1,'last');
        nbd=find(bin_centers(idd)<=bin_edges,1);
    elseif fai_MODE==2
        if isempty(angle_tolerance)
            error('No Angular Tolerance')
        end
        if isempty(main_fault)
            error('No Main Strike Direction')
        end
        faiup=main_fault(i)+angle_tolerance;
        faidown=main_fault(i)-angle_tolerance;
        nbu=find(bin_edges<=main_fault(i),1,'last');
        nbd=find(main_fault(i)<=bin_edges,1);
    end
    
    if faiup>180
        faiup=faiup-180;
        fu=1;
    end
    if faidown<0
        faidown=180+faidown;
        fd=1;
    end 
    uplim=find(faiup<=bin_edges,1);
    downlim=find(bin_edges<=faidown,1,'last');
    if fu==1
        rangeup=[rangeup;[0,bin_edges(uplim);bin_edges(nbu),bin_edges(end)]];
    else 
        rangeup=[rangeup;[bin_edges(nbu),bin_edges(uplim)]];
    end     
    if fd==1
        rangedown=[rangedown;[0,bin_edges(nbd);bin_edges(downlim),bin_edges(end)]];
    else 
        rangedown=[rangedown;[bin_edges(downlim),bin_edges(nbd)]];
    end
end
for i=1:size(rangeup,1)
    id1=find(bin_edges<=rangeup(i,1),1,'last');
    id2=find(bin_edges<=rangeup(i,2),1,'last');
    id2=id2-1;
    lim=[lim,[id1:id2]];
end
for i=1:size(rangedown,1)
    id1=find(bin_edges<=rangedown(i,1),1,'last');
    id2=find(bin_edges<=rangedown(i,2),1,'last');
    id2=id2-1;
    lim=[lim,[id1:id2]];
end
lim=sort(lim);
lim=unique(lim);
binCenters = bin_edges(1:end-1) + diff(bin_edges)/2;
for i=1:f_sort(end)
    FAI=rad2deg(atan2(B(i,1),B(i,2)));
    if FAI<0
        FAI=FAI+180;
    end
    LENGTH=line_nums_1(i,2)-line_nums_1(i,1)+1;
    index=line_nums_1(i,1):line_nums_1(i,2);
    f_sort(index,11:13)=rand(1,3).*ones(LENGTH,3);
    isInAnyRange=false;
    if isInAnyRange==false
        for i1=1:size(rangeup,1)
            if FAI<=rangeup(i1,2) & rangeup(i1,1)<=FAI    
                MEDIAN1=[MEDIAN1;[median(f_sort(index,1:3)),f_sort(index(1,1),14),size(f_sort(index,1:3),1)]];
                f_sort1=[f_sort1;f_sort(index,:)];
                faia=[faia;FAI];
                isInAnyRange=true;
                break
            end
        end
    end
    if isInAnyRange==false
        for i2=1:size(rangedown,1)
            if FAI<=rangedown(i2,2) & rangedown(i2,1)<=FAI 
                MEDIAN1=[MEDIAN1;[median(f_sort(index,1:3)),f_sort(index(1,1),14),size(f_sort(index,1:3),1)]];
                f_sort1=[f_sort1;f_sort(index,:)];
                faia=[faia;FAI];
                isInAnyRange=true;
                break
            end
        end
    end
    if isInAnyRange == false
        MEDIAN0=[MEDIAN0;[median(f_sort(index,1:3)),f_sort(index(1),14),size(f_sort(index,1:3),1)]];
        f_sort0=[f_sort0;f_sort(index,:)];
        faiu=[faiu;FAI];
    end 
end
if ~isempty(MEDIAN0) && ~isempty(MEDIAN1)
    [uniqueRows, iM1] = unique(MEDIAN1, 'rows', 'stable');
    MEDIAN1=MEDIAN1(iM1,:);
    [uniqueRows, iM0] = unique(MEDIAN0, 'rows', 'stable');
    MEDIAN0=MEDIAN0(iM0,:);
    [uniqueRows, if1] = unique(f_sort1, 'rows', 'stable');
    f_sort1=f_sort1(if1,:);
    [uniqueRows, if0] = unique(f_sort0, 'rows', 'stable');
    f_sort0=f_sort0(if0,:);
    f_sort_E2F=f_sort;
    for i=1:size(MEDIAN0,1)
        guiyi1=max(pdist2(MEDIAN0(i,1:3),MEDIAN1(:,1:3)));
        guiyi2=max(MEDIAN0(i,5)./MEDIAN1(:,5));
        [d,idxmerge]=min(pdist2(MEDIAN0(i,1:3),MEDIAN1(:,1:3))./guiyi1+ ...
                 (MEDIAN0(i,5)./MEDIAN1(:,5)./guiyi2)');
        % [d,idxmerge]=min(pdist2(MEDIAN0(i,1:3),MEDIAN1(:,1:3))./guiyi1);
        M1_cho=MEDIAN1(idxmerge,4);
        M1_color_idx=find(f_sort_E2F(:,14)==M1_cho);
        rs=unique(f_sort_E2F(M1_color_idx,11:13),'rows');
        M0_color_idx=find(f_sort_E2F(:,14)==MEDIAN0(i,4));
        f_sort_E2F(M0_color_idx,11)=rs(1,1);
        f_sort_E2F(M0_color_idx,12)=rs(1,2);
        f_sort_E2F(M0_color_idx,13)=rs(1,3);
        f_sort_E2F(M0_color_idx,14)=M1_cho;
        MEDIAN1(idxmerge,:)=[];    
        if isempty(MEDIAN1)
            break
        end
    end
    f_sort_E2F= sortrows(f_sort_E2F,14);
    line_nums_2=[];
    line_2=unique(f_sort_E2F(:,end));
    for i=1:length(line_2)
       indices = find(f_sort_E2F(:,end) == line_2(i));
       line_nums_2(i,1)=indices(1);
       line_nums_2(i,2)=indices(end);
       line_nums_2(i,3)=line_2(i);
    end
end
if isempty(MEDIAN0) || isempty(MEDIAN1)
    [uniqueRows, iM1] = unique(MEDIAN1, 'rows', 'stable');
    MEDIAN1=MEDIAN1(iM1,:);
    [uniqueRows, iM0] = unique(MEDIAN0, 'rows', 'stable');
    MEDIAN0=MEDIAN0(iM0,:);
    [uniqueRows, if1] = unique(f_sort1, 'rows', 'stable');
    f_sort1=f_sort1(if1,:);
    [uniqueRows, if0] = unique(f_sort0, 'rows', 'stable');
    f_sort0=f_sort0(if0,:);
    f_sort_E2F=f_sort;
    f_sort_E2F= sortrows(f_sort_E2F,14);
    line_nums_2=[];
    line_2=unique(f_sort_E2F(:,end));
    for i=1:length(line_2)
       indices = find(f_sort_E2F(:,end) == line_2(i));
       line_nums_2(i,1)=indices(1);
       line_nums_2(i,2)=indices(end);
       line_nums_2(i,3)=line_2(i);
    end
end
figure('Position', [50, 20, 1580, 960]);
sg=sgtitle('Lineament Identification and Classification');
set(sg,'fontsize',20,'fontweight','bold');
subplot(4,3,[1 4 7])
f_sort1=KMtoDEG(f_sort1,la0,lo0);
Scatter3_E2F(hypo,f_sort1,['Accepted'],1);
subplot(4,3,[2 5 8])
f_sort0=KMtoDEG(f_sort0,la0,lo0);
Scatter3_E2F(hypo,f_sort0,['Unaccepted'],1);
if ~isempty(MEDIAN0)
    MEDIAN0=sortrows(MEDIAN0,-5);
end
f_sort_deg=KMtoDEG(f_sort_E2F,la0,lo0);
subplot(4,3,[3 6 9])
Scatter3_E2F(hypo,f_sort_deg,['Integrated'],1);
disp('Hough Finished')
% Angle statistics
subplot(4,3,10)
POLAR_gram(faia,0)
title('Strike of "Accepted"')
subplot(4,3,11)
POLAR_gram(faiu,1)
title('Strike of "Unaccepted"')
subplot(4,3,12)
POLAR_gram(faia,2,faiu)
title('Strike of all')
% myVar
vars = whos;  
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function main_strike_direction_Callback(hObject, eventdata, handles)
global main_fault main_fnum
inputStr = get(hObject, 'String'); 
main_fault = str2double(split(inputStr,',')) 
main_fnum=length(main_fault);
guidata(hObject,handles)
% myVar
vars = whos; 
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function main_strike_direction_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function PBAD_factor_text_ButtonDownFcn(hObject, eventdata, handles)

%% STEP 2
function C_manual_edit_Callback(hObject, eventdata, handles)
global C Cmode
inputStr = get(hObject, 'String'); 
C = [str2double(split(inputStr,','))]'
Cmode=2;
% myVar
vars = whos;  
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function C_manual_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% myVar
vars = whos;  
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function Segment_cluster_Callback(hObject, eventdata, handles)
global C f_sort_E2F line_nums_2 hypo minPoi la0 lo0 colortemplate Kused FE_deg Cmode
Kused=[];figureloc=[1,2,3,4,9,10,11,12;5,6,7,8,13,14,15,16];
if isempty(C)
    error('No C Value')
end
if Cmode==1
    figure('Position', [50, 20, 1580, 960]);
    sg=sgtitle('Fault Segment Clustering')
    set(sg,'fontsize',20,'fontweight','bold')
    annotation('textbox', [0.42 0.01 0.2 0.05], 'String', 'Longitude', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 20);
    annotation('textbox', [0.1 0.44 0.2 0.2], 'String', 'Latitude', ...
        'EdgeColor', 'none', 'Rotation', 90, 'VerticalAlignment', 'middle', 'FontSize', 20);
end
for l=1:length(C)     
    if Cmode==1
        subplot(4,4,[figureloc(1,l),figureloc(2,l)])
    elseif Cmode==2
        figure('Position', [50+l*50, 20, 680, 580]);
        sg=sgtitle('Fault Segment Clustering')
        set(sg,'fontsize',20,'fontweight','bold')
    end
    n=1;FE_km=[];RJA=[];RJB=[];
    nw=1;RabY=[];c=C(l);
    for i=1:f_sort_E2F(end)
        im=1;
        if ismember(i,line_nums_2(:,3))
            k=find(line_nums_2(:,3)==i);
            ftest4=f_sort_E2F(line_nums_2(k,1):line_nums_2(k,2),:);
            while 1
                Mobs=max(ftest4(:,4));
                if size(ftest4,1)>100
                    [a,b,Mthre]=Mmain(ftest4,hypo);
                    Mw=max(Mobs,Mthre);
                else
                    Mw=Mobs;
                end
                % RJA=[RJA,a];
                % RJB=[RJB,b];
                radius=c*radius_model(Mw);
                [~,s4]=sort(ftest4(:,4),'descend');
                ftest4=ftest4(s4,:);
                idx_5=dbscan(ftest4(:,1:3),radius,minPoi);
                nw=nw+1;
                if all(idx_5 == -1) || isempty(idx_5)
                    % RabY=[RabY;i,a,b,size(ftest4,1),Mobs,Mthre,Mw,0,];
                    break
                end
                % RabY=[RabY;i,a,b,size(ftest4,1),Mobs,Mthre,Mw,1];
                Kused=[Kused;Mw,1,l];
                F3=ftest4(find(idx_5==1),:);
                F3(:,11)=colortemplate(i,3*(n-1)+1);F3(:,12)=colortemplate(i,3*(n-1)+2);F3(:,13)=colortemplate(i,3*(n-1)+3);
                F3(:,14)=n;
                FE_km=[FE_km;F3];
                n=n+1;
                ftest4(idx_5==1,:)=[];
                if isempty(ftest4)
                    break
                end
            end
        end
    end 
    if isempty(FE_km)
        disp('C is too small, please increase it !')
        break
    end
    FE_deg=KMtoDEG(FE_km,la0,lo0);
    Scatter3_E2F(hypo,FE_deg,['C= ',num2str(C(l))],0);
    % figure;histogram(RJA,20);title('a')
    % figure;histogram(RJB,20);title('b')
end
disp('Segment Clustering Finished');
% myVar
vars = whos;  
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function Fault_fitting_Callback(hObject, eventdata, handles)
global C f_sort_E2F line_nums_2 hypo minPoi la0 lo0 colortemplate FaultParameters Magrecord ASPECT3D MwD Cmode
ASPECT3D=[];MwD=[];
Azimuth=[];Elevation=[];FaultParameters=[];Magrecord=[];
figureloc=[1,2,3,4,9,10,11,12;5,6,7,8,13,14,15,16];
if isempty(C)
    error('No C value')
end
if Cmode==1
    figure('Position', [50, 20, 1580, 960]);
    sg=sgtitle('Fault Zone Fitting')
    set(sg,'fontsize',20,'fontweight','bold')
    annotation('textbox', [0.42 0.01 0.2 0.05], 'String', 'Longitude', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 20);
    annotation('textbox', [0.1 0.44 0.2 0.2], 'String', 'Latitude', ...
        'EdgeColor', 'none', 'Rotation', 90, 'VerticalAlignment', 'middle', 'FontSize', 20);
end
for l=1:length(C)
    if Cmode==1
        subplot(4,4,[figureloc(1,l),figureloc(2,l)])
    elseif Cmode==2
       figure('Position', [50+l*50, 20, 680, 580]);
        sg=sgtitle('Fault Fitting')
        set(sg,'fontsize',20,'fontweight','bold')
    end
    [MaxMw,FE_km,FE_deg]=zdHouDBS(C,l,f_sort_E2F,line_nums_2,la0,lo0,minPoi,hypo,colortemplate);
    MwD=[MwD;[MaxMw,C(l).*ones(length(MaxMw),1)]];
    EP_km=[];n=0;ell=[];
    line_e=unique(FE_deg(:,end));
    for i=1:length(line_e)
       line_nums_e(i,1)=find(ismember(FE_deg(:,end),line_e(i,1),'rows'),1);
       line_nums_e(i,2)=find(ismember(FE_deg(:,end),line_e(i,1),'rows'),1,'last')+1;
       line_nums_e(i,3)=line_e(i,1);
    end
    for i=1:length(line_e)
        n=n+1;
        sp_km=FE_km(line_nums_e(i,1):line_nums_e(i,2)-1,:);
        m1 = mean(sp_km(:,1:3));
        [a,b,c,p1_u,p1_d,p2_u,p2_d,p3_u,p3_d,az,dip,aspect3d,ellipsoid_points]=plotcov_3d(m1,sp_km);
        ell=[ell;[ellipsoid_points,ones(size(ellipsoid_points,1),1).*n]];
        Azimuth=[Azimuth;[az,l]];
        Elevation=[Elevation;[dip,l]];
        Magrecord=[Magrecord;[sp_km(:,4),ones(size(sp_km,1),1).*l,ones(size(sp_km,1),1).*n]];
        FaultParameters=[FaultParameters;[median(sp_km(:,1:3)),a,b,c,p1_u,p1_d,p2_u,p2_d,p3_u,p3_d,az,dip,aspect3d,C(l),size(sp_km,1),sp_km(1,14)]];
        % EP_km=[EP_km;ellipsoid_points,sp_km(1,11:13).*ones(size(ellipsoid_points,1),3)];
        ASPECT3D=[ASPECT3D;[l,size(sp_km,1),aspect3d]];
    end
    Scatter3_E2F(hypo,FE_deg,['C= ',num2str(C(1,l))],0);hold on
    ellsurf(ell,la0,lo0); 
end
disp('Fitting Finished')
%% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

% --- Executes on button press in Evaluation.
function Evaluation_Callback(hObject, eventdata, handles)
global ASPECT3D RSMratio C mi FaultParameters MwD Cmode Magrecord hypo
if isempty(FaultParameters)
        error('Run Fit first')
end
RSMratio=[];figureloc=[1,2,3,4,9,10,11,12;5,6,7,8,13,14,15,16];RSM=[];
for j=1:length(C)
    SubRLength=[];SubRWidth=[];Moreal=[];k=[];
    num0=find(FaultParameters(:,28)==C(j)); 
    FR=FaultParameters(num0,:);
    for i=1:size(FR,1)
        dist1=pdist2([FR(i,7:9)],[FR(i,10:12)]);
        dist2=pdist2([FR(i,13:15)],[FR(i,16:18)]);
        SubRLength=[SubRLength;dist1];
        SubRWidth=[SubRWidth;dist2];
    end
    numD=find(MwD(:,2)==C(j));
    num1=find(Magrecord(:,2)==j);
    D=MagD(MwD(numD,1));
    Motheo=3*10^10.*(10^3.*SubRLength).*(10^3.*SubRWidth).*10^3.*SubRLength.*D;
    % Motheo=3*10^11.*(10^3.*SubRLength).*(10^3.*SubRWidth).*SubRLength;
    MR=Magrecord(num1,:);
    for i1=1:MR(end,3)       
        num2=find(MR(:,3)==i1);
        Moreal=[Moreal;sum(10.^(1.5.*MR(num2,1)+9.1))];
    end
    RSMratio=[RSMratio;log10(Moreal)./log10(Motheo),j.*ones(size(FR,1),1),Motheo,Moreal];
    RSM=[RSM;log10(Moreal)./log10(Motheo)];
end
figure('Position', [50, 20, 1580, 960]);
sg=sgtitle('Choose a C-value');
set(sg,'fontsize',20,'fontweight','bold');
EvaLF(ASPECT3D,RSMratio,size(hypo,1),C,(10-mi)/10);
%% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

% --- Executes on button press in Dip.
function Dip_Callback(hObject, eventdata, handles)
global FaultParameters C Cmode
if isempty(FaultParameters)
        error('Run Fit first')
end
figureloc=[1,2,3,4,9,10,11,12;5,6,7,8,13,14,15,16];
if Cmode==1
    figure('Position', [50, 20, 1580, 960]);
    sg1=title('Dip');
    set(sg1,'fontsize',20,'fontweight','bold');
    maxel=[];
    for l=1:length(C)
        eldeg=FaultParameters(FaultParameters(:,28)==C(l),26);
        maxel=[maxel;max(histcounts(eldeg.*pi/180,(0:10:90).*pi/180))];
    end
    for l=1:length(C)
        subplot(4,4,[figureloc(1,l),figureloc(2,l)])
        %% Dip
        Dip=FaultParameters(FaultParameters(:,28)==C(l),26);
        ElevationStat(Dip,max(maxel));
        title(sprintf('C: %s',num2str(C(l))));
    end
elseif Cmode==2
    for l=1:length(C)
        figure('Position', [50+l*50, 20, 680, 580]);
        set(gca,'fontsize',20,'fontweight','bold');
        %% Dip
        Dip=FaultParameters(FaultParameters(:,28)==C(l),26);
        ElevationStat(Dip);
        title(sprintf('C: %s (Dip)',num2str(C(l))));
    end
end
%% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

% --- Executes on button press in Azimuth_and_Dip.
function Azimuth_and_Dip_Callback(hObject, eventdata, handles)
%% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function Optimal_C_Callback(hObject, eventdata, handles)
global Optimal_C
inputStr = get(hObject, 'String'); % 获取 edit text 的字符串内容
Optimal_C = str2double(inputStr) 
%% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

% --- Executes during object creation, after setting all properties.
function Optimal_C_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

% --- Executes on button press in Fault_Model.
function Fault_Model_Callback(hObject, eventdata, handles)
global Optimal_C f_sort_E2F line_nums_2 hypo minPoi la0 lo0 colortemplate Optimal_FaultParameters FE_deg_O
if isempty(line_nums_2)
        error('Run Classify first')
end
figure('Position', [50, 20, 880, 960]); 
Optimal_FaultParameters=[];
[MaxMw,FE_km,FE_deg_O]=zdHouDBS(Optimal_C,1,f_sort_E2F,line_nums_2,la0,lo0,minPoi,hypo,colortemplate);
line_e=unique(FE_deg_O(:,end));
for i=1:length(line_e)
   line_nums_e(i,1)=find(ismember(FE_deg_O(:,end),line_e(i,1),'rows'),1);
   line_nums_e(i,2)=find(ismember(FE_deg_O(:,end),line_e(i,1),'rows'),1,'last')+1;
   line_nums_e(i,3)=line_e(i,1);
end
for i=1:length(line_e)
    sp_km=FE_km(line_nums_e(i,1):line_nums_e(i,2)-1,:);
    m1 = mean(sp_km(:,1:3));
    [a,b,c,p1_u,p1_d,p2_u,p2_d,p3_u,p3_d,az,dip,aspect3d,ellipsoid_points]=plotcov_3d(m1,sp_km);
    Optimal_FaultParameters=[Optimal_FaultParameters;[median(sp_km(:,1:3)),a,b,c,p1_u,p1_d,p2_u,p2_d,p3_u,p3_d,az,dip,aspect3d,Optimal_C,size(sp_km,1),sp_km(1,14)]];
end
sgtitle(sprintf(['Fault Structure Modeling\n',sprintf('Optimal C-value = %.4f',Optimal_C)]));
set(gca,'fontsize',20,'fontweight','bold');
Rectangle_points=SUBfault_con(Optimal_C,f_sort_E2F,line_nums_2,hypo,Optimal_FaultParameters,minPoi,colortemplate);
%% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

% --- Executes on button press in out_file.
function out_file_Callback(hObject, eventdata, handles)
global currentPath FE_deg_O Optimal_FaultParameters la0 lo0
header = {
    '#1 Event Easting (km)',...
    '#2 Event Northing (km)',...
    '#3 Depth',...
    '#4 Magnitude',...
    '#5-10 Year Month Day Hour Minute Second',...
    '#11-13 R G B',...
    '#14 Fault ID',...
};
fileID = fopen([currentPath,'/OutputFile/Fault_Segment_Clusterd.txt'],'w');
for i=1:length(header)
    fprintf(fileID, '%s\n', header{i}); 
end
fclose(fileID);

FE_km_O=FE_deg_O;
evlo=FE_km_O(:, 1);evla=FE_km_O(:, 2);
FE_km_O(:, 1) = deg2km(distance([la0 * ones(length(evla), 1), evlo], [la0 * ones(length(evla), 1), lo0 * ones(length(evla), 1)]));
FE_km_O(:, 2) = deg2km(distance([evla, lo0 * ones(length(evla), 1)], [la0 * ones(length(evla), 1), lo0 * ones(length(evla), 1)]));

dlmwrite([currentPath,'/OutputFile/Fault_Segment_Clusterd.txt'], FE_km_O,'delimiter', ' ','precision', 6, '-append'); % 将A写入指定文件，元素间默认以逗号分隔
disp(['Has Written in ',[currentPath,'/OutputFile/Fault_Segment_Clusterd.txt']])

header = {
    '#1 Fault ID',...
    '#2-4 Centroid X (km) Centroid Y (km) Centroid Z (km)',...
    '#5-7 Major Axis Upper Endpoint X (km) Major Axis Upper Endpoint Y (km) Major Axis Upper Endpoint Z (km)',...
    '#8-10 Major Axis Lower Endpoint X (km) Major Axis Lower Endpoint Y (km) Major Axis Lower Endpoint Z (km)',...
    '#11-13 Intermediate Axis Upper Endpoint X (km) Intermediate Axis Upper Endpoint Y (km) Intermediate Axis Upper Endpoint Z (km)',...
    '#14-16 Intermediate Axis Lower Endpoint X (km) Intermediate Axis Lower Endpoint Y (km) Intermediate Axis Lower Endpoint Z (km)',...
    '#17-19 Short Axis Upper Endpoint X (km) Short Axis Upper Endpoint Y (km) Short Axis Upper Endpoint Z (km)',...
    '#20-22 Short Axis Lower Endpoint X (km) Short Axis Lower Endpoint Y (km) Short Axis Lower Endpoint Z (km)',...
    '#23 Strike (deg)',...
    '#24 Dip Angle (deg)',...
    '#25 Optimal C Value',...
    '#26 Number of Events constituting the fault'
};
fileID = fopen([currentPath,'/OutputFile/Fault_Segment_Modeling.txt'], 'w');
for i=1:length(header)
    fprintf(fileID, '%s\n', header{i}); 
end
fclose(fileID);
Fault_Segment_Modeling=[Optimal_FaultParameters(:,30),Optimal_FaultParameters(:,1:3),Optimal_FaultParameters(:,7:24),Optimal_FaultParameters(:,25:26),Optimal_FaultParameters(:,28:29)];
dlmwrite([currentPath,'/OutputFile/Fault_Segment_Modeling.txt'], Fault_Segment_Modeling,'delimiter', ' ','precision', 6, '-append'); 
disp(['Has Written in ',[currentPath,'/OutputFile/Fault_Segment_Modeling.txt']])
%% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

% --- Executes during object creation, after setting all properties.
function main_fault_orientation_text_CreateFcn(hObject, eventdata, handles)

% --- Executes on key press with focus on PBAD_factor and none of its controls.
function PBAD_factor_KeyPressFcn(hObject, eventdata, handles)

% --- Executes on key press with focus on MAD_m2_edit and none of its controls.
function MAD_m2_edit_KeyPressFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function STEP0_text_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function STEP1_text_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function STEP2_text_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function STEPEND_text_CreateFcn(hObject, eventdata, handles)

function angle_tolerance_Callback(hObject, eventdata, handles)
global angle_tolerance
inputStr = get(hObject, 'String'); % 获取 edit text 的字符串内容
angle_tolerance = str2double(inputStr) 
guidata(hObject,handles)
%% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

% --- Executes during object creation, after setting all properties.
function angle_tolerance_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% myVar
vars = whos;  
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function Preferred_Orientation_Callback(hObject, eventdata, handles)
global Optimal_C Optimal_FaultParameters
if isempty(Optimal_FaultParameters)
        error('Run Fault Model first')
end
Azimuth=Optimal_FaultParameters(Optimal_FaultParameters(:,28)==Optimal_C,25);
Azimuth(find(Azimuth<0))=Azimuth(find(Azimuth<0))+180;
[counts, edges, number]=histcounts(Azimuth,'BinEdges',0:10:180)
binCenters = edges(1:end-1) + diff(edges)/2
FRend=Optimal_FaultParameters(:,29);
Aznum=[];
for j=1:size(counts,2)
    Aznum=[Aznum,sum(FRend(find(number==j)))]
end
figure 
title('Fault Orientation');
set(gca,'color',[0.9 0.9 0.9],'fontsize',12,'LineWidth',1.5,'box','on')
for i=1:length(counts)
    hold on
    bar(binCenters(i), Aznum(i), 10, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5);
end
[~,id]=maxk(Aznum,1);
bar(binCenters(id), Aznum(id), 10, 'FaceColor',[1 0 0],'EdgeColor',[0 0 0],'LineWidth',1.5);
xticks(0:10:180)
xlabel('Azimuth')
ylabel('Event counts')
% myVar
vars = whos; 
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over run_MODE1.
function run_MODE1_ButtonDownFcn(hObject, eventdata, handles)

% --- Executes on button press in Strike.
function Strike_Callback(hObject, eventdata, handles)
global FaultParameters C Cmode
if isempty(FaultParameters)
        error('Run Fit first')
end
figureloc=[1,2,3,4,9,10,11,12;5,6,7,8,13,14,15,16];
if Cmode==1
    figure('Position', [50, 20, 1580, 960]);
    sg1=sgtitle('Strike');
    set(sg1,'fontsize',20,'fontweight','bold');
    maxaz=[];
    for l=1:length(C)
        azdeg=FaultParameters(FaultParameters(:,28)==C(l),25);
        maxaz=[maxaz;max(histcounts(azdeg.*pi/180,18))];
    end
    for l=1:length(C)
        subplot(4,4,[figureloc(1,l),figureloc(2,l)])
        %% Azimuth
        Azimuth=FaultParameters(FaultParameters(:,28)==C(l),25);
        AzimuthStat([Azimuth,Azimuth-180],max(maxaz));
        title(sprintf('C: %s',num2str(C(1,l))));
        pos = get(gca, 'Position');
    end
elseif Cmode==2
    for l=1:length(C)
        figure('Position', [50+l*50, 20, 680, 580]);
        set(gca,'fontsize',20,'fontweight','bold');
        %% Azimuth
        Azimuth=FaultParameters(FaultParameters(:,28)==C(l),25);
        AzimuthStat([Azimuth,Azimuth-180]);
        title(sprintf('C: %s (Strike)',num2str(C(1,l))));
    end
end
%% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

% --- Executes on button press in RSM.
function RSM_Callback(hObject, eventdata, handles)
global RSMratio C mi FaultParameters MwD Cmode Magrecord
if isempty(FaultParameters)
        error('Run Fit first')
end
RSMratio=[];figureloc=[1,2,3,4,9,10,11,12;5,6,7,8,13,14,15,16];RSM=[];
for j=1:length(C)
    SubRLength=[];SubRWidth=[];Moreal=[];
    id=find(FaultParameters(:,28)==C(j)); 
    FR=FaultParameters(id,:);
    for i=1:size(FR,1)
        dist1=pdist2([FR(i,7:9)],[FR(i,10:12)]);
        dist2=pdist2([FR(i,13:15)],[FR(i,16:18)]);
        SubRLength=[SubRLength;dist1];
        SubRWidth=[SubRWidth;dist2];
    end
    idD=find(MwD(:,2)==C(j));
    idmag=find(Magrecord(:,2)==j);
    D=MagD(MwD(idD,1));
    Motheo=3*10^10.*(10^3.*SubRLength).*(10^3.*SubRWidth).*10^3.*SubRLength.*D;
    mag=Magrecord(idmag,:);
    for i1=1:mag(end,3)       
        num2=find(mag(:,3)==i1);
        Moreal=[Moreal;sum(10.^(1.5.*mag(num2,1)+9.1))];
    end
    RSMratio=[RSMratio;log10(Moreal)./log10(Motheo),j.*ones(size(FR,1),1),Motheo,Moreal];
    RSM=[RSM;log10(Moreal)./log10(Motheo)];
end
% log10(RSMratio(:,4))./log10(RSMratio(:,3))
% log(RSMratio(:,4))./log(RSMratio(:,3))
if Cmode==1
    figure('Position', [50, 20, 1580, 960]);
    sg1=sgtitle('RSM');
    set(sg1,'fontsize',20,'fontweight','bold');
    annotation('textbox', [0.4 0.01 0.2 0.05], 'String', 'Ellipsoid based Mo', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 20);
    annotation('textbox', [0.12 0.4 0.2 0.2], 'String', 'Cumulative Mo of Events', ...
        'EdgeColor', 'none', 'Rotation', 90, 'VerticalAlignment', 'middle', 'FontSize', 20);
    for l=1:length(C)
        subplot(4,4,[figureloc(1,l),figureloc(2,l)])
        %% RSM
        idxe=find(RSMratio(:,2)==l);
        loglog(RSMratio(idxe,3),RSMratio(idxe,3),'LineWidth',2,'LineStyle','-','Color',[0 0 0]);hold on
        d1=min(RSMratio(idxe,3));
        d2=max(RSMratio(idxe,3));
        dx=(d2-d1)/20;
        loglog(d1:dx:d2, 10^0.6 .* (d1:dx:d2).^0.6, 'LineWidth',1.5,'LineStyle',':','Color',[0 0 0]);
        loglog(d1:dx:d2, 10^1.4 .* (d1:dx:d2).^1.4, 'LineWidth',1.5,'LineStyle','--','Color',[0 0 0]);
        loglog(RSMratio(idxe,3) ,RSMratio(idxe,4),'+','MarkerSize',8,'LineWidth',1.5,'MarkerEdgeColor',[161 217 156]./255);
        set(gca,'FontSize',15,'LineWidth',1);
        legend('RSM=1',['RSM=',num2str(1-mi/10)],['RSM=',num2str(1+mi/10)],'location','southeast','box','off');
        pos = get(gca, 'Position');
        annotation('textbox', [pos(1)+pos(3)*0.08, pos(2)+pos(4)*0.12, pos(3), 0.8*pos(4)], 'LineStyle','none',...
            'String', ['C= ',num2str(C(1,l))], ...
            'FontAngle','normal','FontSize', 15, 'Color', 'black');
    end
elseif Cmode==2
    mi
    for l=1:length(C)
        figure('Position', [50+l*50, 20, 680, 580]);
        sg1=sgtitle('RSM ');
        set(sg1,'fontsize',20,'fontweight','bold');
        %% RSM
        idxe=find(RSMratio(:,2)==l);
        d1=min(RSMratio(idxe,3));
        d2=max(RSMratio(idxe,3));
        dx=(d2-d1)/20;
        loglog(RSMratio(idxe,3),RSMratio(idxe,3),'LineWidth',2,'LineStyle','-','Color',[0 0 0]);hold on
        loglog(d1:dx:d2, 10^0.6 .* (d1:dx:d2).^0.6, 'LineWidth',1.5,'LineStyle',':','Color',[0 0 0]);
        loglog(d1:dx:d2, 10^1.4 .* (d1:dx:d2).^1.4, 'LineWidth',1.5,'LineStyle','--','Color',[0 0 0]);
        loglog(RSMratio(idxe,3) ,RSMratio(idxe,4),'+','MarkerSize',8,'LineWidth',1.5,'MarkerEdgeColor',[161 217 156]./255);
        set(gca,'FontSize',10,'LineWidth',1);
        xlabel('Ellipsoid-based Mo');
        ylabel('Cumulative Mo of Events');
        legend('RSM=1',['RSM=',num2str(1-mi/10)],['RSM=',num2str(1+mi/10)],'location','southeast','box','off');
        pos = get(gca, 'Position');
        annotation('textbox', [pos(1)+pos(3)*0.08, pos(2)+pos(4)*0.12, pos(3), 0.8*pos(4)], 'LineStyle','none',...
            'String', ['C= ',num2str(C(1,l))], ...
            'FontAngle','normal','FontSize', 15, 'Color', 'black'); 
        % figure;histogram(result,10)
    end
end
%% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function BACKGROUMD_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in Aspect.
function Aspect_Callback(hObject, eventdata, handles)
global ASPECT3D C Cmode
if isempty(ASPECT3D)
        error('Run Fit first')
end
figureloc=[1,2,3,4,9,10,11,12;5,6,7,8,13,14,15,16];
if Cmode==1
    figure('Position', [50, 20, 1580, 960]);
    sg1=sgtitle('3-D Aspect Ratio');
    set(sg1,'fontsize',20,'fontweight','bold');
    maxylim=[];
    for l=1:length(C)
        
        id=find(ASPECT3D(:,1)==l);
        maxylim=[maxylim;max(histcounts(ASPECT3D(id,3),'Binedges',0:0.1:1))];
    end
    for l=1:length(C)
        subplot(4,4,[figureloc(1,l),figureloc(2,l)])
        id=find(ASPECT3D(:,1)==l);
        histogram(ASPECT3D(id,3),'BinEdges',0:0.1:1,'LineWidth',2,'FaceColor',[203 238 249]./255,'EdgeColor',[66 146 197]./255)
        ylim([0 max(maxylim)+5])
        set(gca,'FontSize',15,'LineWidth',1.5,'box','on');
        pos = get(gca, 'Position');
        annotation('textbox', [pos(1)+pos(3)*0.08, pos(2)+pos(4)*0.12, pos(3), 0.8*pos(4)], 'LineStyle','none',...
            'String', ['C= ',num2str(C(1,l))], ...
            'FontAngle','normal','FontSize', 15, 'Color', 'black');
    end
elseif Cmode==2
    for l=1:length(C)
        figure('Position', [50+l*50, 20, 680, 580]);
        sg1=sgtitle('3-D Aspect Ratio');
        set(sg1,'fontsize',20,'fontweight','bold');
        id=find(ASPECT3D(:,1)==l);
        histogram(ASPECT3D(id,3),'BinEdges',0:0.1:1,'LineWidth',2,'FaceColor',[203 238 249]./255,'EdgeColor',[66 146 197]./255)   
        set(gca,'FontSize',15,'LineWidth',1.5,'box','on');
        pos = get(gca, 'Position');
        annotation('textbox', [pos(1)+pos(3)*0.08, pos(2)+pos(4)*0.12, pos(3), 0.8*pos(4)], 'LineStyle','none',...
            'String', ['C= ',num2str(C(1,l))], ...
            'FontAngle','normal','FontSize', 15, 'Color', 'black');
    end
end
%% myVar
vars = whos;  % 获取当前函数中的所有变量信息
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end


% --- Executes on button press in degree3.
function degree3_Callback(hObject, eventdata, handles)
global C Cmode
if (get(hObject,'Value') == get (hObject,'Max'))
    C=[50,60,70,80,90,100,110,120]
    Cmode=1;
end
% myVar
vars = whos;  
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end


% --- Executes on button press in degree1.
function degree1_Callback(hObject, eventdata, handles)
global C Cmode
if (get(hObject,'Value') == get (hObject,'Max'))
    C=[1,1.5,2,2.5,3,3.5,4,5]
    Cmode=1;
end
% myVar
vars = whos;  
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

% --- Executes on button press in degree2.
function degree2_Callback(hObject, eventdata, handles)
global C Cmode
if (get(hObject,'Value') == get (hObject,'Max'))
    C=[5,10,15,20,25,30,35,40]
    Cmode=1;
end
% myVar
vars = whos;  
for i = 1:length(vars)
    assignin('base', vars(i).name, eval(vars(i).name));
end

function edit15_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Manual_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on key press with focus on main_strike_direction and none of its controls.
function main_strike_direction_KeyPressFcn(hObject, eventdata, handles)
