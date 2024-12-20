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

% Last Modified by GUIDE v2.5 20-Dec-2024 21:50:08

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


% --- Executes on button press in catalog_file.
function catalog_file_Callback(hObject, eventdata, handles)
[filename,path]=uigetfile('*','Earthquake Catalog')
global R s minPoi Mc hypo evla evlo_km evla_km depth mag hypo_out_points fai sita A B line_nums_1 f_sort la0 lo0 hypotout NEvents currentPath colortemplate abnum
ax1= handles.axes1;
ax2= handles.axes2;
[R,C,s,minPoi,Mc,hypo]=Detect_Region([path filename]);colortemplate=rand(999,999);
lon=hypo(:,9);
lat=hypo(:,8);
depth=hypo(:,10);
mag=hypo(:,11);
hypotout=hypo;
fdata = hypo;
la0=min(lat);lo0=min(lon);
evla = fdata(:,8);evlo = fdata(:,9);evdp = fdata(:,10);
evlo_km = deg2km(distance([la0*ones(length(evla),1),evlo],[la0*ones(length(evla),1),lo0*ones(length(evla),1)]));
evla_km = deg2km(distance([evla,lo0*ones(length(evla),1)],[la0*ones(length(evla),1),lo0*ones(length(evla),1)]));
hypo_out_points=[evlo_km,evla_km,evdp,mag,hypo(:,2),hypo(:,3),hypo(:,4),hypo(:,5),hypo(:,6),hypo(:,7)];
% 3D Hough Transform 
currentPath = pwd;
fault_try_path=[currentPath,'/data/fault_try.dat'];
% Write Hough Transform Result to 'fault_try.dat'
fopen(fault_try_path,'w');
dlmwrite(num2str(fault_try_path), hypo_out_points,'precision', 8);
command = './hough3dlines ./data/fault_try.dat ';
system(command);
% % Write 'fault_try.dat' to 'Fault_predict.txt'
% fopen([num2str(currentPath),'/FaultSegment/Fault_predict.txt'],'w');
% command = './hough3dlines ./data/fault_try.dat > ./FaultSegment/Fault_predict.txt';
% system(command);
% Load Earthquake Data ('sp.output.txt')
f_sort=load([num2str(currentPath),'/sp.output.txt']);
NEvents=size(f_sort,1);
for i=1:f_sort(end)
   indices = find(f_sort(:,end) == i);
   line_nums_1(i,1)=indices(1);
   line_nums_1(i,2)=indices(end);
end
disp(['lines: ',num2str(line_nums_1(:,1)')])
ab_output=readtable([num2str(currentPath),'/ab.output.txt']);
A=table2array(ab_output(:,1:3));
B=table2array(ab_output(:,4:6));
abnum=table2array(ab_output(:,7));
sita=real(acos(sqrt(B(:,1).^2+B(:,2).^2))*180/pi);
fai=rad2deg(atan2(B(:,1),B(:,2)));
%% Spatial location of earthquake catalog
scatter3(evlo_km,evla_km,depth,power(exp(mag),1/5),'filled','MarkerEdgeColor',[153/255 0 0],'MarkerFaceColor',[255/255 51/255 51/255],'Parent',ax1);
view(ax1,-30,60);
xlim(ax1,[min(evlo_km) max(evlo_km)]);xlabel(ax1,'Easting (km)')
ylim(ax1,[min(evla_km) max(evla_km)]);ylabel(ax1,'Northing (km)')
zlim(ax1,[min(depth) max(depth)]);zlabel(ax1,'Depth (km)')
set(ax1,'ZDir','reverse','FontSize',15,'Color',[0.9 0.9 0.9],'LineWidth',1.5,'box','on','XGrid','off','YGrid','off','ZGrid','off');
axis(ax1,'equal')
%% Magnitude information
histogram(mag,'FaceColor',[255 153 51]./255,'EdgeColor',[204 102 0]./255,'LineWidth',1.5,'Parent',ax2);
xlabel(ax2,'Magnitude');ylabel(ax2,'Counts')
set(ax2,'fontsize',12,'LineWidth',1.5);

% --- Executes on button press in run_MODE1.
function run_MODE1_Callback(hObject, eventdata, handles)
global fai_MODE
if (get(hObject,'Value') == get (hObject,'Max'))
    fai_MODE=1;
else
    fai_MODE=2;
end


function PBAD_Multiple_edit_Callback(hObject, eventdata, handles)
global PBAD_Multiple main_fnum
inputStr = get(hObject, 'String');
PBAD_Multiple = str2double(inputStr)
main_fnum=1;

% --- Executes during object creation, after setting all properties.
function PBAD_Multiple_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PBAD_Multiple_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MAD_m2_radio.
function MAD_m2_radio_Callback(hObject, eventdata, handles)
global fai_MODE
if (get(hObject,'Value') == get (hObject,'Max'))
    fai_MODE=2;
else
    fai_MODE=1;
end


function MAD_m2_edit_Callback(hObject, eventdata, handles)
global na2
inputStr = get(hObject, 'String')
na2 = str2double(inputStr)


% --- Executes during object creation, after setting all properties.
function MAD_m2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MAD_m2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Hough_classify.
function Hough_classify_Callback(hObject, eventdata, handles)
global fai PBAD_Multiple fai_MODE main_fnum A B line_nums_1 f_sort la0 lo0 hypotout NEvents line_nums_2 f_sort_E2F main_fault abnum angle_tolerance
ax3= handles.axes3;
cla(ax3)
%% Choose fault azimuth model
lim=[];faiup=[];faidown=[];uplim=[];downlim=[];faia=[];faiu=[];
MEDIAN0=[];MEDIAN1=[];f_sort0=[];f_sort1=[];Values=[];rangedown=[];rangeup=[];
if isempty(fai_MODE)
        error('No MODE Selected')
end
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
            error('No PBAD-Multiple ')
        end
        idd=find(bin_counts==max(bin_counts),1);
        bin_centers(idd)
        PBAD=median(abs(fai-bin_centers(idd)));
        faiup=bin_centers(idd)+PBAD_Multiple*PBAD;
        faidown=bin_centers(idd)-PBAD_Multiple*PBAD;
        nbu=find(bin_edges<=bin_centers(idd),1,'last');
        nbd=find(bin_centers(idd)<=bin_edges,1);
    elseif fai_MODE==2
        if isempty(angle_tolerance)
            error('No Angle Tolerance')
        end
        if isempty(main_fault)
            error('No Main Fault Orientation')
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

% rangedown
% rangeup

binCenters = bin_edges(1:end-1) + diff(bin_edges)/2;
hold(ax3, 'on');
for j=1:size(bin_counts,2)
    Values=[Values,sum(abnum(find(bin_numbers==j)))];
end
for i=1:length(bin_counts)
    if ismember(i,lim)
        % bar(binCenters(i), Values(i), 10, 'FaceColor',[0 0.8 1],'EdgeColor',[0 0.2 0.6],'Parent',ax3,'LineWidth',1.5);
        % bar(binCenters(i), bin_counts(i), 10, 'FaceColor',[0 0.8 1],'EdgeColor',[0 0.2 0.6],'Parent',ax3,'LineWidth',1.5);
        rectangle('Position', [bin_edges(i), 0, 10, bin_counts(i)], 'FaceColor', [0 0.8 1], 'EdgeColor', [0 0.2 0.6], 'LineWidth', 1.5,'Parent',ax3);
    else 
        % bar(binCenters(i), Values(i), 10, 'FaceColor',[1 0.2 0],'EdgeColor',[0.6 0 0],'Parent',ax3,'LineWidth',1.5);
        % bar(binCenters(i), bin_counts(i), 10, 'FaceColor',[1 0.2 0],'EdgeColor',[0.6 0 0],'Parent',ax3,'LineWidth',1.5);
        rectangle('Position', [bin_edges(i), 0, 10, bin_counts(i)], 'FaceColor', [1 0.2 0], 'EdgeColor', [0.6 0 0], 'LineWidth', 1.5,'Parent',ax3);
    end
end
set(ax3,'FontSize',12);
xticks(ax3, 0:10:180);
yticks(ax3, 0:round(max(bin_counts)/5):max(bin_counts));
xlim(ax3,[0 180]);
ylim(ax3,[0 max(bin_counts)]);
for i=1:f_sort(end)
    FAI=rad2deg(atan2(B(i,1),B(i,2)));
    if FAI<0
        FAI=FAI+180;
    end
    LENGTH=line_nums_1(i,2)-line_nums_1(i,1)+1;
    index=line_nums_1(i,1):line_nums_1(i,2);
    f_sort(index,11:13)=rand(1,3).*ones(LENGTH,3);
    isInAnyRange=false;
    if isInAnyRange==false;
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
    if isInAnyRange==false;
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
    if isInAnyRange == false;
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
        [d,idxmerge]=min(pdist2(MEDIAN0(i,1:3),MEDIAN1(:,1:3))./guiyi1);
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
figure('Position', [100, 100, 1880, 1560]);
sg=sgtitle('Lineament Identification and Classification');
set(sg,'fontsize',20,'fontweight','bold');
subplot(4,3,[1 4 7])
f_sort1=KMtoDEG(f_sort1,la0,lo0);
Scatter3_E2F(hypotout,f_sort1,NEvents,['Accepted'],1);
subplot(4,3,[2 5 8])
f_sort0=KMtoDEG(f_sort0,la0,lo0);
Scatter3_E2F(hypotout,f_sort0,NEvents,['Unaccepted'],1);
if ~isempty(MEDIAN0)
    MEDIAN0=sortrows(MEDIAN0,-5);
end
f_sort_deg=KMtoDEG(f_sort_E2F,la0,lo0);
subplot(4,3,[3 6 9])
Scatter3_E2F(hypotout,f_sort_deg,NEvents,['Integrated'],1);
disp('Hough Finished')
%% Azimuth statistics
subplot(4,3,10)
POLAR_gram(faia,0)
title('Accepted Azimuth')
subplot(4,3,11)
POLAR_gram(faiu,1)
title('Unaccepted Azimuth')
subplot(4,3,12)
POLAR_gram(faia,2,faiu)
title('Integrated Azimuth')

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
set(hObject,'ZDir','reverse','FontSize',15,'Color',[0.9 0.9 0.9],'LineWidth',1.5);
view(hObject,-30,60);box on;grid off
xlabel(hObject,'Easting (km)');ylabel(hObject,'Northing (km)');zlabel(hObject,'Depth (km)')
xticks(hObject,[]);yticks(hObject,[]);


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
set(hObject,'FontSize',12,'Box','on','LineWidth',1.5);
xlabel(hObject,'Magnitude');ylabel(hObject,'Counts')
xticks(hObject,[]);yticks(hObject,[]);


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
set(hObject,'FontSize',10,'Box','on','LineWidth',1.5);
xlabel(hObject,'Azimuth');ylabel(hObject,'Fault lines')
xticks(hObject,[]);yticks(hObject,[]);


% --- Executes on button press in C_default_1_radio.
function C_default_1_radio_Callback(hObject, eventdata, handles)
global C Cmode
if (get(hObject,'Value') == get (hObject,'Max'))
    C=[1,1.5,2,2.5,3,3.5,4,5]
    Cmode=1;
end

% --- Executes on button press in C_default_2_radio.
function C_default_2_radio_Callback(hObject, eventdata, handles)
global C Cmode
if (get(hObject,'Value') == get (hObject,'Max'))
    C=[5,10,15,20,25,30,35,40]
    Cmode=1;
end

% --- Executes on button press in C_default_3_radio.
function C_default_3_radio_Callback(hObject, eventdata, handles)
global C Cmode
if (get(hObject,'Value') == get (hObject,'Max'))
    C=[50,60,70,80,90,100,110,120]
    Cmode=1;
end


function C_manual_edit_Callback(hObject, eventdata, handles)
global C Cmode
inputStr = get(hObject, 'String'); % 获取 edit text 的字符串内容
C = str2double(split(inputStr,',')) % 将字符串转换为数值
C=C'
Cmode=2;

% --- Executes during object creation, after setting all properties.
function C_manual_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C_manual_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Segment_cluster.
function Segment_cluster_Callback(hObject, eventdata, handles)
global C f_sort_E2F line_nums_2 hypotout Mc minPoi la0 lo0 NEvents colortemplate Kused FE Cmode
Kused=[];figureloc=[1,2,3,4,9,10,11,12;5,6,7,8,13,14,15,16];
if isempty(C)
    error('No C Value')
end
if Cmode==1
    figure('Position', [100, 100, 1880, 1560]);
    sg=sgtitle('Fault Segement Clustering')
    set(sg,'fontsize',20,'fontweight','bold')
end
for l=1:length(C)     
    if Cmode==1
        subplot(4,4,[figureloc(1,l),figureloc(2,l)])
    elseif Cmode==2
       figure('Position', [100+l*50, 50, 680, 580]);
        sg=sgtitle('Segement Clustering')
        set(sg,'fontsize',20,'fontweight','bold')
    end
    n=1;FE_km=[];
    nw=1;RabY=[];c=C(l);
    for i=1:f_sort_E2F(end)
        im=1;
        if ismember(i,line_nums_2(:,3))
            k=find(line_nums_2(:,3)==i);
            ftest4=f_sort_E2F(line_nums_2(k,1):line_nums_2(k,2),:);
            while 1
                MaxMag=max(ftest4(:,4));
                if im==1
                    [a,b,Y,CauMag]=Mmain(ftest4,Mc,hypotout);
                    im=2;
                elseif im==2
                    CauMag=(log10(size(ftest4,1))-a)/b+Mc;
                end
                Mag=max(MaxMag,CauMag);
                Ns=(log10(size(ftest4,1))-a)/b+Mc;
                radius=c*radius_model(Mag);
                [~,s4]=sort(ftest4(:,4),'descend');
                ftest4=ftest4(s4,:);
                idx_5=dbscan(ftest4(:,1:3),radius,minPoi);
                nw=nw+1;
                if all(idx_5 == -1) || isempty(idx_5)
                    RabY=[RabY;i,a,b,mean(Y(:,1)),size(ftest4,1),MaxMag,CauMag,Ns,Mag,0,];
                    break
                end
                RabY=[RabY;i,a,b,mean(Y(:,1)),size(ftest4,1),MaxMag,CauMag,Ns,Mag,1];
                Kused=[Kused;Mag,1,l];
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
    FE=KMtoDEG(FE_km,la0,lo0);
    Scatter3_E2F(hypotout,FE,NEvents,['C= ',num2str(C(l))],0,l);
end
disp('Segment Clustering Finished');



function main_fault_orientation_edit_Callback(hObject, eventdata, handles)
global main_fault main_fnum
inputStr = get(hObject, 'String'); % 获取 edit text 的字符串内容
main_fault = str2double(split(inputStr,',')) 
main_fnum=length(main_fault);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function main_fault_orientation_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to main_fault_orientation_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over PBAD_Multiple_text.
function PBAD_Multiple_text_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to PBAD_Multiple_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Fault_fitting.
function Fault_fitting_Callback(hObject, eventdata, handles)
global C f_sort_E2F line_nums_2 hypotout Mc minPoi la0 lo0 NEvents colortemplate hypo_out_points s FaultpaRameters Magrecord ECC MwD Cmode
number_EPC=[];number_Eve=[];ECC=[];MwD=[];
ratio=[1;1;1];Azimuth=[];Elevation=[];FaultpaRameters=[];Magrecord=[];
figureloc=[1,2,3,4,9,10,11,12;5,6,7,8,13,14,15,16];
if isempty(C)
    error('No C value')
end
if Cmode==1
    figure('Position', [100, 100, 1880, 1560]);
    sg=sgtitle('Fault Fitting by Ellipsoid')
    set(sg,'fontsize',20,'fontweight','bold')
end
for l=1:length(C)
    if Cmode==1
        subplot(4,4,[figureloc(1,l),figureloc(2,l)])
    elseif Cmode==2
       figure('Position', [100+l*50, 50, 680, 580]);
        sg=sgtitle('Fault Fitting')
        set(sg,'fontsize',20,'fontweight','bold')
    end
    [MaxMw,FE_km,FE]=zdHouDBS(C,l,f_sort_E2F,line_nums_2,la0,lo0,minPoi,Mc,hypotout,colortemplate);
    MwD=[MwD;[MaxMw,C(l).*ones(length(MaxMw),1)]];
    EP=[];EP_km=[];n=0;ell=[];
    line_e=unique(FE(:,end));
    for i=1:length(line_e)
       line_nums_e(i,1)=find(ismember(FE(:,end),line_e(i,1),'rows'),1);
       line_nums_e(i,2)=find(ismember(FE(:,end),line_e(i,1),'rows'),1,'last')+1;
       line_nums_e(i,3)=line_e(i,1);
    end
    for i=1:length(line_e)
        sp_km=FE_km(line_nums_e(i,1):line_nums_e(i,2)-1,:);
        w=power(10.^(11.8+1.5*sp_km(:,4))/max(10.^(11.8+1.5*sp_km(:,4))),1/3);
        L=w*ratio(1,1);W=w*ratio(2,1);H=w*ratio(3,1);weight=[L,W,H];
        m1 = mean(sp_km(:,1:3));
        C1 = cov(sp_km(:,1:3));
        loca=median(hypo_out_points(:,1:3));
        [a,b,c,p1_u,p1_d,p2_u,p2_d,p3_u,p3_d,az,dip,linea,Vol,elliptic_points]=plotcov_3d(C1, m1,sp_km,weight,ratio,0,s,loca,la0,lo0);
        n=n+1;
        ell=[ell;[elliptic_points,ones(size(elliptic_points,1),1).*n]];
        Azimuth=[Azimuth;[az,l]];
        Elevation=[Elevation;[dip,l]];
        mass=KMtoDEG(sp_km,la0,lo0);
        Magrecord=[Magrecord;[sp_km(:,4),ones(size(sp_km,1),1).*l,ones(size(sp_km,1),1).*n]];
        FaultpaRameters=[FaultpaRameters;[median(mass(:,1:3)),a,b,c,p1_u,p1_d,p2_u,p2_d,p3_u,p3_d,az,dip,linea,size(sp_km,1)/Vol,C(l),size(sp_km,1),mass(1,14)]];
        EP_km=[EP_km;elliptic_points,sp_km(1,11:13).*ones(size(elliptic_points,1),3)];
        ECC=[ECC;[linea,l,size(sp_km,1)/Vol,length(sp_km(:,1))]];
    end
    Scatter3_E2F(hypotout,FE,NEvents,['C= ',num2str(C(1,l))],0,l);hold on
    ellsurf(ell,la0,lo0);
    EP=KMtoDEG(EP_km,la0,lo0);  
    number_EPC=[number_EPC;length(unique(EP(:,4:6),'rows'))];
    number_Eve=[number_Eve;length(FE(:,1))];
end
FaultpaRameters(FaultpaRameters(:,28)>1,28)=1;
ECC(ECC(:,3)>1,3)=1;
ratio_Eve=number_Eve./NEvents;
disp('Fitting Finished')


% --- Executes on button press in Evaluation.
function Evaluation_Callback(hObject, eventdata, handles)
global ECC Energyratio NEvents C mi FaultpaRameters MwD Cmode Magrecord
Energyratio=[];figureloc=[1,2,3,4,9,10,11,12;5,6,7,8,13,14,15,16];
for j=1:length(C)
    SubRLength=[];SubRWidth=[];Moreal=[];k=[];
    num0=find(FaultpaRameters(:,29)==C(j));
    FR=FaultpaRameters(num0,:);
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
    Energyratio=[Energyratio;log10(Motheo)./log10(Moreal),j.*ones(size(FR,1),1),Motheo,Moreal];
end
figure('Position', [100, 100, 1880, 1560]);
sg=sgtitle('Fault Evaluation');mi=4;
set(sg,'fontsize',20,'fontweight','bold');
EvaLF(ECC,Energyratio,NEvents,C,(10-mi)/10,2);
if Cmode==1
    figure('Position', [100, 100, 1880, 1560]);
    sg1=sgtitle('RSM');
    set(sg1,'fontsize',20,'fontweight','bold');
    for l=1:length(C)
        subplot(4,4,[figureloc(1,l),figureloc(2,l)])
        %% RSM
        idxe=find(Energyratio(:,2)==l);
        loglog(Energyratio(idxe,4),Energyratio(idxe,4),'LineWidth',2,'LineStyle','-','Color',[0 0 0]);hold on
        d1=min(Energyratio(idxe,4));
        d2=max(Energyratio(idxe,4));
        dx=(d2-d1)/20;
        loglog(d1:dx:d2,10^mi.*(d1:dx:d2),'LineWidth',1.5,'LineStyle','--','Color',[0 0 0]);
        loglog(d1:dx:d2,10^-mi.*(d1:dx:d2),'LineWidth',1.5,'LineStyle',':','Color',[0 0 0]);
        set(gca,'FontSize',15,'LineWidth',1);
        loglog(Energyratio(idxe,4) ,Energyratio(idxe,3),'+','MarkerSize',8,'LineWidth',1.5);
        legend('RSM=1',['RSM=',num2str(1+mi/10)],['RSM=',num2str(1-mi/10)],'location','southeast','box','off');
        pos = get(gca, 'Position');
        annotation('textbox', [pos(1)+pos(3)*0.08, pos(2)+pos(4)*0.12, pos(3), 0.8*pos(4)], 'LineStyle','none',...
            'String', ['C= ',num2str(C(1,l))], ...
            'FontAngle','normal','FontSize', 15, 'Color', 'black');
        % ylim([10^2 10^20])
        % xlim([10^7 10^15])

        % ylim([10^6 10^20])
        % xlim([10^10 10^17])
        % 
        % ylim([10^9 10^25])
        % xlim([10^13 10^20])
        % yticks([10^4 10^12 10^20])
        % xticks([10^8 10^11 10^14])
    end
elseif Cmode==2
    for l=1:length(C)
        figure('Position', [100+l*50, 50, 680, 580]);
        sg1=sgtitle('RSM ');
        set(sg1,'fontsize',20,'fontweight','bold');
        %% RSM
        idxe=find(Energyratio(:,2)==l);
        loglog(Energyratio(idxe,4),Energyratio(idxe,4),'LineWidth',2,'LineStyle','-','Color',[0 0 0]);hold on
        d1=min(Energyratio(idxe,4));
        d2=max(Energyratio(idxe,4));
        dx=(d2-d1)/20;
        loglog(d1:dx:d2,10^mi.*(d1:dx:d2),'LineWidth',1.5,'LineStyle','--','Color',[0 0 0]);
        loglog(d1:dx:d2,10^-mi.*(d1:dx:d2),'LineWidth',1.5,'LineStyle',':','Color',[0 0 0]);
        set(gca,'FontSize',10,'LineWidth',1);
        loglog(Energyratio(idxe,4) ,Energyratio(idxe,3),'+','MarkerSize',8,'LineWidth',1.5);
        legend('RSM=1',['RSM=',num2str(1+mi/10)],['RSM=',num2str(1-mi/10)],'location','southeast','box','off');
        pos = get(gca, 'Position');
        annotation('textbox', [pos(1)+pos(3)*0.08, pos(2)+pos(4)*0.12, pos(3), 0.8*pos(4)], 'LineStyle','none',...
            'String', ['C= ',num2str(C(1,l))], ...
            'FontAngle','normal','FontSize', 15, 'Color', 'black');
    end
end


% --- Executes on button press in RSM_Azimuth_Dip.
function RSM_Azimuth_Dip_Callback(hObject, eventdata, handles)
global FaultpaRameters C Cmode
figureloc=[1,2,3,4,9,10,11,12;5,6,7,8,13,14,15,16];
if Cmode==1
    figure('Position', [100, 100, 1880, 1560]);
    sg1=sgtitle('Dip');
    set(sg1,'fontsize',20,'fontweight','bold');
    for l=1:length(C)
        subplot(4,4,[figureloc(1,l),figureloc(2,l)])
        %% Dip
        Dip=FaultpaRameters(FaultpaRameters(:,29)==C(l),26);
        ElevationStat(Dip);
        title(sprintf('C: %s',num2str(C(l))));
        pos = get(gca, 'Position');
    end
    figure('Position', [100, 100, 1880, 1560]);
    sg1=sgtitle('Azimuth');
    set(sg1,'fontsize',20,'fontweight','bold');
    for l=1:length(C)
        subplot(4,4,[figureloc(1,l),figureloc(2,l)])
        %% Azimuth
        Azimuth=FaultpaRameters(FaultpaRameters(:,29)==C(l),25);
        AzimuthStat([Azimuth,Azimuth-180]);
        title(sprintf('C: %s',num2str(C(1,l))));
        pos = get(gca, 'Position');
    end
elseif Cmode==2
    for l=1:length(C)
        figure('Position', [100+l*50, 50, 680, 580]);
        sg1=sgtitle('Dip and Azimuth');
        set(sg1,'fontsize',20,'fontweight','bold');
        %% Dip
        subplot(1,2,1)
        Dip=FaultpaRameters(FaultpaRameters(:,29)==C(l),26);
        ElevationStat(Dip);
        title(sprintf('C: %s (Dip)',num2str(C(l))));
        %% Azimuth
        subplot(1,2,2)
        Azimuth=FaultpaRameters(FaultpaRameters(:,29)==C(l),25);
        AzimuthStat([Azimuth,Azimuth-180]);
        title(sprintf('C: %s (Azimuth)',num2str(C(1,l))));
    end
end


% --- Executes on button press in Azimuth_and_Dip.
function Azimuth_and_Dip_Callback(hObject, eventdata, handles)


function Optimal_C_edit_Callback(hObject, eventdata, handles)
global Optimal_C
inputStr = get(hObject, 'String'); % 获取 edit text 的字符串内容
Optimal_C = str2double(inputStr) 

% --- Executes during object creation, after setting all properties.
function Optimal_C_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Optimal_C_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Fault_Model.
function Fault_Model_Callback(hObject, eventdata, handles)
global Optimal_C f_sort_E2F line_nums_2 hypotout FaultpaRameters Mc minPoi la0 lo0 colortemplate
figure('Position', [100, 100, 880, 960]); 
sg=sgtitle('Fault Structure Modeling');
set(sg,'fontsize',20,'fontweight','bold');
Rectangle_points=SUBfault_con(Optimal_C,f_sort_E2F,line_nums_2,hypotout,FaultpaRameters,Mc,minPoi,la0,lo0,colortemplate);

% --- Executes on button press in out_file.
function out_file_Callback(hObject, eventdata, handles)
global currentPath FE FaultpaRameters
fopen([currentPath,'/OutputFile/Fault_Segment_Clusterd.txt'],'w');
dlmwrite([currentPath,'/OutputFile/Fault_Segment_Clusterd.txt'], FE,'delimiter', ' ','precision', 8); % 将A写入指定文件，元素间默认以逗号分隔
disp(['Has Written in ',[currentPath,'/OutputFile/Fault_Segment_Clusterd.txt']])
fopen([currentPath,'/OutputFile/Fault_Segment_Modeling.txt'],'w');
Fault_Segment_Modeling=[FaultpaRameters(:,31),FaultpaRameters(:,1:3),FaultpaRameters(:,7:18),FaultpaRameters(:,25:27),FaultpaRameters(:,29:30)];
dlmwrite([currentPath,'/OutputFile/Fault_Segment_Modeling.txt'], Fault_Segment_Modeling,'delimiter', ' ','precision', 8); % 将A写入指定文件，元素间默认以逗号分隔
disp(['Has Written in ',[currentPath,'/OutputFile/Fault_Segment_Modeling.txt']])

% --- Executes on button press in OC_Segment_cluster.
function OC_Segment_cluster_Callback(hObject, eventdata, handles)
global f_sort_E2F line_nums_2 hypotout Mc minPoi la0 lo0 NEvents colortemplate Optimal_C
n=1;FE_km=[];
nw=1;RabY=[];c=Optimal_C;
for i=1:f_sort_E2F(end)
    im=1;
    if ismember(i,line_nums_2(:,3))
        k=find(line_nums_2(:,3)==i);
        ftest4=f_sort_E2F(line_nums_2(k,1):line_nums_2(k,2),:);
        while 1
            MaxMag=max(ftest4(:,4));
            if im==1
                [a,b,Y,CauMag]=Mmain(ftest4,Mc,hypotout);
                im=2;
            elseif im==2
                CauMag=(log10(size(ftest4,1))-a)/b+Mc;
            end
            Mag=max(MaxMag,CauMag);
            Ns=(log10(size(ftest4,1))-a)/b+Mc;
            radius=c*radius_model(Mag);
            [~,s4]=sort(ftest4(:,4),'descend');
            ftest4=ftest4(s4,:);            
            idx_5=dbscan(ftest4(:,1:3),radius,minPoi);
            nw=nw+1;
            if all(idx_5 == -1) || isempty(idx_5)
                RabY=[RabY;i,a,b,mean(Y(:,1)),size(ftest4,1),MaxMag,CauMag,Ns,Mag,0,];
                break
            end
            RabY=[RabY;i,a,b,mean(Y(:,1)),size(ftest4,1),MaxMag,CauMag,Ns,Mag,1];
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

FE=KMtoDEG(FE_km,la0,lo0);
figure
ax6=gca;
if ismember(FE(:,1)<0,1)
    scatter3(abs(FE(:, 1)), FE(:, 2), FE(:, 3), 8, FE(:,11:13),"filled");
    set(ax6,'Zdir','reverse','Xdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12,'color',[0.9 0.9 0.9],'box','on','linewidth',1.5,'XGrid','off','YGrid','off','ZGrid','off');
    ax6.XTickLabel=strcat(ax6.XTickLabel, '°W');
else
    scatter3(FE(:, 1), FE(:, 2), FE(:, 3),8, FE(:,11:13),"filled");
    set(ax6,'Zdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12,'color',[0.9 0.9 0.9],'box','on','linewidth',1.5,'XGrid','off','YGrid','off','ZGrid','off');
    ax6.XTickLabel=strcat(ax6.XTickLabel, '°E');
end
d1=(max(hypotout(:,9))-min(hypotout(:,9)))/50;
d2=(max(hypotout(:,8))-min(hypotout(:,8)))/50;
if min(FE(:,1))<0
    xlim(ax6,[abs(max(hypotout(:,9))+d1) abs(min(hypotout(:,9))-d1)]);
    ylim(ax6,[min(hypotout(:,8))-d2 max(hypotout(:,8))+d2]); 
else
    xlim(ax6,[min(hypotout(:,9))-d1 max(hypotout(:,9))+d1]);
    ylim(ax6,[min(hypotout(:,8))-d2 max(hypotout(:,8))+d2]);
end
view(ax6,0,90)
ax6.YTickLabel=strcat(ax6.YTickLabel, '°N');
pos6 = get(ax6, 'Position');
annotation('textbox', [pos6(1)+0.5 pos6(2), pos6(3), 0.2*pos6(4)], 'LineStyle','none',...
            'String', ['Fault lines=',num2str(length(unique(FE(:,14)))),newline,'Events=',num2str(size(FE,1)),newline,'Remaining=',num2str(NEvents-size(FE,1)) ], ...
            'FontAngle','normal','FontSize', 12, 'Color', 'black'); 
title(ax6,sprintf('Optimal C = %d',c))



% --- Executes during object creation, after setting all properties.
function main_fault_orientation_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to main_fault_orientation_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on key press with focus on PBAD_Multiple_edit and none of its controls.
function PBAD_Multiple_edit_KeyPressFcn(hObject, eventdata, handles)
% global na1
% inputStr = eventdata.Character; % 获取 edit text 的字符串内容
% na1 = str2double(inputStr)% 将字符串转换为数值


% --- Executes on key press with focus on MAD_m2_edit and none of its controls.
function MAD_m2_edit_KeyPressFcn(hObject, eventdata, handles)
% global na2
% inputStr = eventdata.Character; % 获取 edit text 的字符串内容
% na2 = str2double(inputStr)% 将字符串转换为数值


% --- Executes on key press with focus on main_fault_orientation_edit and none of its controls.
function main_fault_orientation_edit_KeyPressFcn(hObject, eventdata, handles)
% global main_f main_fnum
% inputStr = get(hObject, 'String'); % 获取 edit text 的字符串内容
% main_f = str2double(split(inputStr,',')) % 将字符串转换为数值
% main_fnum=length(main_f)


% --- Executes during object creation, after setting all properties.
function axes6_CreateFcn(hObject, eventdata, handles)
global ax6
ax6 = hObject; % 将 hObject 赋值给全局变量 ax1
set(ax6,'ZDir','reverse','FontSize',15,'Color',[0.9 0.9 0.9],'LineWidth',1.5);
box on;grid off
xticks(ax6,[]);yticks(ax6,[]);


% --- Executes during object creation, after setting all properties.
function STEP0_text_CreateFcn(hObject, eventdata, handles)
% parentBackgroundColor = get(hObject.Parent, 'BackgroundColor');
% set(hObject, 'BackgroundColor', 'none');


% --- Executes during object creation, after setting all properties.
function STEP1_text_CreateFcn(hObject, eventdata, handles)
% set(hObject, 'BackgroundColor', 'none');


% --- Executes during object creation, after setting all properties.
function STEP2_text_CreateFcn(hObject, eventdata, handles)
% set(hObject, 'BackgroundColor', 'none');


% --- Executes during object creation, after setting all properties.
function STEPEND_text_CreateFcn(hObject, eventdata, handles)
% set(hObject, 'BackgroundColor', 'none');



function angle_tolerance_edit_Callback(hObject, eventdata, handles)
global angle_tolerance
inputStr = get(hObject, 'String'); % 获取 edit text 的字符串内容
angle_tolerance = str2double(inputStr) 
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function angle_tolerance_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_tolerance_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Preferred_Orientation.
function Preferred_Orientation_Callback(hObject, eventdata, handles)
global Optimal_C FaultpaRameters
Azimuth=FaultpaRameters(FaultpaRameters(:,29)==Optimal_C,25);
Azimuth(find(Azimuth<0))=Azimuth(find(Azimuth<0))+180;
[counts, edges, number]=histcounts(Azimuth,'BinEdges',0:10:180)
binCenters = edges(1:end-1) + diff(edges)/2
abnum=FaultpaRameters(:,30);
Aznum=[];
for j=1:size(counts,2)
    Aznum=[Aznum,sum(abnum(find(number==j)))]
end
figure 
axPO=gca;
% title('Preferred Fault Orientation');
set(axPO,'color',[0.9 0.9 0.9],'fontsize',12,'LineWidth',1.5,'box','on')
for i=1:length(counts)
    hold on
    bar(binCenters(i), Aznum(i), 10, 'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5);
end
xticks(0:10:180)
xlabel('Azimuth')
ylabel('Event counts')
