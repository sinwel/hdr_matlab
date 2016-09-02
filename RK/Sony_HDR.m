function varargout = Sony_HDR(varargin)
% SONY_HDR MATLAB code for Sony_HDR.fig
%      SONY_HDR, by itself, creates a new SONY_HDR or raises the existing
%      singleton*.
%
%      H = SONY_HDR returns the handle to a new SONY_HDR or the handle to
%      the existing singleton*.
%
%      SONY_HDR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SONY_HDR.M with the given input arguments.
%
%      SONY_HDR('Property','Value',...) creates a new SONY_HDR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Sony_HDR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Sony_HDR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Sony_HDR

% Last Modified by GUIDE v2.5 11-Jul-2016 14:59:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Sony_HDR_OpeningFcn, ...
                   'gui_OutputFcn',  @Sony_HDR_OutputFcn, ...
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


% --- Executes just before Sony_HDR is made visible.
function Sony_HDR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Sony_HDR (see VARARGIN)

% Choose default command line output for Sony_HDR
handles.output = hObject;


% handles.zigzagpattern = ['L', 'S', 'L', 'L'; 'S', 'S', 'L' 'S'];
% 
% set(handles.uitable1, 'Data', handles.zigzagpattern);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Sony_HDR wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Sony_HDR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.raw;*.pgm', 'All Image Files';...
          '*.*','All Files' });
set(handles.text7, 'String', [pathname filename]);
% Update handles structure
guidata(hObject, handles);      

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
param = update_parameters(hObject, eventdata, handles);
fid = fopen(get(handles.text7, 'String'));
image = fread(fid, [param.width param.height], 'uint16')';
fclose(fid);

%
% image_cut = image(param.height/2 - 64:param.height/2+63,param.width/2-128:param.width/2+127);
% fidw = fopen('raw_256x128.raw','wb');
% fwrite(fidw, image_cut, 'uint16')';
% fclose(fidw);

image = image(param.height/2 - 64:param.height/2+63,param.width/2-224:param.width/2+31);
param.height = 128;
param.width  = 256;
image = int16(image);
Ximage = zigzag_hdr2(image, param, int16(param.exptimes), handles);
% Ximage = zigzag_hdr2_float(image, param, param.exptimes, handles);


Y = demosaic(uint16(Ximage*65535/param.exptimes), get(get(handles.uipanel1, 'SelectedObject'), 'String'));
Y = double(Y)/65535;
figure;imshow(Y.^0.45);
if get(handles.checkbox4, 'Value')
%     b=zeros(48,545);fig=figure('visible','off');imshow(b);text(0,24,'RockChip Demo for Asus', 'color','w','fontsize',48);s=getframe;mask=frame2im(s);close(fig);
%     mim = mask(:,:,1)*param.exptimes;
%     Ximage(fix(size(Ximage,1)*7/8):fix(size(Ximage,1)*7/8)+size(mim,1)-1, fix(size(Ximage,2)*5/8):fix(size(Ximage,2)*5/8)+size(mim,2)-1) = Ximage(fix(size(Ximage,1)*7/8):fix(size(Ximage,1)*7/8)+size(mim,1)-1, fix(size(Ximage,2)*5/8):fix(size(Ximage,2)*5/8)+size(mim,2)-1) + mim;
%     Ximage = min(Ximage, param.exptimes);
    imwrite(uint16(Ximage*1023/param.exptimes), 'hdr_out.pgm');
end




function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cfg = update_parameters(hObject, eventdata, handles)
cfg.width       = str2num(get(handles.edit3, 'String'));
cfg.height      = str2num(get(handles.edit4, 'String'));
cfg.blacklevel  = str2num(get(handles.edit1, 'String'));
cfg.exptimes    = str2num(get(handles.edit2, 'String'));
cfg.bits        = str2num(get(handles.edit5, 'String'));
cfg.noise       = str2num(get(handles.edit10, 'String'));
cfg.wdrgain     = str2num(get(handles.edit11, 'String'));

cfg.zigzagpattern(1,1) = strcmp(get(handles.text9, 'String'), 'L');
cfg.zigzagpattern(1,2) = strcmp(get(handles.text10, 'String'), 'L');
cfg.zigzagpattern(1,3) = strcmp(get(handles.text11, 'String'), 'L');
cfg.zigzagpattern(1,4) = strcmp(get(handles.text12, 'String'), 'L');
cfg.zigzagpattern(2,1) = strcmp(get(handles.text13, 'String'), 'L');
cfg.zigzagpattern(2,2) = strcmp(get(handles.text14, 'String'), 'L');
cfg.zigzagpattern(2,3) = strcmp(get(handles.text15, 'String'), 'L');
cfg.zigzagpattern(2,4) = strcmp(get(handles.text16, 'String'), 'L');
cfg.zigzagpattern(3,1) = strcmp(get(handles.text17, 'String'), 'L');
cfg.zigzagpattern(3,2) = strcmp(get(handles.text18, 'String'), 'L');
cfg.zigzagpattern(3,3) = strcmp(get(handles.text19, 'String'), 'L');
cfg.zigzagpattern(3,4) = strcmp(get(handles.text20, 'String'), 'L');
cfg.zigzagpattern(4,1) = strcmp(get(handles.text21, 'String'), 'L');
cfg.zigzagpattern(4,2) = strcmp(get(handles.text22, 'String'), 'L');
cfg.zigzagpattern(4,3) = strcmp(get(handles.text23, 'String'), 'L');
cfg.zigzagpattern(4,4) = strcmp(get(handles.text24, 'String'), 'L');


% 
% function out = zigzag_hdr2(in, param, times, handles)
% 
% if get(handles.checkbox4, 'Value')
%     in = defect_pixel_processhdr(in);
% end
% 
% in = in - param.blacklevel;
% in = in/(2^param.bits - param.blacklevel);
% 
% F = repmat(param.zigzagpattern, [size(in,1)/4, size(in,2)/4]);
% 
% 
% l_image = in.*F;
% s_image = in.*(1-F);
% 
% image = l_image+s_image*times;
% 
% tl1 = [1 0 0; 0 0 0; 0 0 -1];
% tl2 = [1 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 0];
% tl3 = [0 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 1];
% 
% tr1 = [0 0 1; 0 0 0; -1 0 0];
% tr2 = [0 0 0 0 1; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 0];
% tr3 = [0 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 1 0 0 0 0];
% 
% thv1 = [1 0 -1];
% thv2 = [1 0 0 0 -1];
% 
% pl = abs(conv2(image, tl1, 'same')) + abs(conv2(image, tl2, 'same')) + abs(conv2(image, tl3, 'same'));
% pr = abs(conv2(image, tr1, 'same')) + abs(conv2(image, tr2, 'same')) + abs(conv2(image, tr3, 'same'));
% ph = abs(conv2(image, thv1,  'same')) + abs(conv2(image, thv2,  'same'));
% pv = abs(conv2(image, thv1', 'same')) + abs(conv2(image, thv2', 'same'));
% 
% hl = [1 0 0; 0 0 0; 0 0 1]/2;
% hr = [0 0 1; 0 0 0; 1 0 0]/2;
% hh = [1 0 0 0 1]/2;
% 
% fl = conv2(image, hl, 'same');
% fr = conv2(image, hr, 'same');
% fh = conv2(image, hh, 'same');
% fv = conv2(image, hh', 'same');
% 
% FF = pl>pr;
% fg = FF.*fr+(1-FF).*fl;
% FF = ph>pv;
% frb= FF.*fv+(1-FF).*fh;
% 
% TF = repmat([1 0; 0 1], [size(in,1)/2, size(in,2)/2]);
% f_image = fg.*TF + frb.*(1-TF);
% 
% l_image = l_image + f_image.*(1-F);
% s_image = s_image*times + f_image.*(F);
% 
% if get(handles.checkbox2, 'Value')
%     imwrite(uint16(l_image*1023), 'hdr_long_exp.pgm');
%     imwrite(uint16(s_image*1023/times), 'hdr_short_exp.pgm');
% end
% 
% if get(handles.checkbox1, 'Value')
%     figure;imshow((double(demosaic(uint8(l_image*255), 'grbg'))/255).^0.45);
%     figure;imshow((double(demosaic(uint8(255*s_image/times),'grbg'))/255).^0.45);
% end
% 
% t = [1 2 1; 2 4 2; 1 2 1]/16;
% fs_image = conv2(s_image, t, 'same');
% fl_image = conv2(l_image, t, 'same');
% 
% 
% D = (s_image-l_image);
% D = abs(D);
% D = D*(2^param.bits/param.noise/param.exptimes);
% D = min(D,1);
% 
% d = 1;
% b = 3+d;
% c = -2-2*d;
% 
% % x = 0:0.01:1;
% % figure;plot(b*x.^2+c*x.^3+d*x.^4);
% D = b*D.^2+c*D.^3+d*D.^4;
% 
% D = min(D,1);
% D = max(D,0);
% 
% D = max(D, double(s_image>0.9));
% D = min(D, 1-double(s_image<0.85/times));
% 
% t=fspecial('gaussian',[9 9], 2);
% t = ones(3,3)/9;
% t = [1 0 1 0 1; 0 0 0 0 0; 1 0 1 0 1; 0 0 0 0 0; 1 0 1 0 1]/9;
% D = conv2(D, t, 'same');
% 
% % figure;imshow(D);
% 
% out = l_image.*(1-D)+s_image.*D;
% 
% if get(handles.checkbox5, 'Value')
%     L = ordfilt2(out, 9, ones(3,3));
%     L = L/(param.exptimes);
%     L = 1-(1-L).^4;
%     L = min(L,1);
%     t = fspecial('gaussian',[5 5], 0.9);
% %     L = conv2(L,t,'same');
% 
%     fL = bilateralFilter(L,L,0,1,8,0.125);
% %     t = fspecial('gaussian',[1 300], 90);
% %     fL = conv2(L, t, 'same');
% %     fL = conv2(fL, t', 'same');
% %     fL = conv2(fL,t,'same');
% %     fL = L;
%     X = min(fL, L+0.08);
%     X = max(X, L-0.08);
% %     figure;imshow(abs(X));
% %     figure;imshow(abs(L-fL));
% 
%     
%     bx = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
%     mx = [0 16 15 14 13 12 11.4 10.5 9.7 9.2 8.7 8.2 7.7 7.2 6.78 6.44 6.18 5.93 5.72 5.55 5.4 5.25 5.1 4.95];
%     e = [0 -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18 -19 -20 -21 -22 -23];
% 
%     x = 0:1/(1024-64):1;
%     d = bx-1-2*e;
%     c = 2-2*bx+e;    
%     i = fix(param.exptimes);
% 
%     if param.wdrgain
%         i = min(i,24);
%         i = min(i,param.wdrgain);
%     else
%         i = fix(param.exptimes/2);
%     end
% 
%     b = bx(i);
%     p = mx(i);
%     S = b+c(i)*X.^(p/16)+d(i)*X.^(2*p/16)+e(i)*X.^(3*p/16);    
% 
% % x = 0:0.01:1;
% % figure;plot(b*x+c(i)*x.^(1+p/16)+d(i)*x.^(1+2*p/16)+e(i)*x.^(1+3*p/16));
%     
% %     figure;imshow(X);
% %     S = X.^0.45./X;
% %     S = (1-(1-X).^param.exptimes)./X;
% %     S(isnan(S)) = 1;
% 
%     out = out.*S;
% end
% 


% 
% function out = defect_pixel_processhdr(I)
% X = I;
% %P = [0 0 1 0 0; 0 1 0 1 0; 1 0 1 0 1; 0 1 0 1 0; 0 0 1 0 0];
% P = [
%     0 0 0 0 1 0 0 0 0;
%     0 0 0 0 0 0 0 0 0;
%     0 0 1 0 0 0 1 0 0;
%     0 0 0 0 0 0 0 0 0;
%     1 0 0 0 1 0 0 0 1;
%     0 0 0 0 0 0 0 0 0;
%     0 0 1 0 0 0 1 0 0;
%     0 0 0 0 0 0 0 0 0;
%     0 0 0 0 1 0 0 0 0;
%     ];
% c1 = ordfilt2(X, 1, P);
% c2 = ordfilt2(X, 2, P);
% c8 = ordfilt2(X, 6, P);
% c9 = ordfilt2(X, 9, P);
% dis = c8-c2;
% t1 = (c8-c1) > 1.2*dis;
% t2 = (c9-c2) > 1.2*dis;
% t1 = t1.*(X == c1);
% t2 = t2.*(X == c9);
% t = double(t1|t2);
% X = X.*(1-t)+t1.*c2+t2.*c8;
% Y = X;
% 
% X = I;
% P = [
%     1 0 1 0 1;
%     0 0 0 0 0;
%     1 0 1 0 1;
%     0 0 0 0 0;
%     1 0 1 0 1;
%     ];
% c1 = ordfilt2(X, 1, P);
% c2 = ordfilt2(X, 2, P);
% c8 = ordfilt2(X, 6, P);
% c9 = ordfilt2(X, 9, P);
% dis = c8-c2;
% t1 = (c8-c1) > 1.2*dis;
% t2 = (c9-c2) > 1.2*dis;
% t1 = t1.*(X == c1);
% t2 = t2.*(X == c9);
% t = double(t1|t2);
% X = X.*(1-t)+t1.*c2+t2.*c8;
% 
% T = repmat([1 0; 0 1], size(I,1)/2, size(I,2)/2);
% out = X.*T+Y.*(1-T);
% 
% 

function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
