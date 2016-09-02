
function out = zigzag_hdr2(in, param, times, handles)

if get(handles.checkbox4, 'Value')
    in = defect_pixel_processhdr(in);
end



% 插值
[l_image, s_image] = conv2_test_ceva_xm4(in,param);



if get(handles.checkbox2, 'Value')
    imwrite(uint16(l_image*1023), 'hdr_long_exp.ppm');
    imwrite(uint16(s_image*1023/times), 'hdr_short_exp.ppm');
end

if get(handles.checkbox1, 'Value')
    figure;imshow((double(demosaic(uint8(l_image*255), 'grbg'))/255).^0.45);
    figure;imshow((double(demosaic(uint8(255*s_image/times),'grbg'))/255).^0.45);
end

if 0;%get(handles.checkbox1, 'Value')
    figure;imshow((double(demosaic(uint8(l_image/4), 'grbg'))/255).^0.45);
    figure;imshow((double(demosaic(uint8(s_image/(4*times)),'grbg'))/255).^0.45);
end

% t = [1 2 1; 2 4 2; 1 2 1]/16;
% t = [1 2 1; 2 4 2; 1 2 1];
% fs_image = conv2(double(s_image), double(t), 'same');
% fl_image = conv2(double(l_image), double(t), 'same');
%% 查表 调整
D = (s_image-l_image);
gain = 2^param.bits/(param.noise*param.exptimes);
D = D*gain;

out = tonning_map_ceva_xm4(D,l_image,s_image,param,Max_datain,0.9, 0.85,times);
figure;imshow(double(out)/Max_datain);
%% WDR modual.
if get(handles.checkbox5, 'Value')
    L = ordfilt2(out, 9, ones(3,3));
    L = L/(param.exptimes);
    L = 1-(1-L).^4;
    L = min(L,1);
    t = fspecial('gaussian',[5 5], 0.9);
%     L = conv2(L,t,'same');

    fL = bilateralFilter(L,L,0,1,8,0.125);
%     t = fspecial('gaussian',[1 300], 90);
%     fL = conv2(L, t, 'same');
%     fL = conv2(fL, t', 'same');
%     fL = conv2(fL,t,'same');
%     fL = L;
    X = min(fL, L+0.08);
    X = max(X, L-0.08);
%     figure;imshow(abs(X));
%     figure;imshow(abs(L-fL));

    
    bx = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
    mx = [0 16 15 14 13 12 11.4 10.5 9.7 9.2 8.7 8.2 7.7 7.2 6.78 6.44 6.18 5.93 5.72 5.55 5.4 5.25 5.1 4.95];
    e = [0 -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18 -19 -20 -21 -22 -23];

    x = 0:1/(1024-64):1;
    d = bx-1-2*e;
    c = 2-2*bx+e;    
    i = fix(param.exptimes);

    if param.wdrgain
        i = min(i,24);
        i = min(i,param.wdrgain);
    else
        i = fix(param.exptimes/2);
    end

    b = bx(i);
    p = mx(i);
    S = b+c(i)*X.^(p/16)+d(i)*X.^(2*p/16)+e(i)*X.^(3*p/16);    

% x = 0:0.01:1;
% figure;plot(b*x+c(i)*x.^(1+p/16)+d(i)*x.^(1+2*p/16)+e(i)*x.^(1+3*p/16));
    
%     figure;imshow(X);
%     S = X.^0.45./X;
%     S = (1-(1-X).^param.exptimes)./X;
%     S(isnan(S)) = 1;

    out = out.*S;
end

