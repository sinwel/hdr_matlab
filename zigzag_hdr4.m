function out = zigzag_hdr4(in, param, times, handles)

% if get(handles.checkbox4, 'Value')
%     in = defect_pixel_processhdr(in);
% end

in = in - param.blacklevel;
in = in/(2^param.bits - param.blacklevel);

F = repmat(param.zigzagpattern, [size(in,1)/4, size(in,2)/4]);


l_image = in.*F;
s_image = in.*(1-F);

image = l_image+s_image*times;

tl1 = [1 0 0; 0 0 0; 0 0 -1];
tl2 = [1 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 0];
tl3 = [0 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 1];

tr1 = [0 0 1; 0 0 0; -1 0 0];
tr2 = [0 0 0 0 1; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 0];
tr3 = [0 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 1 0 0 0 0];

thv1 = [1 0 -1];
thv2 = [1 0 0 0 -1];

pl = abs(conv2(image, tl1, 'same')) ;%+ abs(conv2(image, tl2, 'same')) + abs(conv2(image, tl3, 'same'));
pr = abs(conv2(image, tr1, 'same')) ;% + abs(conv2(image, tr2, 'same')) + abs(conv2(image, tr3, 'same'));
ph = abs(conv2(image, thv2,  'same'));% + abs(conv2(image, thv1,  'same'));
pv = abs(conv2(image, thv2', 'same'));% + abs(conv2(image, thv1', 'same'));

hl = [1 0 0; 0 0 0; 0 0 1]/2;
hr = [0 0 1; 0 0 0; 1 0 0]/2;
hh = [1 0 0 0 1]/2;

fl = conv2(image, hl, 'same');
fr = conv2(image, hr, 'same');
fh = conv2(image, hh, 'same');
fv = conv2(image, hh', 'same');

FF = pl>pr;
fg = FF.*fr+(1-FF).*fl;
FF = ph>pv;
frb= FF.*fv+(1-FF).*fh;

TF = repmat([1 0; 0 1], [size(in,1)/2, size(in,2)/2]);
f_image = fg.*TF + frb.*(1-TF);

l_image = l_image + f_image.*(1-F);
s_image = s_image*times + f_image.*(F);

if get(handles.checkbox2, 'Value')
    imwrite(uint16(l_image*1023), 'hdr_long_exp.pgm');
    imwrite(uint16(s_image*1023/times), 'hdr_short_exp.pgm');
end

if get(handles.checkbox1, 'Value')
    figure;imshow((double(demosaic(uint8(l_image*255), 'grbg'))/255).^0.45);
    figure;imshow((double(demosaic(uint8(255*s_image/times),'grbg'))/255).^0.45);
end

t = [1 2 1; 2 4 2; 1 2 1]/16;
t = [1];
fs_image = conv2(s_image, t, 'same');
fl_image = conv2(l_image, t, 'same');


D = (s_image-l_image);
D = abs(D);
D = D*(2^param.bits/(param.noise*param.exptimes));
D = min(D,1);

d = 1;
b = 3+d;
c = -2-2*d;

% x = 0:0.01:1;
% figure;plot(b*x.^2+c*x.^3+d*x.^4);
D = b*D.^2+c*D.^3+d*D.^4;

D = min(D,1);
D = max(D,0);

D = max(D, double(s_image>0.9));
D = min(D, 1-double(s_image<0.85/times));

t=fspecial('gaussian',[9 9], 2);
t = ones(3,3)/9;
t = [1 0 1 0 1; 0 0 0 0 0; 1 0 1 0 1; 0 0 0 0 0; 1 0 1 0 1]/9;
t = [1];
D = conv2(D, t, 'same');
D = ordfilt2(D, 9, ones(3,3));
% figure;imshow(D);

out = l_image.*(1-D)+s_image.*D;

D = double(f_image.*(F) > 1);

out = image.*(1-D) + f_image.*D;

if get(handles.checkbox5, 'Value')
    L = ordfilt2(out, 9, ones(3,3));
    L = L/(param.exptimes);
    L = 1-(1-L).^4;
    L = min(L,1);
    t = fspecial('gaussian',[5 5], 0.9); 
%     L = conv2(L,t,'same');

%     fL = bilateralFilter(L,L,0,1,8,0.125);
    t = fspecial('gaussian',[1 256], 96);
    fL = conv2(L, t, 'same');
    fL = conv2(fL, t', 'same');
%     fL = conv2(fL,t,'same');
%     fL = L;
    X = min(fL, L+0.125);
    X = max(X, L-0.125);
%     figure;imshow(abs(X));
%     figure;imshow(abs(L-fL));

    
    bx = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
    mx = [16 16 15 14 13 12 11.4 10.5 9.7 9.2 8.7 8.2 7.7 7.2 6.78 6.44 6.18 5.93 5.72 5.55 5.4 5.25 5.1 4.95];
    e = [0 -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -16 -17 -18 -19 -20 -21 -22 -23];

    x = 0:1/(1024-64):1;
    d = bx-1-2*e+(1-bx)./(2*bx.*mx/16);
    c = 2-2*bx+e-(1-bx)./(2*bx.*mx/16);    
    
%     d = bx-1-2*e;
%     c = 2-2*bx+e;  
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
