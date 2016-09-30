function out = zigzag_hdr3_float(in, param, times, handles)

if 0% get(handles.checkbox4, 'Value')
    in = defect_pixel_processhdr(in);
end

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
pr = abs(conv2(image, tr1, 'same'))  ;%+ abs(conv2(image, tr2, 'same')) + abs(conv2(image, tr3, 'same'));
ph = abs(conv2(image, thv2,  'same')) ;%+ abs(conv2(image, thv1,  'same'));
pv = abs(conv2(image, thv2', 'same')) ;%+ abs(conv2(image, thv1', 'same'));

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

if 0% get(handles.checkbox2, 'Value')
    imwrite(uint16(l_image*1023), 'hdr_long_exp.pgm');
    imwrite(uint16(s_image*1023/times), 'hdr_short_exp.pgm');
end

if 0% get(handles.checkbox1, 'Value')
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
avg_h = [1 2 1; 2 4 2; 1 2 1]/16;
D_filtered = conv2(D, avg_h, 'same');

D = conv2(D, t, 'same');
D = ordfilt2(D, 9, ones(3,3));
% figure;imshow(D);

out = l_image.*(1-D)+s_image.*D;

out_thumb = l_image.*(1-D_filtered)+s_image.*D_filtered;

figure;imshow(double(demosaic(uint8(out*255), 'grbg'))/255);
out_demosic = demosaic(uint16(out*1023), 'grbg');
imwrite(uint16(out_demosic*64), 'out16bit_float.png');

figure;imshow(double(demosaic(uint8(out_thumb*255), 'grbg'))/255);
if 1% get(handles.checkbox5, 'Value')
    L = ordfilt2(out, 9, ones(3,3));
    L = L/(param.exptimes);
    L = 1-(1-L).^4;
    L = min(L,1);
    %%
%     t = fspecial('gaussian',[5 5], 0.9);
%     L = conv2(L,t,'same');

%     fL = bilateralFilter(L,L,0,1,8,0.125);
    t = fspecial('gaussian',[1 32], 96);
    fL = conv2(L, t, 'same');
    fL = conv2(fL, t', 'same');
    fL = imresize(fL, 0.0625);
    fL = imresize(fL, size(L));
%     fL = conv2(fL,t,'same');
%     fL = L;
    X = min(fL, L+0.125);
    X = max(X, L-0.125);
%     figure;imshow(abs(X));
%     figure;imshow(abs(L-fL));
figure;imshow(double(demosaic(uint8(X*255), 'grbg'))/255);

   
    
    k = fix(param.exptimes);

    if param.wdrgain
        k = min(k,24);
        k = min(k,param.wdrgain);
    else
        k = fix(param.exptimes/2);
    end

    n= 0.05;

    x0 = 1/k/32;
    d = (k-1) / ((1+n)*x0^n - 2*(1+2*n)*x0^(2*n) + (1+3*n)*x0^(3*n));
    S = 1+d*X.^(n)-2*d*X.^(2*n)+d.*X.^(3*n);

    out1 = out.*double(S);
    figure;imshow(double(demosaic(uint8(out1*255), 'grbg'))/255);
    out_demosic1 = demosaic(uint16(out1*1023), 'grbg');
    imwrite(uint16(out_demosic1*64), 'outWDR_16bit_float.png');
    %%
    % 16x16 block downscale
    scale_gap = 16;
    for m=1:scale_gap:1504
        for n=1:scale_gap:2016
            out_blk = out_thumb(m:m+scale_gap-1,n:n+scale_gap-1);
            thumb(floor(m/scale_gap)+1,floor(n/scale_gap+1)) = mean(out_blk(:));
        end
    end
    thumb_S = thumb/(param.exptimes);
    thumb_S = 1-(1-thumb_S).^4;
    thumb_S = min(thumb_S,1);
    t = fspecial('gaussian',[1 3], 96);

    thumb_S = conv2(thumb_S, t, 'same');
    thumb_S = conv2(thumb_S, t', 'same');

    thumb_S = imresize(thumb_S, size(out));
    thumb_S  = min(thumb_S, thumb_S+0.125);
    thumb_S1 = max(thumb_S, thumb_S-0.125);
    
figure;imshow(double(demosaic(uint8(thumb_S1*255), 'grbg'))/255);
    k = fix(param.exptimes);

    if param.wdrgain
        k = min(k,24);
        k = min(k,param.wdrgain);
    else
        k = fix(param.exptimes/2);
    end

    n= 0.05;

    x0 = 1/k/32;
    d = (k-1) / ((1+n)*x0^n - 2*(1+2*n)*x0^(2*n) + (1+3*n)*x0^(3*n));
    S = 1+d*thumb_S1.^(n)-2*d*thumb_S1.^(2*n)+d.*thumb_S1.^(3*n);
    out2 = out.*double(S);
    figure;imshow(double(demosaic(uint8(out2*255), 'grbg'))/255);
end