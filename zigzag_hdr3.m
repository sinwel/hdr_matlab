function out = zigzag_hdr3(in, param, times)%, handles)
% debug
if nargin < 1
%      in = floor(rand(4,48)*1023);     
in = [650,	655,	875	,889,	632,	123,	296	,1017	,668	,479	,229	,854	,389	,57	  371,	9	 ,  517,	650,	401,	269,	97	,997,	902	,845	,745	,223	,315	,603	,232	,151	,754	,156	,508	53	,305	,389	,316	,924	,939	,520	,965	,895	,600	,100	,33	  ,317	,620	,779,  
718,	548,	48	,382,	771,	477,	851	,528	,178	,656	,260	,266	,160	,6	  244,	576,	964,	625,	777,	786,	970	,768,	373	,249	,703	,1000	,251	,688	,789	,931	,903	,441	,144	1003,	274	,766	,858	,731	,944	,731	,740	,44	  ,727	,936	,1018	,300	,96  	,271,  
784,	533,	419	,729,	494,	299,	351	,585	,700	,793	,136	,907	,530	,283	919,	322,	901,	136,	850,	947,	501	,644,	767	,490	,184	,578	,440	,278	,918	,185	,84	  ,675	,827	616	,846	,322	,644	,92	  ,666	,275	,919	,920	,499	,640	,429	,516	,26	  ,629,  
348,	972,	830	,682,	407,	923,	67	,678	,171	,293	,389	,238	,381	,171	227,	870,	587,	625,	143,	160,	796	,598,	1018,	39	,526	,335	,278	,994	,125	,177	,622	,474	,293	546	,938	,611	,167	,820	,619	,474	,208	,791	,42	  ,230	,379	,779	,599	,311,  
];
 
     param.zigzagpattern = [1 0 1 1;
                     0 0 1 0;
                     1 1 1 0;
                     1 0 0 0;];
    param.bits      = 10;
    param.noise     = 64;
    param.exptimes  = 8;
    times           = param.exptimes;
    param.blacklevel = 64;
    
    in = in + param.blacklevel;
end
if 0%get(handles.checkbox4, 'Value')
    in = defect_pixel_processhdr(in);
end
% -------------------------------------------------------------------------
% get fushion image by interpolation with R/B & G bayer.
in = in - param.blacklevel;
normlizeValue = 2^param.bits - param.blacklevel;
% in = in/(2^param.bits - param.blacklevel);



F = repmat(param.zigzagpattern, [size(in,1)/4, size(in,2)/4]);


l_image = in.*F;
s_image = in.*(1-F);

image = l_image+s_image*times;

tl1 = [1 0 0; 0 0 0; 0 0 -1];
% tl2 = [1 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 0];
% tl3 = [0 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 1];

tr1 = [0 0 1; 0 0 0; -1 0 0];
% tr2 = [0 0 0 0 1; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 0];
% tr3 = [0 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 1 0 0 0 0];

% thv1 = [1 0 -1];
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

l_ref_image = l_image + f_image.*(1-F);
s_ref_image = s_image*times + f_image.*(F);

% % debug %%
if 1
    image_debug = in;
    bcm_r = repmat([0 1; 0 0], [size(image_debug,1)/2, size(image_debug,2)/2]);
    bcm_b = repmat([0 0; 1 0], [size(image_debug,1)/2, size(image_debug,2)/2]);
    bcm_g = repmat([1 0; 0 1], [size(image_debug,1)/2, size(image_debug,2)/2]);
    % get r/b
    I_r = image_debug.*bcm_r;
    I_r_  = I_r';
%     rIdx = find(I_r_~=0);
%     I_r_packed = reshape(I_r_(rIdx), size(image_debug,2)/2, size(image_debug,1)/2)';
    I_b = image_debug.*bcm_b;
    I_b_  = I_b';
%     bIdx = find(I_b_~=0);
%     I_b_packed = reshape(I_b_(bIdx), size(image_debug,2)/2, size(image_debug,1)/2)';
    % get g
    I_g  = image_debug.*(bcm_g);
    I_g_  = I_g';
%     gIdx = find(I_g_~=0);
%     I_g_packed = reshape(I_g_(gIdx), size(image,2)/2, size(image,1))';
    % swap 1,3 line, long time togethter, short time together.
%     tmp = I_g_packed(2,:);
%     I_g_packed(2,:) = I_g_packed(3,:);
%     I_g_packed(3,:) = tmp;
    
   
    thv2_ = [1 0 0 0 -1];
    hh_ = [1 0 0  0 1];
    
    tl1_ = [1 0  0; 0 0 0; 0  0  -1 ];    
    tr1_ = [0 0  1 ; 0 0 0; -1 0  0]; 
    
    hl_ = [1 0 0  ; 0 0 0;  0 0 1  ];
    hr_ = [ 0 0 1 ; 0 0 0;  1 0 0];
    % red
    r_mode_h = abs(conv2(I_r, thv2_,  'same')) ;%+ abs(conv2(image, thv1,  'same'));
    r_mode_v = abs(conv2(I_r, thv2_', 'same')) ;%+ abs(conv2(image, thv1', 'same'));

    
    r_h = conv2(I_r, hh_, 'same');
    r_v = conv2(I_r, hh_', 'same');
	r_choose_v = r_mode_h>r_mode_v;
    fushion_r= r_choose_v.*r_v+(1-r_choose_v).*r_h;
    
    % blue
    b_mode_h = abs(conv2(I_b, thv2_,  'same')) ;%+ abs(conv2(image, thv1,  'same'));
    b_mode_v = abs(conv2(I_b, thv2_', 'same')) ;%+ abs(conv2(image, thv1', 'same'));

    
    b_h = conv2(I_b, hh_, 'same');
    b_v = conv2(I_b, hh_', 'same');
	b_choose_v = b_mode_h>b_mode_v;
    fushion_b= b_choose_v.*b_v+(1-b_choose_v).*b_h;
    
    % green
    g_mode_l = abs(conv2(I_g, tl1_, 'same')) ;%+ abs(conv2(image, thv1,  'same'));
    g_mode_r = abs(conv2(I_g, tr1_, 'same')) ;%+ abs(conv2(image, thv1', 'same'));

    
    g_l = conv2(I_g, hl_, 'same');
    g_r = conv2(I_g, hr_, 'same');
	g_choose_r = g_mode_l > g_mode_r;
    fushion_g= g_choose_r.*g_r + (1-g_choose_r).*g_l;
    
    
    fushion = fushion_r + fushion_b + fushion_g;
    PP = 8*(F) + (1-F); % 长曝光插值是用short得到的需要倍乘8，短曝光是long插值得到反而不需要。
    find((fushion.*PP)~=(f_image*2))


    ll_image = l_image + (fushion/2).*(1-F); % orginal long + orignal(short)but interpolation long.
    ss_image = (s_image + (fushion/2).*(F))*times;

    find(ll_image~=l_ref_image)
    find(ss_image~=s_ref_image)

end

% -------------------------------------------------------------------------

if 0%get(handles.checkbox2, 'Value')
    imwrite(uint16(l_image*1023), 'hdr_long_exp.pgm');
    imwrite(uint16(s_image*1023/times), 'hdr_short_exp.pgm');
end

if 0%get(handles.checkbox1, 'Value')
    figure;imshow((double(demosaic(uint8(l_ref_image), 'grbg'))/255).^0.45);
    figure;imshow((double(demosaic(uint8(s_ref_image/times),'grbg'))/255).^0.45);
end

% t = [1 2 1; 2 4 2; 1 2 1]/16;
% t = [1];
% fs_image = conv2(s_image, t, 'same');
% fl_image = conv2(l_image, t, 'same');

% -------------------------------------------------------------------------
% use difference for tonning the image
% 
% 
D = (s_ref_image-l_ref_image);
D = abs(D);
D = D*(2^param.bits/(param.noise*param.exptimes));
D = min(D,normlizeValue);

%% fix point
% debug 
if 1
    tone_map = importdata('tone_mapping_961.dat');    
    tone_map = reshape(tone_map',1,numel(tone_map));
    table_map = tone_map(1:961);
    D_scale = table_map(D+1);

end
%% Float data
d = 1;
b = 3+d;
c = -2-2*d;

% x = 0:0.01:1;
% figure;plot(b*x.^2+c*x.^3+d*x.^4);
% D = b*D.^2+c*D.^3+d*D.^4;
D_float = D/normlizeValue;
D_scale_ = (b*D_float.^2+c*D_float.^3+d*D_float.^4)*normlizeValue;
D_scale_ = min(D_scale_,normlizeValue);
D_scale_ = max(D_scale_,0);





D_scale = max(D_scale, double(s_ref_image>0.9*normlizeValue)*normlizeValue );
D_scale = min(D_scale, (1 - double(s_ref_image < (0.85*normlizeValue)/times))*normlizeValue);

% t=fspecial('gaussian',[9 9], 2);
% t = ones(3,3)/9;
% t = [1 0 1 0 1; 0 0 0 0 0; 1 0 1 0 1; 0 0 0 0 0; 1 0 1 0 1]/9;
% t = [1];
% D = conv2(D, t, 'same');
D_adj = ordfilt2(D_scale, 9, ones(3,3));
% figure;imshow(D);
% out = l_ref_image.*(1-(D_adj/normlizeValue))+s_ref_image.*(D_adj/normlizeValue);
out = (l_ref_image.*(normlizeValue-D_adj)+s_ref_image.*D_adj)/normlizeValue;
out_demosic = demosaic(uint16(out), 'grbg');
imwrite(uint16(out_demosic*64), 'out16bit.png');
figure;imshow(double(demosaic(uint8(out*255/1023), 'grbg'))/255);
if 0%get(handles.checkbox5, 'Value')
    L = ordfilt2(out, 9, ones(3,3));
    L = L/(param.exptimes);
    L = normlizeValue*(1-(1 - L/normlizeValue).^4);
    L = min(L,normlizeValue);
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

    out = out.*S;
end

