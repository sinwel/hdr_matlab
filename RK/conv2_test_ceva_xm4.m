function [fg_fushion, frb_fushion] = conv2_test_ceva_xm4(in,param);
% clear;
% clc;

% debug
in = int16([meshgrid(1:16,1:16) ,meshgrid(1:16,1:16) ] + 64);
param.width = 32;
param.height = 16;
padding = 2;
param.blacklevel = 64;
param.bits = 10;
param.exptimes = 8;
param.zigzagpattern = [ 1 0 1 1;
                        0 0 1 0
                        1 1 1 0
                        1 0 0 0];
in = in - param.blacklevel;
Max_datain = 2^param.bits - param.blacklevel;
% in = in/(2^param.bits - param.blacklevel); % 归一化 0 -1
max(in(:))

F = repmat(param.zigzagpattern, [size(in,1)/4, size(in,2)/4]);
width = size(in,2);
height= size(in,1);

l_image = in.*int16(F);% long exposue pic
s_image = in.*int16(1-F);% short exposure pic

image = l_image+s_image*param.exptimes;
max(image(:))


image_padding = zeros(param.height+2*padding, param.width+2*padding);

tl1 = [1 0 0; 0 0 0; 0 0 -1];
tl2 = [1 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 0];
tl3 = [0 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 1];

tr1 = [0 0 1; 0 0 0; -1 0 0];
tr2 = [0 0 0 0 1; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 0];
tr3 = [0 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 1 0 0 0 0];

thv1 = [1 0 -1];
thv2 = [1 0 0 0 -1];


% padding = 2;
% image = floor(rand(8,8)*255);
% image =  meshgrid(1:width,1:height);

pl = int16(abs(conv2(double(image), tl1, 'same')) + abs(conv2(double(image), tl2, 'same')) + abs(conv2(double(image), tl3, 'same')));
pr = int16(abs(conv2(double(image), tr1, 'same')) + abs(conv2(double(image), tr2, 'same')) + abs(conv2(double(image), tr3, 'same')));
G_choose_r  =  pl > pr;
ph = int16(abs(conv2(double(image), thv1,  'same')) + abs(conv2(double(image), thv2,  'same')));
pv = int16(abs(conv2(double(image), thv1', 'same')) + abs(conv2(double(image), thv2', 'same')));
RB_choose_v =  ph > pv;

hl = [1 0 0; 0 0 0; 0 0 1]/2;
hr = [0 0 1; 0 0 0; 1 0 0]/2;
hh = [1 0 0 0 1]/2;

fl = uint16(conv2(double(image), hl, 'same'));
fr = uint16(conv2(double(image), hr, 'same'));
fh = uint16(conv2(double(image), hh, 'same'));
fv = uint16(conv2(double(image), hh', 'same'));

%% 卷积有旋转的操作，导致left对应了right

fg = int16(G_choose_r).*int16(fr)+int16(1-G_choose_r).*int16(fl);
 
frb= int16(RB_choose_v).*int16(fv)+int16(1-RB_choose_v).*int16(fh);


if 1
% read the data, and get index address to calculate.

% rotate the kernel.
tl1_mpy = [-1 0 0; 0 0 0; 0 0 1];
tl2_mpy = [0 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 1];
tl3_mpy = [1 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 0];


tr1_mpy = [0 0 -1; 0 0 0; 1 0 0];
tr2_mpy = [0 0 0 0 0; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 1 0 0 0 0];
tr3_mpy = [0 0 0 0 1; 0 0 0 0 0; 0 0 -1 0 0; 0 0 0 0 0; 0 0 0 0 0];

thv1_mpy = [-1 0 1];
thv2_mpy = [-1 0 0 0 1];

hl_mpy =  [0 0 1; 0 0 0; 1 0 0];
hr_mpy =  [1 0 0; 0 0 0; 0 0 1];
hh_mpy =  [1 0 0 0 1];
% first padding zeros.

image_padding(padding+1:padding+height,padding+1:padding+width) = image;
% get 3x3 block and 5x5 block for left and right
for i=padding+1:padding+height
    for j=padding+1:padding+width
        single_block = image_padding(i-2:i+2,j-2:j+2);
        % G 
        G_left_mode(i-padding,j-padding) =  abs( single_block(2,2)*-1 + single_block(4,4)*1  ) + ...
                                            abs( single_block(3,3)*-1 + single_block(5,5)*1  ) + ...
                                            abs( single_block(3,3)*-1 + single_block(1,1)*1 );
        G_left(i-padding,j-padding) =  uint16( ( single_block(2,2)*1 + single_block(4,4)*1  ));                                
                                        
        G_right_mode(i-padding,j-padding) = abs( single_block(2,4)*-1 + single_block(4,2)*1  ) + ...
                                            abs( single_block(3,3)*-1 + single_block(5,1)*1  ) + ...
                                            abs( single_block(3,3)*-1 + single_block(1,5)*1 );    
                                        
        G_right(i-padding,j-padding) =  uint16(  ( single_block(2,4)*1 + single_block(4,2)*1  )); 
        
        G_choose_r_ = G_left_mode > G_right_mode;


        
        % R B
        G_hori_mode(i-padding,j-padding) =   abs( single_block(3,2)*-1 + single_block(3,4)*1  ) + ...
                                            abs( single_block(3,1)*-1 + single_block(3,5)*1  );
        G_hori(i-padding,j-padding)     =  uint16( ( single_block(3,1)*1 + single_block(3,5)*1  ));                                
                                        
        G_veri_mode(i-padding,j-padding) = abs( single_block(2,3)*-1 + single_block(4,3)*1  ) + ...
                                            abs( single_block(1,3)*-1 + single_block(5,3)*1  );    
                                        
        G_vert(i-padding,j-padding) =  uint16( ( single_block(1,3)*1 + single_block(5,3)*1  )); 
        
        RB_choose_v_ = G_hori_mode > G_veri_mode;
       
    end    
end
% find(G_choose_r~=G_choose_r_)
% find(fl~=G_left)
% find(fr~=G_right)
% 
% find(RB_choose_v~=RB_choose_v_)
% find(fh~=G_hori)
% find(fv~=G_vert)

fg_fushion = int16(G_choose_r_).*int16(G_right) + int16(1-G_choose_r_).*int16(G_left);
frb_fushion = int16(RB_choose_v_).*int16(G_vert) + int16(1-RB_choose_v_).*int16(G_hori);

find(fg ~= fg_fushion) 
find(frb ~= frb_fushion)
else
    fg_fushion = fg; 
	frb_fushion = frb;
end

 
TF = int16(repmat([1 0; 0 1], [size(image,1)/2, size(image,2)/2]));
f_image = fg_fushion.*TF + frb_fushion.*(1-TF);
% f_image 的向右移位操作 在这里不需要做，放到后级
f_image_ = floor(double(f_image)/2);

l_image = l_image + f_image.*int16(1-F);
s_image = s_image*param.exptimes + f_image.*int16(F);
