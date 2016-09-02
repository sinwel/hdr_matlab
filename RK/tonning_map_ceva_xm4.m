function out = tonning_map_ceva_xm4(D,l_image, s_image, param, Max_datain,up_factor, down_factor,times)

Qvalue = Max_datain;
D = min(D,Qvalue); % 限幅,  D 是短曝光时间 * times 可能会大于 最大值 Max_datain.
% D = D*(2^param.bits/param.noise/param.exptimes);
% D = min(D,1); % 限幅 

d = 1;
b = 3+d;
c = -2-2*d;
% loop up table
% x = 0:0.01:1;
% figure;plot(b*x.^2+c*x.^3+d*x.^4);
if 0
    
else
    D_float = double(D)/double(Qvalue);
    D = int16((b*D_float.^2+c*D_float.^3+d*D_float.^4)*Qvalue);
    % D = min(D,1);
    D = min(D,Qvalue); % 限幅,  D 是短曝光时间 * times 可能会大于 最大值 Qvalue.
    D = max(D,0);
end


D_bak = D;
% D = max(D, double(s_image>0.9)); % 取大于0.9倍s_image做为
% D = min(D, 1-double(s_image<0.85/times));% s_image 值太小就截断为0
Upbound = int16(up_factor*double(Qvalue));
% D = max(D, int16( Qvalue*( s_image > Upbound  ) ));
Lowbound = int16(down_factor*double(Qvalue)/times);
% D = min(D, int16( Qvalue*( 1 - (s_image < Lowbound) ) ) );

IDX_UP = s_image > Upbound;
D_bak(IDX_UP) = Qvalue;
IDX_DOWN = s_image < Lowbound ;
D_bak(IDX_DOWN) = 0;

% t=fspecial('gaussian',[9 9], 2);
% t = ones(3,3)/9;
t = [1 0 1 0 1; 0 0 0 0 0; 1 0 1 0 1; 0 0 0 0 0; 1 0 1 0 1]/9;
% t = [1 0 1 0 1; 0 0 0 0 0; 1 0 1 0 1; 0 0 0 0 0; 1 0 1 0 1];
D = int16(conv2(double(D), t, 'same'));

% figure;imshow(D);
out = (int32(l_image).*int32(Qvalue-D)+int32(s_image).*int32(D)) / int32(Qvalue) ;
