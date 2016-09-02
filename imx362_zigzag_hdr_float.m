function imx362_zigzag_hdr_float(filename);
if nargin < 1
%     filename = './RK/hdrTest1.raw';
    filename = 'raw_2016x1504.raw';
end
cfg.width = 1504;
cfg.height = 2016;
% cfg.width = 4032;
% cfg.height = 3024;

cfg.blacklevel  = 64;
cfg.exptimes    = 8;
cfg.bits        = 10;
cfg.noise       = 64;
cfg.wdrgain     = 0;

cfg.zigzagpattern = [1 0 1 1;
                     0 0 1 0;
                     1 1 1 0;
                     1 0 0 0;];


fid = fopen(filename);
image = fread(fid, [cfg.width cfg.height], 'uint16')';
fclose(fid);


image = double(image');
figure;imshow(image,[]);
% imgclip = image(757:757+1053,1009:1009+2015);
% fidw = fopen('raw_2016x1504.raw','wb');
% fwrite(fidw, uint16(imgclip), 'uint16')';
% fclose(fidw);
zigzag_hdr3_float(image, cfg, cfg.exptimes,gcf);
% zigzag_hdr3(image, cfg, cfg.exptimes);