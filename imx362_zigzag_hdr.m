function imx362_zigzag_hdr(filename);
if nargin < 1
%     filename = './RK/hdrTest1.raw';
%     filename = 'raw_2016x1504.raw';
    filename_tr = 'raw_2016x1504_new.raw';
end
cfg.width = 2016;
cfg.height = 1504;
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


% fid = fopen(filename);
% image = fread(fid, [cfg.width cfg.height], 'uint16')';
% fclose(fid);
% image = double(image');
% figure;imshow(image,[]);

fid_tr = fopen(filename_tr);
image_tr = fread(fid_tr, [cfg.width cfg.height], 'uint16')';
fclose(fid_tr);

figure;imshow(image_tr,[]);
image  = image_tr;
% imgclip = image(757:757+1503,1009:1009+2015);
% figure;imshow(imgclip,[]);

% fidw = fopen('raw_2016x1504_new.raw','wb');
% fwrite(fidw, uint16(imgclip)', 'uint16')';
% fclose(fidw);
% zigzag_hdr3_float(image, cfg, cfg.exptimes,gcf);
% data_p = zeros(36,132);
% data_p(5:end,3:130) = image(1:32,1:128);
% fidw = fopen('132x36.dat','w');
% for i=1:36
%     for j=1:132
%         fprintf(fidw,'%-4d, ',data_p(i,j));
%     end
%     fprintf(fidw,'\n',data_p(i,j));
% end
% fidw = fopen('128x32.dat','w');
% 
% block_clip = image(1:32,1:128);
% for i=1:32
%     for j=1:128
%         fprintf(fidw,'%-4d, ',block_clip(i,j));
%     end
%     fprintf(fidw,'\n',data_p(i,j));
% end
% fclose(fidw);
zigzag_hdr3(image, cfg, cfg.exptimes);