fid_vecc = fopen('G:\CEVA-XM4_v1.6.1\IMX362HDR\data\raw_2016x1504_Out.raw');
image_vecc = fread(fid_vecc, [2016 1504], 'uint16')';

image_vecc_align = image_vecc;
fclose(fid_vecc);


image_vecc_align_demosic = demosaic(uint16(image_vecc_align), 'grbg');
imwrite(uint16(image_vecc_align_demosic*64), 'G:\HDR\image_vecc_align16bit.png');
figure;imshow(double(demosaic(uint8(image_vecc_align*255/1023), 'grbg'))/255);