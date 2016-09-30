clear
fid_vecc = fopen('G:\CEVA-XM4_v1.6.1\IMX362HDR\data\raw_2016x1504_Out.raw');
image_vecc = fread(fid_vecc, [2016 1504], 'uint16')';
fclose(fid_vecc);
param.exptimes = 8;
% thumb for WDR
scale_gap = 16;
for m=1:scale_gap:1504
    for n=1:scale_gap:2016
        out_blk = image_vecc(m:m+scale_gap-1,n:n+scale_gap-1);
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
imwrite(uint16(image_vecc_align_demosic*64), 'image_vecc_align16bit.png');
