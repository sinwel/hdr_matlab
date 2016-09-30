
function GenerateTabandCopyBanks(param);
if nargin < 1
   param.bits           = 10;
   param.blacklevel     = 64;
end

d = 1;
b = 3+d;
c = -2-2*d;

% x = 0:0.01:1;
% figure;plot(b*x.^2+c*x.^3+d*x.^4);


for tidx = 1:961
    D = tidx/961;
    D = b*D.^2+c*D.^3+d*D.^4;
    D_int(tidx) = uint16(D*255);
end
fid = fopen('longshort_mapping_single.dat','w');
for i=1:961
    if i < 961            
        fprintf(fid,'%-4d, ',D_int(i));
    else
        fprintf(fid,'%-4d ',D_int(i));
    end
   	if  0==rem(i,32)
        fprintf(fid,' \n');   
    end
end
%     fprintf(fid,'} ');
% end
fclose(fid);
% creat 16 banks copy.
AAt = [];
ttt = [D_int,0];
t2  = reshape(ttt',2,481)';
for i=1:16    
   AAt = [AAt , t2] ;
end

[m,n] = size(AAt);

fid = fopen('longshort_mapping_16banks.dat','w');
for i=1:m
    for j=1:n              
        fprintf(fid,'%-4d, ',AAt(i,j));
    end
    fprintf(fid,' \n');       
end
%     fprintf(fid,'} ');
% end
fclose(fid);
%
