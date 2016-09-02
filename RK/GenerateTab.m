
function GenerateTab(param);
if nargin < 1
   param.bits           = 10;
   param.blacklevel     = 64;
end
tabidx = 0;
fid = fopen('G:\HDR\RK\tone_mapping_961.dat','w');
% for exptimes = 1:16
    tabidx = tabidx + 1;
    %bits = 10 + ceil( log2(exptimes)); % ×î´ó14bit
%     RefValue = 2^bits;
    % 
    RefValue = 2^param.bits - param.blacklevel; % 64 is blacklevel.
    precison = 1/(RefValue);
    idx = 0;
%     fprintf(fid,'{ ');
    for i=0:precison:1
        d = 1;
        b = 3+d;
        c = -2-2*d;
        % loop up table
        % x = 0:0.01:1;
        % figure;plot(b*x.^2+c*x.^3+d*x.^4);
        D_float = i;
        idx = idx + 1;
        toneTab(tabidx,idx) = int16((b*D_float.^2+c*D_float.^3+d*D_float.^4)*RefValue);
        
        if i==1-precison
            fprintf(fid,'%-4d ',toneTab(tabidx,idx));
        else
            fprintf(fid,'%-4d, ',toneTab(tabidx,idx));
        end
        if mod(idx,64)==0
            fprintf(fid,' \n');
        end
    end
%     fprintf(fid,'} ');
% end
fclose(fid);
%
