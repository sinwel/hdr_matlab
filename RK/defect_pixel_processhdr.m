
function out = defect_pixel_processhdr(I)
X = I;
%P = [0 0 1 0 0; 0 1 0 1 0; 1 0 1 0 1; 0 1 0 1 0; 0 0 1 0 0];
P = [
    0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 1 0 0;
    0 0 0 0 0 0 0 0 0;
    1 0 0 0 1 0 0 0 1;
    0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 1 0 0;
    0 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0;
    ];
c1 = ordfilt2(X, 1, P);
c2 = ordfilt2(X, 2, P);
c8 = ordfilt2(X, 6, P);
c9 = ordfilt2(X, 9, P);
dis = c8-c2;
t1 = (c8-c1) > 1.2*dis;
t2 = (c9-c2) > 1.2*dis;
t1 = t1.*(X == c1);
t2 = t2.*(X == c9);
t = double(t1|t2);
X = X.*(1-t)+t1.*c2+t2.*c8;
Y = X;

X = I;
P = [
    1 0 1 0 1;
    0 0 0 0 0;
    1 0 1 0 1;
    0 0 0 0 0;
    1 0 1 0 1;
    ];
c1 = ordfilt2(X, 1, P);
c2 = ordfilt2(X, 2, P);
c8 = ordfilt2(X, 6, P);
c9 = ordfilt2(X, 9, P);
dis = c8-c2;
t1 = (c8-c1) > 1.2*dis;
t2 = (c9-c2) > 1.2*dis;
t1 = t1.*(X == c1);
t2 = t2.*(X == c9);
t = double(t1|t2);
X = X.*(1-t)+t1.*c2+t2.*c8;

T = repmat([1 0; 0 1], size(I,1)/2, size(I,2)/2);
out = X.*T+Y.*(1-T);


