function y = funct(s,t,r)
global B;
global bm;
global bn;
global tm;
global tn;
global target;

Y = B(s+1:s+tm, t+1:t+tm);
Y = imrotate(Y,r);
%Y = imrotate(Y,r,'crop');

size(target);
size(B);
[ym,yn]=size(Y);
target1=B(1:ym,1:yn);
size(target1);

%f = zeros(1,bm);

%cost = zeros(bm-tm, bn-tn);

        cost = 0.0;
        for i=1:tm
            for j=1:tn
                if(Y(i,j)~=0)
                cost = cost + (target1(i,j)-Y(i,j))^2;
                end
            end;
        end;
 
y = cost;
end