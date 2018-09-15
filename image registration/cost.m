function y = cost(s,t,r)
global B;
global bm;
global bn;
global tm;
global tn;
global target;

target1 = imrotate(target,30)
[t1m,t1n]=size(target1);

Y = B(s+1:s+tm, t+1:t+tn);
Y = imrotate(Y,r);
%Y = imrotate(Y,r,'crop');
[ym,yn]=size(Y);



if(ym>t1m)
    im=t1m;
else im = ym;
end

if(yn>t1n)
    in=t1n;
else in = yn;
end
%f = zeros(1,bm);

%cost = zeros(bm-tm, bn-tn);

        cost = 0.0;
        for i=1:im
            for j=1:in
              %  if(Y(i,j)~=0&&target1(i,j)~=0)
                cost = cost + (target1(i,j)-Y(i,j))^2;
             %   end
            end;
        end;
 
y = cost;
end