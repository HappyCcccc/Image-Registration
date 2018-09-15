function y = cost2(s,t,r)
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

target1 = imrotate(target1,30);
[t1m,t1n]=size(target1);

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
%{
        cost = 0.0;
        
        for i=1:im
            for j=1:in
                if(Y(i,j)~=0&&target1(i,j)~=0)
                cost = cost + (target1(i,j)-Y(i,j))^2;
                end
            end;
        end;
 %}
W = Y(1:im,1:in);
target2 = target1(1:im,1:in);

cost = sum(sum(    ((W&target2).*(target2-W)).^2   ));
y = cost;
end