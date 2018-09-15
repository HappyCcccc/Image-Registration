function y = regimage1(s,t,r)
global B;
global bm;
global bn;
global tm;
global tn;
global target;

NB = imrotate(B,30,'bilinear');
[sn1,tn1]=size(NB);
target1 = NB(sn1/2+1:sn1/2+tm, tn1/2+1:tn1/2+tn);

Y = B(s+1:s+tm, t+1:t+tn);
Y = imrotate(Y,r,'bilinear');
%Y = imrotate(Y,r,'crop');
[ym,yn]=size(Y);

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