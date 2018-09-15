B = imread('NM10.jpg');
B = rgb2gray(B);
size(B);
%imshow(B)

%pause(10)

%A = imread('NC1NHAP020407050.tif');
%B = rgb2gray(A);
xoffset = 3000;
yoffset = 2500;
dsize = 500;
tsize = 50;
B = B(xoffset:xoffset+dsize,yoffset:yoffset+dsize);

%imwrite(B,'map.jpg');

[bm,bn] = size(B)
%imshow(B)
%pause(10)


T = imread('NM10.jpg');
T = rgb2gray(T);
tx = dsize/2;

target = T(xoffset+tx:xoffset+tx +tsize,yoffset+tx:yoffset+tx+tsize);
[tm,tn] = size(target);

%subplot(1,2,1);
%imshow(B, 'InitialMagnification', 100)
%subplot(1,2,2);
%imshow(target, 'InitialMagnification', 100)
%pause(10)


B = double(B)./256.0;
target = double(target)./256.0;



%imshow(target)

f = zeros(1,bm);
cost = zeros(bm-tm, bn-tn);
for s=1:bm-tm
    for t=1:bn-tn
      %  s
        cost(s,t) = 0.0;
        for i=1:tm
            for j=1:tn
                cost(s,t) = cost(s,t) + (target(i,j)-B(s+i,t+j))^2;
            end;
        end;
    end;
end;
cost(25,25);
cc = 20;
[X,Y] = meshgrid(1:cc:bn-tn, 1:cc:bm-tm);
%[X,Y] = meshgrid(size(cost));

size(cost)

size([X,Y])

for i=1:length(X)
    for j=1:length(Y)
        %Gamma(i,j) = parab([X(1,i) Y(j,1)]);
        %Gamma(i,j) = var(SigmaInv, pointsp,  [X(1,i) Y(j,1)]);
        %Gamma(i,j) = mu(SigmaInv, pointsp,vals,  [X(1,i) Y(j,1)]);
    end;
end;
colormap hot
%mesh(X,Y,cost)

%[X,Y] = meshgrid(-3:.25:3);
%Z = peaks(X,Y);
scost = cost(1:cc:bn-tn,1:cc:bm-tm);

[XI,YI] = meshgrid(1:1:bn-tn, 1:1:bm-tm);
%ZI = interp2(X,Y,scost,XI,YI,'spline');
ZI = interp2(X,Y,scost,XI,YI);
%mesh(X,Y,cost), hold, 
subplot(1,1,1);
mesh(XI,YI,ZI);

%contour(XI,YI,ZI)

%hold off
%axis([-3 3 -3 3 -5 20])