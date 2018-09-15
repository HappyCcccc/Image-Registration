% 86 132 133
global B ;
global bm;
global bn;
global tm;
global tn;
global target;
global b4m;
global b4n;
B= imread('NM01.jpg');
B = rgb2gray(B);

[bm,bn]=size(B);
T = imread('NM01.jpg');
T = rgb2gray(T);


xoffset = 100;
yoffset = 150;
dsize = 50;
tsize = 20;
B = B(xoffset:xoffset+dsize,yoffset:yoffset+dsize);
B1 = imrotate(B,30,'bilinear');
B4 = imrotate(B,45,'bilinear');
[b4m,b4n]=size(B4);

[bm1,bn1] = size(B1);
bm1 = round(bm1/2);
bn1 = round(bn1/2);
%target =B1(bm1:bm1+tsize,bn1:bn1+tsize);
tx = dsize/2;


target = T(xoffset+tx:xoffset+tx +tsize,yoffset+tx:yoffset+tx+tsize);
%target = T(150:175,150:170);
[tm,tn] = size(target);


B = double(B)./256.0;
target = double(target)./256.0;

tic
d = 3;
n = 1;
Corner = zeros(1,d);% matrix for corner node
Width = zeros(1,d);% matrix for width and length
Center = zeros(1,d);% matrix for center node
area = ones(1,1);% matrix for retrangle area
f = zeros(1,1);% function value for center points
getrho = zeros(1,1);
a = zeros(1,d);
b = zeros(1,d);
V = zeros(1,n);
M = zeros(1,n);

for i = 1:d
    Width(1,i)=1.0;
    Center(1,i)=0.5;
end

maxWidth = 0.0;
maxIndex = 1;
maxD = 0;

Mn = intmax('int64');
Vn = intmax('int64');
Center(maxIndex,:);


for i = 1:n
    
maxWidth = 0.0;
% split the retrangle
    % decide which axis to split on
for j = 1:d
        if(Width(maxIndex,j)>maxWidth)
            maxWidth = Width(maxIndex,j);
            maxD = j;
        end
end

%update corner node and width/length information for new retrangles
%NC1 = zeros(d);

Width(maxIndex, maxD)=1/3*maxWidth;   
Center = Corner + Width/2;

area(maxIndex)=1;
for j=1:d
    area(maxIndex) = area(maxIndex)*Width(maxIndex, j);
end

a = Center(maxIndex,:);
f(maxIndex,1)= funct1(a,d);


if(area(maxIndex)<Vn)
    Vn = area(maxIndex);
    flag = 1;
end

if(f(maxIndex,1)<Mn)
    Mn = f(maxIndex,1);
    flag = 1;
    b = a;
end

%{
if (area(maxIndex) ~= Vn || f(maxIndex) ~= Mn)
    flag = 1;
else 
    flag = 0;
end
%}
NC1 = zeros(1,d);
NC2 = zeros(1,d);
NW1 = zeros(1,d); 
NW2 = zeros(1,d);
NCE1 = zeros(1,d);
NCE2 = zeros(1,d);
NA1= ones(1,1);
NA2= ones(1,1);
NF1= zeros(1,1);
NF2 = zeros(1,1);

for j = 1:d
    if (j == maxD)
        NC1(j) = 1/3*maxWidth + Corner(maxIndex,j);
        NC2(j) = 2/3*maxWidth + Corner(maxIndex,j);
        NW1(j) = 1/3*maxWidth;
        NW2(j) = 1/3*maxWidth;
    else 
        NC1(j) = Corner(maxIndex,j);
        NC2(j) =  Corner(maxIndex,j);
        NW1(j) = Width(maxIndex,j);
        NW2(j) =  Width(maxIndex,j);
    end
    NA1 = NA1*NW1(1,j);
    NA2 = NA2*NW2(1,j);
end
NCE1=NC1+NW1/2;
NCE2=NC2+NW2/2;
NF1=funct1(NCE1,d);
NF2 = funct1(NCE2,d);

if(NF1<Mn)
    Mn = NF1;
    b = NCE1;
    flag = 1;
end
 
if(NF2<Mn)
    Mn = NF2;
    b = NCE2;
    flag = 1;
end

if(NA1<Vn)
    Vn = NA1;
    flag = 1;
end

if(NA2<Vn)
    Vn = NA2;
    flag = 1;
end


%update corner and width information
Corner = [Corner;NC1;NC2]; 
Width = [Width;NW1;NW2];
Center = [Center;NCE1;NCE2];
area = [area;NA1;NA2];
f = [f;NF1;NF2];

[m,d] = size(Corner);

%Vn = min(area);
% Mn = min(f);

% get GVn
GVn = d*(Vn*log(1/Vn)).^(4/d);

if (i==1||flag)
%if (i==1||V(i)~=V(i-1) || M(i)~= M(i-1))
    for j = 1 : m
    % getrho(j) = area(j,1).^(2/d)/(f(j)-Mn+d*GVn);
   
    % getrho(j) = area(j,1).^(2/d)*log(1.01+(f(j)-Mn)/GVn).^(2/d)/(f(j)-Mn+d*GVn);
     %getrho(j) = area(j,1).^(2/d)*log(1.01+(f(j)-Mn)/GVn).^(2/d)/(f(j)-Mn+GVn);
     getrho(j) = area(j,1).^(2/d)*log(1.01+(f(j)-Mn)/GVn).^(2/d)/(f(j)-Mn+GVn).^(4/d);
    end
    
        [tmp, ind] = sort(getrho,'descend');
        Center = Center(ind,:);
        Corner = Corner(ind,:);
        Width = Width(ind,:);
        f = f(ind);
        area = area(ind);
        maxIndex = 1;  
else
    maxIndex = maxIndex+1;
end
flag = 0;
Center(maxIndex, :);
scatter3(Center(:,1),Center(:,2),Center(:,3),'*');
%plot3(Center(:,1),Center(:,2),Center(:,3));
%plot(Center(:,1), Center(:,2),'.');
end
b
Mn
toc

function y = qufun(c,d)
%n = size(c);
sum=0;
for i = 1:d
    sum = sum+ (c(i)-0.5).^2;
end
y = sum;
end

%interpolators interp2
function y = funct1(c,d)
global bm;
global tn;

%{
for i = 1
s = round((bm-tn)*c(i));
t = round((bm-tn)*c(i+1));
%r = round(360*c(i+2));
r = 360*c(i+2);
end
%}
s = c(1);
t = c(2);
r = c(3);

%y = s+t;
%{
s = 250;
t = 250;
r = 30;
%}
y = regimage(s,t,r);
end


function y = funct(s,t,r)
global B;
global bm;
global bn;
global tm;
global tn;
global target;

target1 = imrotate(target,200,'bilinear');

Y = B(s+1:s+tm, t+1:t+tn);
Y = imrotate(Y,r,'bilinear');
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

function y = regimage(s,t,r)
global B;
global bm;
global bn;
global tm;
global tn;
global target;

target1 = imrotate(target,30,'bilinear');
[t1m,t1n]=size(target1);

Y = B(s+1:s+tm, t+1:t+tn);
%X = B(1:tm, 1:tn);
%Z = imrotate(X,r,'bilinear');
%Y = imtranslate(Z,[s,t],'linear');
Y = imrotate(Y,r,'bilinear');
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
                if(Y(i,j)~=0&&target1(i,j)~=0)
                cost = cost + (target1(i,j)-Y(i,j))^2;
                end
            end;
        end;
 
y = cost;
end

function y = regimage1(s,t,r)
global B;
global bm;
global bn;
global tm;
global tn;
global target;

NB = imrotate(B,30,'bilinear');
[sn1,tn1]=size(NB);
sn1 = round(sn1/2);
tn1 = round(tn1/2);
target1 = NB(sn1+1:sn1+tm, tn1+1:tn1+tn);

Y = B(s+1:s+tm, t+1:t+tn);
Y = imrotate(Y,r,'bilinear');
%Y = imrotate(Y,r,'crop');
[ym,yn]=size(Y);

 cost = 0.0;
        for i=1:tm
            for j=1:tn
             %   if(Y(i,j)~=0)
                cost = cost + (target1(i,j)-Y(i,j))^2;
            %    end
            end;
        end;
        
y = cost;
end

function y = regimage2(s,t,r)
global B;
global bm;
global bn;
global tm;
global tn;
global target;
global b4m;
global b4n;


B2 = imrotate(B,r*360,'bilinear');
%[b2m,b2n]=size(B2);
B3 = imtranslate(B2,[-s*b4m,-t*b4n],'linear');


%f = zeros(1,bm);

%cost = zeros(bm-tm, bn-tn);

        cost = 0.0;
        for i=1:tm
            for j=1:tn 
                cost = cost + (target(i,j)-B3(i,j))^2;
            end;
        end;
 
y = cost;
end

function y = funct2(s,t,r)
global B;
global bm;
global bn;
global tm;
global tn;
global target;

Y = B(s+1:s+tm, t+1:t+tm);
Y = imrotate(Y,r,'crop');



%f = zeros(1,bm);
size(target);
size(Y);
cost = sum(sum((target-Y).*(target-Y)));
                             
y = cost;

end


function y = funct3(s,t,r)
global B;
global bm;
global bn;
global tm;
global tn; 
global target;

Y = B(s+1:s+tm, t+1:t+tm);

Y = imrotate(Y,r,'bilinear','crop');

%f = zeros(1,bm);
size(target);
size(Y);

%indrow = double(Y(:))+1
%indcol = double(target(:))+1

[~,~,indrow] = unique(Y(:));
[~,~,indcol] = unique(target(:));

jointHistogram = accumarray([indrow indcol],1);
jointProb = jointHistogram / numel(indrow);

indNoZero = jointHistogram ~=0;
jointProb1DNoZero = jointProb(indNoZero);

jointEntropy = -sum(jointProb1DNoZero.*log2(jointProb1DNoZero));

histogramImage1 = sum(jointHistogram,1);
histogramImage2 = sum(jointHistogram,2);

indNoZero = histogramImage1 ~= 0;

prob1NoZero = histogramImage1(indNoZero)/numel(histogramImage1);
entropy1 = -sum(prob1NoZero.*log2(prob1NoZero));

indNoZero = histogramImage2 ~=0;
prob2NoZero = histogramImage2(indNoZero)/numel(histogramImage2);
entropy2 = -sum(prob2NoZero.*log2(prob2NoZero));

mutualinformation = -(entropy1+entropy2-jointEntropy);
                             
y = mutualinformation;

end
