global B ;
global bm;
global bn;
global tm;
global tn;
global target;
B= imread('NM10.jpg');
B = rgb2gray(B);


T = imread('NM10.jpg');
T = rgb2gray(T);


xoffset = 3000;
yoffset = 2500;
dsize = 500;
tsize = 50;
B = B(xoffset:xoffset+dsize,yoffset:yoffset+dsize);


[bm,bn] = size(B);

tx = dsize/2;

target = T(xoffset+tx:xoffset+tx +tsize,yoffset+tx:yoffset+tx+tsize);
[tm,tn] = size(target);


B = double(B)./256.0;
target = double(target)./256.0;

tic
d = 2;
n = 5000;
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
end

if(f(maxIndex,1)<Mn)
    Mn = f(maxIndex,1);
end


if (area(maxIndex) ~= Vn || f(maxIndex) ~= Mn)
    flag = 1;
else 
    flag = 0;
end

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
    if(NF2<Mn)
        Mn = NF2;
    end
end

if(NA1<Vn)
    Vn = NA1;
    if(NA2<Vn)
        Vn = NA2;
    end
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
GVn = d*(Vn*log(1/Vn)).^(2/d);

if (i==1||flag)
%if (i==1||V(i)~=V(i-1) || M(i)~= M(i-1))
    for j = 1 : m
    % getrho(j) = area(j,1).^(2/d)/(f(j)-Mn+d*GVn);
   
    % getrho(j) = area(j,1).^(2/d)*log(1.01+(f(j)-Mn)/GVn).^(2/d)/(f(j)-Mn+d*GVn);
     getrho(j) = area(j,1).^(2/d)*log(1.01+(f(j)-Mn)/GVn).^(2/d)/(f(j)-Mn+d*GVn);
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

Center(maxIndex, :);
plot(Center(:,1), Center(:,2),'.');
end

toc

function y = funct1(c,d)
global bm;
global tn;

for i = 1:d-1
s = round((bm-tn)*c(i));
t = round((bm-tn)*c(i+1));
end
%y = s+t;
y = funct3(s,t);
end


function y = funct(s,t)
global B;
global bm;
global bn;
global tm;
global tn;
global target;
%f = zeros(1,bm);

%cost = zeros(bm-tm, bn-tn);

        cost = 0.0;
        for i=1:tm
            for j=1:tn
                cost = cost + (target(i,j)-B(s+i,t+j))^2;
            end;
        end;
 
y = cost;
end

function y = funct2(s,t)
global B;
global bm;
global bn;
global tm;
global tn;
global target;

Y = B(s+1:s+tm, t+1:t+tm);
%f = zeros(1,bm);
size(target);
size(Y);
cost = sum(sum((target-Y).*(target-Y)));
                             
y = cost;

end


function y = funct3(s,t)
global B;
global bm;
global bn;
global tm;
global tn;
global target;

Y = B(s+1:s+tm, t+1:t+tm);
%f = zeros(1,bm);
size(target);
size(Y);


%indrow = double(Y(:))+1
%indcol = double(target(:))+1

%[~,~,indrow] = unique(Y(:));
%[~,~,indcol] = unique(target(:));


        C = Y;
T = target;

C_norm = C - min(C(:)) +1; 
T_norm = T - min(T(:)) +1;

[~,~,indrow] = unique(C_norm(:));
[~,~,indcol] = unique(T_norm(:));
h= accumarray([indrow indcol],1);
%{
matAB(:,1) = unique(C_norm(:));
matAB(:,2) = unique(T_norm(:));
h = accumarray(matAB+1, 1); % joint histogram
%}
hn = h./sum(h(:)); % normalized joint histogram
y_marg=sum(hn,1); 
x_marg=sum(hn,2);

Hy = - sum(y_marg.*log2(y_marg + (y_marg == 0))); % Entropy of Y
Hx = - sum(x_marg.*log2(x_marg + (x_marg == 0))); % Entropy of X

arg_xy2 = hn.*(log2(hn+(hn==0)));
h_xy = sum(-arg_xy2(:)); % joint entropy
mutualinformation = (Hx + Hy - h_xy); % mutual information
                             
y = mutualinformation;

end


