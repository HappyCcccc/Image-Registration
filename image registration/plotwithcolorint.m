clear all; close all; clc;


%z =0.0;
% Compute the value of the quadratic function.

d = 2;
n = 1000;
Corner = zeros(1,d);% matrix for corner node
Width = zeros(1,d);% matrix for width and length
Center = zeros(2*n+1,d);%matrix for center node
area = ones(2*n+1,1); % matrix for retrangle area
f = zeros(1,1); % function value for center points
garf = zeros(1,1);% gradents for center points
getrho = zeros(1,1);
a = zeros(1,d);

Width(1,1)=1.0;
Width(1,2)=1.0;
Center(1,1)=0.5;
Center(1,2)=0.5;

maxWidth = 0.0;
maxIndex = 1;
maxD = 0;
MinArea = 1.0;
Center(maxIndex,:);

for i = 1:n
    
maxWidth = 0.0;
% split the retrangle
    %decide which axis to split on
for j = 1:d
        if(Width(maxIndex,j)>maxWidth)
            maxWidth = Width(maxIndex,j);
            maxD = j;
      %      maxIndex = j;
        end
end

%update corner node and width/length information for new retrangles
%NC1 = zeros(d);

Width(maxIndex, maxD)=1/3*maxWidth;

NC1 = zeros(1,d);
NC2 = zeros(1,d);
NW1 = zeros(1,d); 
NW2 = zeros(1,d);  

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
    
    %Center(j, maxIndex) = Corner(j,maxIndex)+ 1/2 * Width(j,maxIndex);
    % area(maxIndex+1) = 1/3 * area(maxIndex);
    % area(maxIndex+2) = 1/3 * area(maxIndex);
    % area = area * Width(j,maxIndex);
end

%update corner and width information
Corner = [Corner;NC1;NC2];
Width = [Width;NW1;NW2];
[m,d] = size(Corner);

%calculate center nodes and retrangle area
for j = 1:m
    for k = 1:d
        Center(j,k)=Corner(j,k)+1/2*Width(j,k);
%        area(j,1)=area(j,1)*Width(j,k);
%        f(j) = funct(Center);
    end
end

%{
%get the smallest function value among all vertices
Mn = intmax('int64');
for j = 1:m
    for k = 1:d
        a(k) = Center(j,k);
    end
        f(j) = funct1(a,d);
        if(f(j)<Mn)
            Mn = f(j);
        end
end
Mn;
%}

%get the smallest estimated function value among all vertices
Wid = zeros(1,1);
Mn = intmax('int64');
for j = 1:m
    for k = 1:d
        a(k) = Center(j,k);
    %    Wid = Wid + Width(1,k);
    end
        f(j) = funct1(a,d);
        
        for k =1
        garf(j)= f(j)- abs(nfunct1(a,d)*Width(1,k)/2+nfunct2(a,d)*Width(1,k+1)/2);
        end
        
        if(garf(j)<Mn)
            Mn = garf(j);
        end
end
Mn;

% get the volume of smallest rectangle:travsal all rectangles
area = ones(2*n+1,1);
for j = 1:m
    for k = 1:d
        area(j,1)= area(j,1)*Width(j,k);
    end
end

for j =1:m
    area(j,1);
end

Vn = intmax('int64');
for j = 1 : m
   % for k = 1 : d
   %     area(j,1)=area(j,1)*Width(j,k)
        if (area(j,1) < Vn)
            Vn = area(j,1);
        end
   % end
end
Vn;

% get GVn
GVn = d*(Vn*log(1/Vn)).^(2/d);

%choose the retrangle which have the largest Rho
rho = intmin;
for j = 1 : m
    for k = 1:d
        a(k) = Center(j,k);
    end
        getrho(j) = area(j,1).^(2/d)/(garf(j)-Mn+d*GVn);
        if (getrho(j) > rho)
            rho = getrho(j);
            maxIndex = j; 
            plot(Center(maxIndex,1),Center(maxIndex,2),'.')
            hold on
            Center(maxIndex,:);
            pause(0.00000000000001);
        end
end


% Center(maxIndex, :);
% plot(Center(:,1), Center(:,2),'.');

end



function y = funct(c,d)
%n = size(c);
sum=0;
for i = 1:d
    sum = sum+ (c(i)-0.5).^2;
end
y = sum;
end

function y = funct1(c,d)
%n = size(c);
nc=zeros(d);
sum=20;
B=4.0;
for i = 1:d
    nc(i)=2*B*c(i)-B-0.1;
    sum = sum+ nc(i).^2-10*cos(2*pi*nc(i));
end
%  y=sum;
y = sum*0.001;
end

function y = nfunct1(c,d)
%n = size(c);
nc=zeros(d);
B=4.0;
sum = 0;
 i = 1;
    nc(i)=2*B*c(i)-B-0.1;
    sum = sum+ 2*nc(i)+ 20*pi*sin(2*pi*nc(i));

  y=sum;
end


function y = nfunct2(c,d)
%n = size(c);
nc=zeros(d);
B=4.0;
sum = 0;
 i = 2;
    nc(i)=2*B*c(i)-B-0.1;
    sum = sum+ 2*nc(i)+ 20*pi*sin(2*pi*nc(i));

  y=sum;
end


function y = funct2(c,d)
%n = size(c);
nc=zeros(d);
sum=0;
B=5.1/(4*pi.^2);
C=5/pi;R=6;T=1/(8*pi);


for i = 1:d-1
    nc(i)=15*c(i)-5;
    nc(i+1)=15*c(i+1);
    
    sum = (1/51.95)*((nc(i+1)-B*nc(i).^2+C*nc(i)-R).^2+10*(1-T)*cos(nc(i))-44.81);
end
y = sum*0.001;
end

