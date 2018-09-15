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
figure
imshow(B)

[bm,bn] = size(B);

tx = dsize/2;

target = T(xoffset+tx:xoffset+tx +tsize,yoffset+tx:yoffset+tx+tsize);
[tm,tn] = size(target);


B = double(B)./256.0;
target = double(target)./256.0;


tic
d = 2;
n = 65;
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
fileID = fopen('cos.txt','w');
i=0;
%xbest = rand(d,1)
xbest = [0,0]


fmin = funct1(xbest,d);
fprintf(fileID,'x = %f ',xbest(1));
fprintf(fileID,'y = %f ', xbest(2));
fprintf(fileID,'cost = %f\n', fmin);

for i = 0:n
    for j = 0:n
        for k = 0
            x = [i/n,j/n];
             f = funct1(x,d);
            fprintf(fileID,'x = %f ',x(1));
            fprintf(fileID,'y = %f ', x(2));
            fprintf(fileID,'cost = %f\n',f);
           
            if f < fmin
                xbest = x;
                
              %  plot(xbest(1),xbest(2),f)
                hold on
                fmin = f;
            end
        end
    end
end

x = xbest
cost = funct1(x,d)


fclose(fileID);

function y = funct1(c,d)
global bm;
global tn;

for i = 1:d-1
s = round((bm-tn)*c(i));
t = round((bm-tn)*c(i+1));
end
%y = s+t;
y = funct2(s,t);
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
