% simulated annealing 2-d

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

fileID = fopen('costofs.txt','w');
for i=1:100
tic
x0 = [0.5,0.5];
objectivefun = @funct1;
lb = [0 0 0];
ub = [1 1 1];

options = optimoptions(@simulannealbnd,'MaxStallIterations',4000,'MaxIterations',4000);
        [x, fval, exitflags, output] = simulannealbnd(objectivefun,x0,lb,ub,options);

toc
b =x

s = x(1);
t = x(2);

xtrans = s*(bm-tm)
ytrans = t*(bn-tn)

cost=funct1(x)
fprintf(fileID, 'cost = %f\n',cost);
%dlmwrite('costofsa.txt',cost)
end
fclose(fileID);


function y = funct1(c)
global bm;
global tn;

for i = 1:2-1
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