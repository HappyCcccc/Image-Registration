global B ;
global bm;
global bn;
global tm;
global tn;
global target;
global Mn;
global optimage;



%A = imread('dog2.png');
%A = rgb2gray(A);
%figure
%imshow(A);
%A = double(A)./256.0;
B = imread('dog.jpg');
B = rgb2gray(B);
figure
imshow(B);

B = double(B)./256.0;
[bm,bn] = size(B);

BR = imrotate(B,90,'bilinear');
[bm1,bn1]=size(BR);

target = BR(300:350,300:350);
%target = BR(100:150,150:200);
[tm,tn] = size(target);

Z3 = imtranslate(BR,[-299,-199]);
%imshowpair(target, Z3(1:tm,1:tn), 'montage')
fileID = fopen('costofexhdog.txt','w');
%pause(50)
for loop=1:100
tic
d = 3;
n = 16;
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
u = -0.03125+0.00625*rand(1,1);
i=0;
xbest = rand(d,1);
%xbest = [0,0,0]
fmin = funct(xbest,d);
n= 16;
for i = 1:n
    for j = 1:n
        for k = 1:n
            x = [i/n+u,j/n+u,k/n+u];
            f = funct(x,d);
            
            if f < fmin
                xbest = x;
                fmin = f;
            end
        end
    end
end

x = xbest;


b =x

% termination tolerance
tol = 1e-6;

% maximum number of allowed iterations
maxiter = 5000;

% minimum allowed perturbation
dxmin = 1e-6;

% step size ( 0.33 causes instability, 0.2 quite accurate)
alpha = 0.001;

% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf; x = b'; niter = 0; dx = inf;

% gradient descent algorithm:
while and(gnorm>=tol, and(niter <= maxiter, dx >= dxmin))
    % calculate gradient:
    g = gradient(x);
    gnorm = norm(g);
    % take step:
    xnew = x - alpha*g;
    newalpha = alpha;
    
    while (any(xnew(:)<0|xnew(:)>1))|(funct(xnew)>funct(x))
       newalpha = 0.5*newalpha;
       xnew = x -newalpha*g;
    end
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    % plot current point
    %plot([x(1) xnew(1)],[x(2) xnew(2)],'ko-')
    refresh
    % update termination metrics
    niter = niter + 1
    dx = norm(xnew-x);
    x = xnew;
    
end

xopt = x



toc
s = x(1);
t = x(2);
r = x(3);
Z = imrotate(B,r*180,'bilinear');
xtrans = s*(bm-tm)
ytrans = t*(bn-tn)
rangle = r*180
cost=funct(x,d)

fprintf(fileID, '%f\n',cost);

end;

fclose(fileID);
%Z1 = imtranslate(Z,[-s*(bm-tm),-t*(bn-tn)]);
%Z2 = Z1(1:tm,1:tn);
%imshowpair(target, Z2, 'montage')
%pause(5)
%Z = imrotate(B,30,'bilinear');
%Z3 = imtranslate(BR,[-299,-199]);
%imshowpair(target, Z3(1:tm,1:tn), 'montage')

%regimage2(299.0/(bm-tm),299.0/(bn-tn),30.0/180.0)
%regimage2(b(1),b(2),b(3))

function y = funct(c,d)
global bm;
global tm;

%jmc

s = c(1);
t = c(2);
r = c(3);

%{
s = 250;
t = 250;
r = 30;
%}

y = regimage(s,t,r);
end


function y = gradient(c)
global bm;
global tm;
global bn;
global tn;

%jmc

s = c(1);
t = c(2);
r = c(3);

%{
s = 250;
t = 250;
r = 30;
%}

y = [(regimage(s+1/(bm-tm),t,r)-regimage(s,t,r))/(1/(bm-tm))
     (regimage(s,t+1/(bn-tn),r)-regimage(s,t,r))/(1/(bn-tn))
     (regimage(s,t,r+1/180)-regimage(s,t,r))/(1/(1/180))];
end

function y = regimage(s,t,r)

global B;
global target;
global tm;
global tn;
%jmc
global bm;
global bn;
global Mn;
global optimage;

Z = imrotate(B,r*180,'bilinear');
%Z1 = imtranslate(Z,[-s*(bm-tm),-t*(bn-tn)]);
Z1 = imtranslate(Z,[-s*(bm-tm),-t*(bn-tn)]);
Z2 = Z1(1:tm,1:tn);


[cost,nn] = sumsqr(target-Z2);
cost = sqrt(cost);

if(cost<Mn)
   optimage = Z2; 
end
%size(target)
%size(Z2)

y = cost;

end
function y = regimage2(s,t,r)

global B;
global target;
global tm;
global tn;
%jmc
global bm;
global bn;
global optimage;

s*(bm-tm)
t*(bn-tn)
r*180


Z = imrotate(B,r*180,'bilinear');
%Z2=Z(1+s*(bm-tm):tm+s*(bm-tm),1+t*(bn-tn):tn+t*(bn-tn));
Z1 = imtranslate(Z,[-s*(bm-tm),-t*(bn-tn)]);
Z2 = Z1(1:tm,1:tn);

%{
Z = imrotate(B,30,'bilinear');
Z1 = imtranslate(Z,[-299,-199]);
Z2 = Z1(1:tm,1:tn);
%}


imshowpair(target, optimage, 'montage')

[cost,nn] = sumsqr(target-Z2);
cost = sqrt(cost);

diff = target -Z2;
%diff
%target(1:5,1:5)
%Z2(1:5,1:5)

%size(target)
%size(Z2)

y = cost;

end

