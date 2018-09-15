% deal with target image and floating image

%{
compare fullimagecompare.m using gradient with a random start point
%}

global C;
global B;
global bm;
global bn;
global tm;
global tn;
global target;

C = imread('Brain.png');
%C=imread('rs-trail1.png');
%C=imread('a.png');
%C=imread('bark.000.tiff');
C = double(C)./256.0;
%C =C(1:300,1:300);
%jmc

B = imread('Brain1.png');
%B = imread('rs-trail2.png');
%B = imread('b.png');
%B = imread('bark.030.tiff');
B = double(B)./256.0;
%B = B(1:300,1:300);
[bm,bn] = size(C);

%BR = imrotate(B,30,'bilinear');

%target = BR(100:150,150:200);
%[tm,tn] = size(target);
[tm,tn] = size(target);

%Z3 = imtranslate(BR,[-299,-199]);
%imshowpair(target, Z3(1:tm,1:tn), 'montage')
%pause(50)

tic
d = 3;
n = 200;
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
fileID = fopen('costofbg.txt','w');
for j = 1:100
i=0;


while i<2*n+1
xbest = rand(d,1);
fbest = funct(xbest,d);
% termination tolerance
tol = 1e-6;

% maximum number of allowed iterations
maxiter = 5000;

% minimum allowed perturbation
dxmin = 1e-6;

% step size ( 0.33 causes instability, 0.2 quite accurate)
alpha = 0.001;

% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf; x = rand(d,1); niter = 0; dx = inf;

% gradient descent algorithm:
while and(gnorm>=tol, and(niter <= maxiter, dx >= dxmin))
    % calculate gradient:
    g = gradient(x);
    gnorm = norm(g);
    % take step:
    xnew = x - alpha*g;
    newalpha = alpha;
    
    while (any(xnew(:)<0|xnew(:)>1))|(funct(xnew,d)>funct(x,d))
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
    niter = niter + 1;
    dx = norm(xnew-x);
    x = xnew;
    
end

xopt = x;
fb = funct(xopt,d);
if(fb<fbest)
    fbest = fb;
    xbest = xopt;
end
i = i+niter;

end
x = xbest;
toc
% s = x(1);
% t = x(2);
% r = x(3);
% %Z = imrotate(B,r*180,'bilinear');
% xtrans = s*(bm-tm)
% ytrans = t*(bn-tn)
% rangle = r*180
 cost=funct(x,d)
%regimage2(x(1),x(2),x(3))
fprintf(fileID, 'cost = %f\n',cost);
%dlmwrite('costofsa.txt',cost)
end
fclose(fileID);

function y = funct(c,d)
%jmc
s = c(1);
t = c(2);
r = c(3);
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

global C;
global B;
global target;
global tm;
global tn;
%jmc
global bm;
global bn;

Z = imrotate(B,r*180,'bilinear');
[zm,zn]=size(Z);
Z1 = imtranslate(Z,[-s*0.5*zm,-t*0.5*zn]);
%Z2 = Z1(1:bm,1:bn);

[cm,cn]=size(C);
%[cost,nn] = sumsqr(C-Z1);
[cost,nn] = sumsqr(C(1:min(cm,zm),1:min(cn,zn))-Z1(1:min(cm,zm),1:min(cn,zn)));
cost = cost/(min(cm,zm)*min(cn,zn));
cost = sqrt(cost);

%size(target)
%size(Z2)

y = cost;

end
function y = regimage2(s,t,r)
global C;
global B;

global tm;
global tn;
%jmc
global bm;
global bn;

s*bm*0.5
t*bn*0.5
r*180


Z = imrotate(B,r*180,'bilinear');
[zm,zn]=size(Z);
Z1 = imtranslate(Z,[-s*0.5*zm,-t*0.5*zn]);
%Z2 = Z1(1:bm,1:bn);

[cm,cn]=size(C);
%[cost,nn] = sumsqr(C-Z1);
[cost,nn] = sumsqr(C(1:min(cm,zm),1:min(cn,zn))-Z1(1:min(cm,zm),1:min(cn,zn)));
cost = cost/(min(cm,zm)*min(cn,zn));
cost = sqrt(cost);


%{
Z = imrotate(B,30,'bilinear');
Z1 = imtranslate(Z,[-299,-199]);
Z2 = Z1(1:tm,1:tn);
%}


%imshowpair(C, Z1, 'montage')
figure
imshowpair(C, Z1);

figure
imshowpair(C, Z1, 'montage');




%diff
%target(1:5,1:5)
%Z2(1:5,1:5)

%size(target)?
%size(Z2)

y = cost;

end

