% deal with target image and floating image

%{
choose a piece of image (250:300,250:300)(name it target, but not real target). 
Then rotated target 30 degree into target1. 
Target1 is the real target we want and we are supposed to calculate the cost 
between target1 and B. I choose same size with target in B and name it Y. 
Rotate Y with degree r into Z. Calculate the cost between target1 and Z.

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
d=3;
fileID = fopen('costofbs.txt','w');
for i=1:100
tic
%x0 = [0.5,0.5,0.5];
x0 = rand(d,1);
fbest = funct(x0,d);
objectivefun = @funct;
lb = [0 0 0];
ub = [1 1 1];

options.FunctionTolerance = 0;
options = optimoptions(@simulannealbnd,'MaxStallIterations',4000,'MaxIterations',4000);
        [x, fval, exitflags, output] = simulannealbnd(objectivefun,x0,lb,ub,options);

toc
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

% s = x(1);
% t = x(2);
% r = x(3);
% Z = imrotate(B,r*180,'bilinear');
% xtrans = s*(bm-tm)
% ytrans = t*(bn-tn)
% rangle = r*180
cost=funct(x)
fprintf(fileID, 'cost = %f\n',cost);
%dlmwrite('costofsa.txt',cost)
end
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

%jmc

s = c(1);
t = c(2);
r = c(3);

y = regimage(s,t,r);
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

