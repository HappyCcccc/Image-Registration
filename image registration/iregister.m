% read two graphs
A = imread('Brain.png');
%A = rgb2gray(A);
A = double(A)./256.0;
B = imread('Brain1.png');
%B = rgb2gray(B);
B = double(B)./256.0;

%A = A(1:300,1:300);
%B = B(1:300,1:300);
%fixed = dicomread('knee1.dcm');
%moving = dicomread('knee2.dcm');
fixed = A;
moving =B;
figure
imshow(fixed);
figure
imshow(moving);

%{
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

fixed = B;
moving = target;
%}

% view the misaligned images
%imshowpair (fixed, moving, 'Scaling','joint');

% create the optimizer and metric
%[optimizer, metric] = imregconfig('multimodal');
[optimizer, metric] = imregconfig('monomodal');
%[optimizer, metric] = imregconfig();

%tune the properties of the optimizer
%{
optimizer.InitialRadius = 0.009;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;
%}


optimizer = registration.optimizer.RegularStepGradientDescent();
%optimizer
metric = registration.metric.MeanSquares();
optimizer.MaximumIterations = 3000;
optimizer.MinimumStepLength = 5e-30;

%perform the registration

movingRegistered = imregister(moving, fixed, 'affine',optimizer, metric,'DisplayOptimization',true);
%movingRegistered = imregister(moving, fixed, 'affine', optimizer, metric);
% view the registerd images
figure
imshowpair(fixed, movingRegistered)

figure
imshowpair(fixed, movingRegistered,'montage')
%imshowpair (fixed, movingRegistered,'Scaling','joint')
%MI_GG(fixed, moving);

% how do we get MI by two images input
% after get MI, How does the algorithm works
% what vaule you need to calculate to get next step 

% f(0) f(1) f(1/2) after n evulations

% what's the domain you have, and how do you divide the
% rectangle into three parts


% define G(n)

%calculate Rho(n)






