clear all;
% 1. Establish bounds for variables
bounds = [-2 2;-2 2];

% 2. Send options to Direct
%    We tell DIRECT that the globalmin = 3
%    It will stop within 0.01% of solution
options.testflag  = 1; options.globalmin = 3; options.showits   = 1;
options.tol       = 0.01;

% 2a. NEW!
% Pass Function as part of a Matlab Structure

%by cc
%Problem.f = 'gp';
Problem.f = 'gp';

% 3. Call DIRECT
[fmin,xmin,hist] = direct(Problem,bounds,options);

% 4. Plot iteration statistics
plot(hist(:,2),hist(:,3))
xlabel('Fcn Evals');
ylabel('f_{min}');
title('Iteration Statistics for GP test Function');