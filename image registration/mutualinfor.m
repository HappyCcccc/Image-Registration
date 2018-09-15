B= imread('NM10.jpg');
B = rgb2gray(B);


T = imread('NM10.jpg');
T = rgb2gray(T);


xoffset = 3000;
yoffset = 2500;
dsize = 500;
tsize = 50;
%B = B(xoffset:xoffset+dsize,yoffset:yoffset+dsize);
B = B(xoffset:xoffset+dsize,yoffset:yoffset+dsize);
T = T(xoffset:xoffset+tsize,yoffset:yoffset+tsize);

[bm,bn] = size(B)
[tm,tn] = size(T)
bm-tm
bn-tn

for s=1:bm-tm
    for t=1:bn-tn
        C = B(s:tsize+s,t:tsize+t);
        
        %{
        indrow = double(C(:))+1;
        indcol = double(T(:))+1;

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

        mutualInformation(s,t) = -(entropy1+entropy2-jointEntropy);
        %}
        C = double(C);
T = double(T);

C_norm = C - min(C(:)) +1; 
T_norm = T - min(T(:)) +1;

matAB(:,1) = C_norm(:);
matAB(:,2) = T_norm(:);
h = accumarray(matAB+1, 1); % joint histogram

hn = h./sum(h(:)); % normalized joint histogram
y_marg=sum(hn,1); 
x_marg=sum(hn,2);

Hy = - sum(y_marg.*log2(y_marg + (y_marg == 0))); % Entropy of Y
Hx = - sum(x_marg.*log2(x_marg + (x_marg == 0))); % Entropy of X

arg_xy2 = hn.*(log2(hn+(hn==0)));
h_xy = sum(-arg_xy2(:)); % joint entropy
M(s,t) = -(Hx + Hy - h_xy); % mutual information
        
    end
end
cc = 20;

scost = M(1:cc:bn-tn,1:cc:bm-tm);

[X,Y] = meshgrid(1:cc:bn-tn, 1:cc:bm-tm);

[XI,YI] = meshgrid(1:1:bn-tn, 1:1:bm-tm);
%ZI = interp2(X,Y,scost,XI,YI,'spline');
ZI = interp2(X,Y,scost,XI,YI);

subplot(1,1,1);
mesh(XI,YI,ZI)



%B = double(B)./256.0;
%T = double(T)./256.0;























