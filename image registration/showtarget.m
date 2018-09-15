 B = imread('dog.jpg');
 B = rgb2gray(B);
 B = imrotate(B,90,'bilinear');
 target = B(300:350,300:350);
 imshow(target);
%  imshow(B);
%  
%  hold on
%  plot(300,300,'r*');
%  plot(350,350,'r*');
%  plot(300,350,'r*');
%  plot(350,300,'r*');
 
 %RGB = insertMarker(B,[325 325]);
 %pos = [300 300;300 350;350 300;350 350];
% 
 %color = {'red','white','green','magenta'};
 %RGB = insertMarker(RGB, pos,'x','color',color,'size',10);
% 
