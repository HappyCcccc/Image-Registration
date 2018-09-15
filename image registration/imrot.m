I = imread('circuit.tif');
J = imrotate(I,30,'bilinear');

K = imrotate(J,-30,'bilinear');
L = imrotate(I, 35, 'loose','bilinear');

figure
imshow(I);
figure
imshow(J);
figure
imshow(K);

figure
%{
imshow(J)
pause(3);
imshow(K)
pause(3);
imshow(L)
%}
imshowpair(I,K,'montage')
%imshowpair(I,J,'montage')
axis off


