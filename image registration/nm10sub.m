B= imread('NM10.jpg');
B = rgb2gray(B);


xoffset = 3000;
yoffset = 2500;
dsize = 500;
tsize = 50;
B = B(xoffset:xoffset+dsize,yoffset:yoffset+dsize);
figure
imshow(B)

 hold on
 plot(250,250,'r*');
 plot(300,300,'r*');
 plot(300,250,'r*');
 plot(250,300,'r*');