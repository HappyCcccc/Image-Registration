%scatter3(cos1,cos2,cos,'*')
%surf(cos1,cos2,[cos1 cos2 cos])
%z = [cos1 cos2 cos]

figure
[xi,yi] = meshgrid(0:0.01:1, 0:0.01:1);
zi = griddata(cos1,cos2,cos,xi,yi);
surf(xi,yi,zi);
xlabel('Longitude')
ylabel('Latitude')
zlabel('Depth in Feet')