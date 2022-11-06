function plot_contours(inputfile)

clc
close all

n = 50;

fileID    = H5F.open(inputfile,'H5F_ACC_RDONLY','H5P_DEFAULT');
datasetID = H5D.open(fileID,'/nodes','H5P_DEFAULT');

nodes    = H5D.read(datasetID,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT');
x_coords = nodes(1,:)';
x_coords = x_coords(1:n:end);
y_coords = nodes(2,:)';
y_coords = y_coords(1:n:end);

datasetID = H5D.open(fileID,'/v','H5P_DEFAULT');
v         = H5D.read(datasetID,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT');
v         = v';
v_magn    = sqrt(v(:,1).^2 + v(:,2).^2);

H5D.close(datasetID);
H5F.close(fileID);

[X, Y] = meshgrid(unique(x_coords), unique(y_coords));

Z = griddata(x_coords, y_coords, v_magn(1:n:end), X, Y, 'natural');
levels_v = [1.1 : 0.1 : 2];
contourf(X,Y,Z, levels_v)
xlim([9.5 10.5])
ylim([9.5 10.5])
% xlabel('x [m]','Fontsize',24)
% ylabel('z [m]','Fontsize',24)
% ax = gca;
% ax.XAxis.FontSize = 20;
% ax.YAxis.FontSize = 20;

end