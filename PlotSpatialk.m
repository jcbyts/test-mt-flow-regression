function [gridx, gridy, dx, dy, h] = PlotSpatialk(spatial_kernel, nyx, samplingrate, center, arrowcolor)
% Usage: PlotSpatialk(spatial_kernel, lens, <size>, <samplingrate>, <type>)
% Function: plot 2D Spatial Kernel
% Input:
%       spatial_kernel(1:2*lens^2): spatial_kernel(1:lens^2) is
%       the x component, spatial_kernel(lens:2*lens^2) is the y
%       component
%       lens: # of point at each dimension
%       samplingratex, samplingratey: samplingrate sample/pixel at 2 dim
%

if nargin<5
    arrowcolor='k';
    if nargin < 4
        center = [0 0];
        if nargin < 3
            samplingrate = 2;
            if nargin < 2
                if mod(numel(map(1).staSK)/2, 2)==0
                    nyx = [numel(map(1).staSK)/2;numel(map(1).staSK)/2];
                else
                    disp('you need to pass in the dimensions of your image')
                end
            end
        end
    end
end

if numel(nyx)==1
    nyx = ones(2,1)*nyx;
end

if numel(samplingrate)==1
    samplingrate = ones(2,1)*samplingrate;
end


spacefilt=reshape(spatial_kernel, length(spatial_kernel), 1)/norm(spatial_kernel);

x = linspace(-nyx(2)/2,nyx(2)/2, nyx(2))*samplingrate(1);
y = linspace(-nyx(1)/2,nyx(1)/2, nyx(1))*samplingrate(2);
xax = x + center(1);
yax = y + center(2);

[gridx, gridy] = meshgrid(xax,yax);


dx = spacefilt(1:end/2);
dy = spacefilt(end/2+1:end);

h = quiver(gridx(:),gridy(:),dx,dy,arrowcolor);
xlim([xax(1) xax(end)])
ylim([yax(1) yax(end)])
set(gca, 'YDir', 'normal');