function [fig_handle, r, z] = plot_flux_surfaces(surfaceIndices, indexInc)

if nargin < 2
    indexInc = 1;
end

surfaceIndices = sort(surfaceIndices);
numFiles = length(surfaceIndices);
for ii = 1:numFiles
    filenames{ii} = ['LineFollow_Coord_Data_surface_' num2str(surfaceIndices(ii))];
    %filenames{ii} = ['HSXLineFollow_Coord_Data_surface_' num2str(surfaceIndices(ii))];
end

% loop through each surface (one surface per file) and load the puncture
% plot data points from each file
for ii = 1:numFiles
    try
        disp(['Loading surface ' num2str(surfaceIndices(ii))]);
        %filedata = load([foldername '/' filenames{ii}]);
        filedata = load([filenames{ii}]);
        if (indexInc == 1)
            if (filedata.phiIncInDegrees ~= 90)
                error('what the heck are you doing?  phiInc ~= 90');
            end
        end
        r{ii} = filedata.coords(1:indexInc:end,1);
        z{ii} = filedata.coords(1:indexInc:end,2);
        
    catch
        warning(['Failed to load '  filenames{ii}]);
        r{ii} = [];
        z{ii} = [];
    end
end

fig_handle = figure;
for ii = 1:numFiles
    figure(fig_handle);
    plot(r{ii}, z{ii}, '.', 'MarkerSize', 2);
    hold on;
end
figure(fig_handle);
ylabel('Z');
xlabel('R');
    
axis equal
    
    