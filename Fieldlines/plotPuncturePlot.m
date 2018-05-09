function [fig_handle, r, z] = plotPuncturePlot(IOTA_DIRECTORY_PATH, surfaceIndices, indexInc, highlightSurfaces, colorSet, INTERACTIVE)
% function [fig_handle] = plotPuncturePlot(mode, surfaceIndices, indexInc, highlightSurfaces, colorSet, INTERACTIVE);


if (nargin < 6)
    INTERACTIVE = 0;
end
if (nargin < 5)
    colorSet = 'rmkgcy';
end
if (nargin < 4)
    highlightSurfaces = [];
end
if (nargin < 3)
    indexInc = 1;
end

foldername = [IOTA_DIRECTORY_PATH '/' ];

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
        filedata = load([foldername '/' filenames{ii}]);
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
    
jj = 1;
if isempty(highlightSurfaces)
%     highlightSurfaces = 1;
end

for ii = highlightSurfaces
    figure(fig_handle);
    plot(r{ii}, z{ii}, [colorSet(mod(jj,length(colorSet))+1) '.'], 'MarkerSize', 6);
    jj = jj+1;
    hold on;
end
    
surf_num_user = 1; %surfaceIndices(1);

while INTERACTIVE
    surf_user = plot(r{(surf_num_user)}, ...
        z{(surf_num_user)}, 'r.', 'MarkerSize', 8);
    surf_num_user = input('Choose surface.  -1 to leave highlighted.  0 to exit.');
    if isempty(surf_num_user)
        surf_num_user = 1;
    end
    if surf_num_user > length(surfaceIndices);
        surf_num_user = length(surfaceIndices);
        disp('Maximum surface selected');
    end
    if surf_num_user > length(surfaceIndices);
        surf_num_user = length(surfaceIndices);
        disp('Maximum surface selected');
    end
    if surf_num_user == -1
        surf_num_user = 1;
    else
        set(surf_user, 'Visible', 'off');
    end
    if surf_num_user == 0
        INTERACTIVE = 0;
    end
        
end
    