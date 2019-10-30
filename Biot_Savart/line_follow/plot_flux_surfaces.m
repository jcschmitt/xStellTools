function [figureHandleOut, r, z] = plot_flux_surfaces(surfaceIndices, ...
    indexStart, indexInc, indexEnd, figureHandleIn, colorin)

if nargin < 6
    colorin = 'b';
end
if nargin < 5
    figureHandleIn = NaN;
end
if nargin < 4
    indexEnd = NaN;
end
if nargin < 3
    indexInc = 1;
end
if nargin < 1
    indexStart = 1;
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
            filedata.phiIncInDegrees
            %if (filedata.phiIncInDegrees ~= 90)
            %    error('what the heck are you doing?  phiInc ~= 90');
            %end
        end
        if (isnan(indexEnd))
            indexEnd = length(filedata.coords(:,1));
        end
        r{ii} = filedata.coords(indexStart:indexInc:indexEnd,1);
        z{ii} = filedata.coords(indexStart:indexInc:indexEnd,2);
        
    catch
        warning(['Failed to load '  filenames{ii}]);
        r{ii} = [];
        z{ii} = [];
    end
end

if (isnan(figureHandleIn))
    figureHandleOut = figure;
else
    figureHandleOut= figure(figureHandleIn);end

for ii = 1:numFiles
    figure(figureHandleOut);
    plot(r{ii}, z{ii}, 'Color', colorin, 'Marker', '.', 'MarkerSize', 8, ...
        'LineStyle', 'None');
    hold on;
end
figure(figureHandleOut);
ylabel('Z');
xlabel('R');

axis equal

