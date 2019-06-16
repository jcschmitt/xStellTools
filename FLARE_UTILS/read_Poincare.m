function [R, Z, T, PolFlux] = read_Poincare(filename)

% Verfiy input
if (nargin == 0)
    disp('K==== read_Poincare requires filename.');
    return
end

fid=fopen(filename, 'r');

% First pass gets matrix size
disp('<---- Getting Matrix Size.');

while ~feof(fid);
    line = fgetl(fid);
    if ~(strncmpi(strtrim(line), '#', 1))
        num_lines = 1;
        break;
    end
end

while ~feof(fid);
    line = fgetl(fid);
    num_lines =num_lines + 1;
end
fclose(fid);

R = zeros(1, num_lines);
Z = zeros(1, num_lines);
T = zeros(1, num_lines);
PolFlux = zeros(1, num_lines);

% Second pass gets the data
disp(' - Reading Poincare Data');
fid = fopen(filename,'r');
while ~feof(fid);
    line = fgetl(fid);
    if ~(strncmpi(strtrim(line), '#', 1))
        break;
    end
end

i = 1;

while ~feof(fid)
    vals_line = sscanf(line, '%e', 4)';
    R(i) = vals_line(1)/100; % cm -> m
    Z(i) = vals_line(2)/100; % cm ->
    T(i) = vals_line(3)*pi/180; % degrees -> radians
    PolFlux(i) = vals_line(4);  % unchecked 
    line = fgetl(fid);
    i = i + 1;
end
fclose(fid);




