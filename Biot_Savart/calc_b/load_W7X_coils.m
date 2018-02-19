function W7X_coils = load_W7X_coils()
% function LTX_coils = load_W7X_coils()
%
% Loads the information for the details of the W7X field coils into the
% W7X_coils structure.
%
% JC Schmitt, 2017

% The variable 'in_foldername' is a string that points to the location
% where this Biot-Savart code is located.  For instance, the following line
% works on my Macbook:
% in_foldername = '/Users/jschmitt/Documents/MATLAB/BiotSavart/';

% Initially, in_foldername will be empty, forcing you to fix it.
%in_foldername = '';
in_foldername = '';

in_filename = 'W7X_coils_v1.mat';   

in_filename_complete = [in_foldername in_filename];

disp(['Loading W7X coils from: ' in_filename_complete]);

W7X_coils = load(in_filename_complete);
