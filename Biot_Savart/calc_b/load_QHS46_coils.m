function QHS46_coils = load_QHS46_coils()
% function QHS46_coils = load_QHS46_coils()
%
% Loads the information for the details of the QHS46 field coils into the
% QHS46_coils structure.
%
% JC Schmitt, 2018

% The variable 'in_foldername' is a string that points to the location
% where this Biot-Savart code is located.  For instance, the following line
% works on my Macbook:
% in_foldername = '/Users/jschmitt/Documents/MATLAB/BiotSavart/';

% Initially, in_foldername will be empty, forcing you to fix it.
%in_foldername = '';
in_foldername = '';

in_filename = 'QHS46_34742_v6a.mat';   

in_filename_complete = [in_foldername in_filename];

disp(['Loading QHS46 coils from: ' in_filename_complete]);

QHS46_coils = load(in_filename_complete);
