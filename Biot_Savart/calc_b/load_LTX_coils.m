function LTX_coils = load_LTX_coils()
% function LTX_coils = load_LTX_coils()
%
% Loads the information for the details of the LTX field coils into the
% LTX_coils structure.
%
% JC Schmitt, 2012,2013

% The variable 'in_foldername' is a string that points to the location
% where this Biot-Savart code is located.  For instance, the following line
% works on my Macbook:
% in_foldername = '/Users/jschmitt/Documents/MATLAB/BiotSavart_LTX/';

% Initially, in_foldername will be empty, forcing you to fix it.
%in_foldername = '';
in_foldername = '';

in_filename = 'LTX_coils_v20140408.mat';   

in_filename_complete = [in_foldername in_filename];

disp(['Loading LTX coils from: ' in_filename_complete]);

LTX_coils = load(in_filename_complete);
