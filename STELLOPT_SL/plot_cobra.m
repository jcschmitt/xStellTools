function f = plot_cobra(cobra_extensions, title_prefix)
%READ_COBRA(filename) Reads a COBRA output file into a COBRA structure.
%
%   READ_COBRA(filename) Reads a COBRA output file into a COBRA structure.
%
%   Example:
%       cobra_data=read_cobra('grate_cobra.test');
%
%   Written by: Samuel Lazerson (lazerson@pppl.gov)
%   Version:    1.5
%   Date:       09/28/2011

for ii = 1:length(cobra_extensions)
    cobradata{ii} = read_cobra(['cobra_grate.aten_scbs_bootsj_beta_1p7.' cobra_extensions{ii}]);
end

for ii = 1:length(cobra_extensions)
figure;box on; hold on; plot([1:49]/50, cobradata{ii}.grate', 'o');
title([title_prefix, ' ', cobra_extensions{ii}]);
end