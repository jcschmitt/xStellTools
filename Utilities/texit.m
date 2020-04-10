function str_out = texit(str_in)
str_temp = regexprep(str_in, '\', '\\\');
str_out = regexprep(str_temp, '_', '\\_');

