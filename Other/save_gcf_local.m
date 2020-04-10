function save_gcf_local(prefix, force_flag)

if nargin < 2
    force_flag = 0;
end

fh = gcf();
%hgexport(fh, [prefix '.eps']);

if (exist([prefix '.fig']) && ~force_flag)
    disp (['K====' prefix '.fig file exists.']);
    prompt = ['K===Overwrite this file and others (.eps, .pdf, /jpg, .png) Y/n? (default = Y):'];
    str = input(prompt, 's');
    if isempty(str)
        str = 'Y';
    end
else
    str = 'Y';
end

if str == 'Y'
    savefig([prefix '.fig'])
    print([prefix '.eps'], '-deps')
    print([prefix '.pdf'], '-dpdf')
    print([prefix '.jpg'], '-djpeg')
    print([prefix '.png'], '-dpng')
else
    disp('<----Not saving files');
end

