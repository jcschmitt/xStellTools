filebase = 'LCFS_MirFL14_5p';
title_str = 'Mirror Flip14 5%';

make_descur_file_HSX([filebase '_1_503.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_504.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_505.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_506.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_507.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_508.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_509.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_510.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_511.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_512.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_513.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_514.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_515.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_516.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_517.mat'], '', 50, 800);
make_descur_file_HSX([filebase '_1_518.mat'], '', 50, 800);
figure(401);
legend('R_{LCFS} = 1.503', '1.504', '1.505', '1.506', '1.507', '1.508', '1.509', '1.510', '1.511', '1.512', '1.513', '1.514', '1.515', '1.516', '1.517', '1.518');grid on; axis equal
axis([1.375 1.475 .165 .25])
title(title_str)
make_my_plot_pretty3
set(gcf, 'Position', [200 200 750 550])

