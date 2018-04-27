filebase1 = 'LCFS_MirFL14_4p';
filebase2 = 'outcurve_MirFL14_4p_RLCFS';
title_str = 'Mirror Flip14 4%';


check_descur_output_HSX([filebase1 '_1_503.mat'], [filebase2 '_1_503']);
check_descur_output_HSX([filebase1 '_1_504.mat'], [filebase2 '_1_504']);
check_descur_output_HSX([filebase1 '_1_505.mat'], [filebase2 '_1_505']);
check_descur_output_HSX([filebase1 '_1_506.mat'], [filebase2 '_1_506']);
check_descur_output_HSX([filebase1 '_1_507.mat'], [filebase2 '_1_507']);
check_descur_output_HSX([filebase1 '_1_508.mat'], [filebase2 '_1_508']);
check_descur_output_HSX([filebase1 '_1_509.mat'], [filebase2 '_1_509']);
check_descur_output_HSX([filebase1 '_1_510.mat'], [filebase2 '_1_510']);
check_descur_output_HSX([filebase1 '_1_511.mat'], [filebase2 '_1_511']);
check_descur_output_HSX([filebase1 '_1_512.mat'], [filebase2 '_1_512']);
check_descur_output_HSX([filebase1 '_1_513.mat'], [filebase2 '_1_513']);
check_descur_output_HSX([filebase1 '_1_514.mat'], [filebase2 '_1_514']);
check_descur_output_HSX([filebase1 '_1_515.mat'], [filebase2 '_1_515']);
check_descur_output_HSX([filebase1 '_1_516.mat'], [filebase2 '_1_516']);
check_descur_output_HSX([filebase1 '_1_517.mat'], [filebase2 '_1_517']);
check_descur_output_HSX([filebase1 '_1_518.mat'], [filebase2 '_1_518']);


figure(1);

grid on; axis equal
axis([1.375 1.475 .165 .25])
title(title_str)
make_my_plot_pretty3
set(gcf, 'Position', [200 200 750 550])

