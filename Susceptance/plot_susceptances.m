function plot_susceptances()

S_QHS = load('susceptance_matrix_QHS_Rstart_1_513_AXIS_v2_VMEC_8p46.mat');
S_MirFL14_10P = load('susceptance_matrix_MirFL14_10p_1_508_20110315.mat');
S_Well_11P = load('susceptance_matrix_Well11p_20110315.mat');
S_Hill_11P = load('susceptance_matrix_Hill11p_20110315.mat');
S_QHS_equiv_tok = load('susceptance_matrix_qhs_equivtok_vacuum.mat');
S_qps14 = load('susceptance_matrix_qps14.mat');
S_qo = load('susceptance_matrix_qo.mat');
S_qohb = load('susceptance_matrix_qohb.mat');
S_ncsx_fixed = load('susceptance_matrix_ncsx_fixed.mat');
S_ncsx_qas3 = load('susceptance_matrix_ncsx_qas3.mat');
S_w7x_kissliner_ideal_v1_201 = load('Kisslinger_Ideal_v1_ns_201_xstelltools.mat');
legend_text = {'QHS', 'MirFL1r 10%', 'Well 11%', 'Hill 11%', 'QHS Equiv. Tok.', ...
    'QPS14', 'W7-X', 'W7-X High Beta', 'NCSX Fixed Bound', 'NCSX qas13', 'KI_v1_201'};

S_QHS.linestyle = 'b-';
S_MirFL14_10P.linestyle = 'k--';
S_Well_11P.linestyle = 'r:';
S_Well_11P.linestyle = 'r.';
S_Hill_11P.linestyle = 'r--';
S_QHS_equiv_tok.linestyle = 'k--';
S_qps14.linestyle = 'gx';
S_qo.linestyle = 'cs';
S_qohb.linestyle = 'c+';
S_ncsx_fixed.linestyle = 'mp';
S_ncsx_qas3.linestyle = 'ms';
S_w7x_kissliner_ideal_v1_201.linestyle = ':';

figure
subplot(2,2,1);box on; hold on;
plot(100*S_QHS.r_eff, S_QHS.S11, S_QHS.linestyle);
plot(100*S_MirFL14_10P.r_eff, S_MirFL14_10P.S11, S_MirFL14_10P.linestyle);
plot(100*S_Well_11P.r_eff, S_Well_11P.S11, S_Well_11P.linestyle);
plot(100*S_Hill_11P.r_eff, S_Hill_11P.S11, S_Hill_11P.linestyle);
plot(100*S_QHS_equiv_tok.r_eff, S_QHS_equiv_tok.S11, S_QHS_equiv_tok.linestyle);
plot(100*S_qps14.r_eff, S_qps14.S11, S_qps14.linestyle);
plot(100*S_qo.r_eff, S_qo.S11, S_qo.linestyle);
plot(100*S_qohb.r_eff, S_qohb.S11, S_qohb.linestyle);
plot(100*S_ncsx_fixed.r_eff, S_ncsx_fixed.S11, S_ncsx_fixed.linestyle);
plot(100*S_ncsx_qas3.r_eff, S_ncsx_qas3.S11, S_ncsx_qas3.linestyle);
plot(100*S_w7x_kissliner_ideal_v1_201.r_eff, S_w7x_kissliner_ideal_v1_201.S11, S_w7x_kissliner_ideal_v1_201.linestyle);
xlim([0 12.5])
ylim([-.2 0.02])

subplot(2,2,2);box on; hold on;
plot(100*S_QHS.r_eff, S_QHS.S12, S_QHS.linestyle);
plot(100*S_MirFL14_10P.r_eff, S_MirFL14_10P.S12, S_MirFL14_10P.linestyle);
plot(100*S_Well_11P.r_eff, S_Well_11P.S12, S_Well_11P.linestyle);
plot(100*S_Hill_11P.r_eff, S_Hill_11P.S12, S_Hill_11P.linestyle);
plot(100*S_QHS_equiv_tok.r_eff, S_QHS_equiv_tok.S12, S_QHS_equiv_tok.linestyle);
plot(100*S_qps14.r_eff, S_qps14.S12, S_qps14.linestyle);
plot(100*S_qo.r_eff, S_qo.S12, S_qo.linestyle);
plot(100*S_qohb.r_eff, S_qohb.S12, S_qohb.linestyle);
plot(100*S_ncsx_fixed.r_eff, S_ncsx_fixed.S12, S_ncsx_fixed.linestyle);
plot(100*S_ncsx_qas3.r_eff, S_ncsx_qas3.S12, S_ncsx_qas3.linestyle);
plot(100*S_w7x_kissliner_ideal_v1_201.r_eff, S_w7x_kissliner_ideal_v1_201.S12, S_w7x_kissliner_ideal_v1_201.linestyle);
xlim([0 12.5])
ylim([-.02 0.2])

subplot(2,2,3);box on; hold on;
plot(100*S_QHS.r_eff, S_QHS.S21, S_QHS.linestyle);
plot(100*S_MirFL14_10P.r_eff, S_MirFL14_10P.S21, S_MirFL14_10P.linestyle);
plot(100*S_Well_11P.r_eff, S_Well_11P.S21, S_Well_11P.linestyle);
plot(100*S_Hill_11P.r_eff, S_Hill_11P.S21, S_Hill_11P.linestyle);
plot(100*S_QHS_equiv_tok.r_eff, S_QHS_equiv_tok.S21, S_QHS_equiv_tok.linestyle);
plot(100*S_qps14.r_eff, S_qps14.S21, S_qps14.linestyle);
plot(100*S_qo.r_eff, S_qo.S21, S_qo.linestyle);
plot(100*S_qohb.r_eff, S_qohb.S21, S_qohb.linestyle);
plot(100*S_ncsx_fixed.r_eff, S_ncsx_fixed.S21, S_ncsx_fixed.linestyle);
plot(100*S_ncsx_qas3.r_eff, S_ncsx_qas3.S21, S_ncsx_qas3.linestyle);
plot(100*S_w7x_kissliner_ideal_v1_201.r_eff, S_w7x_kissliner_ideal_v1_201.S21, S_w7x_kissliner_ideal_v1_201.linestyle);
xlim([0 12.5])
ylim([-.02 0.2])

legend(legend_text);

subplot(2,2,4);box on; hold on;
plot(100*S_QHS.r_eff, S_QHS.S22, S_QHS.linestyle);
plot(100*S_MirFL14_10P.r_eff, S_MirFL14_10P.S22, S_MirFL14_10P.linestyle);
plot(100*S_Well_11P.r_eff, S_Well_11P.S22, S_Well_11P.linestyle);
plot(100*S_Hill_11P.r_eff, S_Hill_11P.S22, S_Hill_11P.linestyle);
plot(100*S_QHS_equiv_tok.r_eff, S_QHS_equiv_tok.S22, S_QHS_equiv_tok.linestyle);
plot(100*S_qps14.r_eff, S_qps14.S22, S_qps14.linestyle);
plot(100*S_qo.r_eff, S_qo.S22, S_qo.linestyle);
plot(100*S_qohb.r_eff, S_qohb.S22, S_qohb.linestyle);
plot(100*S_ncsx_fixed.r_eff, S_ncsx_fixed.S22, S_ncsx_fixed.linestyle);
plot(100*S_ncsx_qas3.r_eff, S_ncsx_qas3.S22, S_ncsx_qas3.linestyle);
plot(100*S_w7x_kissliner_ideal_v1_201.r_eff, S_w7x_kissliner_ideal_v1_201.S22, S_w7x_kissliner_ideal_v1_201.linestyle);
xlim([0 12.5])
ylim([-150 0])

make_my_plot_pretty3
