function  merge_bnorm_and_regcoil_nescout( bnorm_filename, ...
    regcoil_nescout_in_filename, nescin_out_filename)
% Merge bnorm output with regcoil_nescout output to make a nescin file that
% can be used for

nfp = 4


bnorm_data = read_bnorm(bnorm_filename);
regnes_data = read_regcoil_nescout( regcoil_nescout_in_filename, nfp );

% quantities used later
m_all = [];
n_all = [];
rmnc_all = [];
zmns_all = [];

for ii = 1:numzetal
    for jj = 1:numtheta
        for kk = 1:data.mn_surf
            m = data.xm_surf(kk);
            n = data.xn_surf(kk);
            % stellarator symmetry only for now...
            rmnc = data.rmnc_surf(kk);
            zmns = data.zmns_surf(kk);
            angle = m*theta(jj) + n*zetal(ii);
            % angle = m*theta(jj) + data.nfp*zetal(ii);
            cosangle = cos(angle);
            sinangle = sin(angle);
            coszeta = cos(zetal(ii));
            sinzeta = sin(zetal(ii));
            surf_x_allfp(ii, jj) = surf_x_allfp(ii,jj) + ...
                rmnc * cosangle * coszeta;
            surf_y_allfp(ii, jj) = surf_y_allfp(ii,jj) + ...
                rmnc * cosangle * sinzeta;
            surf_z_allfp(ii, jj) = surf_z_allfp(ii,jj) + ...
                zmns * sinangle;
            if ( (ii == 1) && (jj == 1))
                m_all(kk) = m;
                n_all(kk) = n;
                rmnc_all(kk) = rmnc;
                zmns_all(kk) = zmns;
            end
        end
    end
end
 
    
fout_handle = fopen(nescin_out_filename, 'r');

     
    
    
    
    
    
    