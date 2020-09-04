function data = read_regcoil_nescout( filename, nfp )
%READ_REGCOIL_NESCOUT
% For parsing the diagnostics output from stellopt/regcoil/nescout runs
%
% Note: The column order should be the same as in NESCIN files

fid = fopen(filename,'r');
if 1 % strfind(filename,'regcoil_nescout')
    fgetl(fid);
    fgetl(fid);
    line=fgetl(fid); temp = sscanf(line,'%d',1);
    data.nfp = temp(1);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    line=fgetl(fid); temp = sscanf(line,'%d',1);
    data.mn_surf = temp(1);
    fgetl(fid);
    fgetl(fid);
    temp = fscanf(fid,'%d %d %e %e %e %e',[6 data.mn_surf]);
    data.xm_surf = temp(1,:);
    data.xn_surf = temp(2,:).*data.nfp;
    data.rmnc_surf = temp(3,:)';
    data.zmns_surf = temp(4,:)';
    data.rmns_surf = temp(5,:)';
    data.zmnc_surf = temp(6,:)';
    if max(abs(data.zmns_surf)) == 0
        disp('<----Looks like an old version of regcoil_nescout data')
        data.rmnc_surf = temp(3,:)';
        data.rmns_surf = temp(4,:)';
        data.zmnc_surf = temp(5,:)';
        data.zmns_surf = temp(6,:)';
    end
end
fclose(fid);
end

