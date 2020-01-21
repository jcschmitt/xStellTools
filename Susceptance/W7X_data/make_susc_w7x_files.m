function sdfdf()
% Inputs
%   1: The path to the VMEC output 'wout' file
%   2: The path to the space-delimited text file from 'xtest_stuff', a.k.a.
%   test_stuff.F90, which contains the susceptance matrix
%
% Outputs
%  The following data is stored in a .mat file.
%    Bsquared     array              
%    S11          array             
%    S12        array             
%    S21          array          
%    S22       array            
%    Vprime        array           
%    r_eff      array

make_figure = 1;

ii = 1;
wout_file_list{ii} = 'Kisslinger_ideal_v1/STELLOPT/wout_vmec_8x8_ns_51.nc';
text_file_list{ii} = 'Kisslinger_ideal_v1/STELLOPT/SUS_51_mod.txt';
output_filenames{ii} = 'Kisslinger_Ideal_v1_ns_51_stellopt.mat';

ii = ii + 1;
wout_file_list{ii} = 'Kisslinger_ideal_v1/STELLOPT/wout_vmec_8x8_ns_101.nc';
text_file_list{ii} = 'Kisslinger_ideal_v1/STELLOPT/SUS_101_mod.txt';
output_filenames{ii} = 'Kisslinger_Ideal_v1_ns_101_stellopt.mat';

ii = ii + 1;
wout_file_list{ii} = 'Kisslinger_ideal_v1/STELLOPT/wout_vmec_8x8_ns_201.nc';
text_file_list{ii} = 'Kisslinger_ideal_v1/STELLOPT/SUS_201_mod.txt';
output_filenames{ii} = 'Kisslinger_Ideal_v1_ns_201_stellopt.mat';

ii = ii + 1;
wout_file_list{ii} = 'Kisslinger_ideal_v2/STELLOPT/wout_vmec_8x8_ns_51.nc';
text_file_list{ii} = 'Kisslinger_ideal_v2/STELLOPT/SUS_51_mod.txt';
output_filenames{ii} = 'Kisslinger_Ideal_v2_ns_51_stellopt.mat';

ii = ii + 1;
wout_file_list{ii} = 'Kisslinger_ideal_v2/STELLOPT/wout_vmec_8x8_ns_101.nc';
text_file_list{ii} = 'Kisslinger_ideal_v2/STELLOPT/SUS_101_mod.txt';
output_filenames{ii} = 'Kisslinger_Ideal_v2_ns_101_stellopt.mat';

ii = ii + 1;
wout_file_list{ii} = 'Kisslinger_ideal_v2/STELLOPT/wout_vmec_8x8_ns_201.nc';
text_file_list{ii} = 'Kisslinger_ideal_v2/STELLOPT/SUS_201_mod.txt';
output_filenames{ii} = 'Kisslinger_Ideal_v2_ns_201_stellopt.mat';


number_of_configs = length(wout_file_list);

sfig = figure;

for ii = 1:number_of_configs
    try
        disp(['<----Beginning iteration ', num2str(ii), ' of ' ...
            num2str(number_of_configs)]);
        disp(['<----Opening vmec file: ' wout_file_list{ii}]);
        VMEC_data = read_vmec(wout_file_list{ii});
        disp(['<----Opening S file: ' text_file_list{ii}]);
        S_data = dlmread(text_file_list{ii});
        
        disp('<----Initializing data');
        B = 0;
        Bsquared = 0;
        S11 = 0;
        S12 = 0;
        S21 = 0;
        S22 = 0;
        Vprime = 0;
        r_eff = 0;
        
        
        disp('<----Filling data');
        Bsquared = VMEC_data.bdotb;
        disp('<----Not filling in B');
        B = NaN * Bsquared;

        S11 = S_data(:, 8);
        S12 = S_data(:, 9);
        S21 = S_data(:, 10);
        S22 = S_data(:, 11);
        Vprime = S_data(:, 3);
        
        % This is a check
        Vprime2 = VMEC_data.vp;
        
        % Minor radius * rho
        r_eff = VMEC_data.Aminor * S_data(:, 2);
        
        disp(['<----Saving data: ' output_filenames{ii}]);
        
        save(output_filenames{ii}, 'B', 'Bsquared', 'S11', 'S12', 'S21', ...
            'S22', 'Vprime', 'r_eff');
    

        if make_figure
            figure(sfig);
            subplot(2,4,1);hold on; box on
            plot(r_eff, S11, '.-')
            title('S11')
            subplot(2,4,2);hold on; box on
            plot(r_eff, S12, '.-')
            title('S12')
            subplot(2,4,5);hold on; box on
            plot(r_eff, S21, '.-')
            title('S21')
            subplot(2,4,6);hold on; box on
            plot(r_eff, S22, '.-')
            title('S22')
            
            subplot(1,2,2);hold on; box on
            plot(r_eff, -S12 ./ S11, '.-')  
            title('iota');
        end
        
    catch
        disp(['<----Something failed on iteration ii=' num2str(ii)]);
    end
    
    figure(sfig);
    subplot(1,2,2);
    hold on;
    legend(wout_file_list);
    ylim([.8 1.2])
    
end






