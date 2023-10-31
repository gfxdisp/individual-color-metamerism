% Check and Correct
AllData = readtable('ObserverMetamerism\ObserverMetamerism\Data\YC_AllData.csv'); % 4884*5
MLE_data = readtable('mle_scaling\color_similarities_scaled.csv');
AllData(4621:4752,:) = [];
ciexyz31_1 = table2array(readtable('ObserverMetamerism\ObserverMetamerism\Data\AuxData/ciexyz31_1.csv')); % 471*4
Obs_1000_CMF_struct = load("ObserverMetamerism\ObserverMetamerism\Data\AuxData/IndividualObs_2deg.mat");
Obs_1000_CMFs = Obs_1000_CMF_struct.xyz_CMFs;
original_wavelengths = 390:5:780;
desired_wavelengths = 380:1:780;
Obs_1000_CMFs_1nm = zeros(length(desired_wavelengths), 3, 1000);
for i = 1:1000
    for j = 1:3
        Obs_1000_CMFs_1nm(:, j, i) = interp1(original_wavelengths, Obs_1000_CMFs(:, j, i), desired_wavelengths, 'linear', 'extrap');
    end
end


display_1_spd = table2array(readtable('ObserverMetamerism\ObserverMetamerism\Data\Spectra/C2_Spectra.csv'));
display_2_spd = table2array(readtable('ObserverMetamerism\ObserverMetamerism\Data\Spectra/X310_Spectra.csv'));
display_3_spd = table2array(readtable('ObserverMetamerism\ObserverMetamerism\Data\Spectra/Projector_Spectra.csv'));
display_4_spd = table2array(readtable('ObserverMetamerism\ObserverMetamerism\Data\Spectra/VG246_Spectra.csv'));
display_spds = {display_1_spd; display_2_spd; display_3_spd; display_4_spd};

exp_num = size(AllData,1);

color_rgb_s = load("ObserverMetamerism\ObserverMetamerism\Data/all_11_colors_rgb.csv");
color_rgb_s = color_rgb_s ./ 256.;

display_patterns = [[1 2];[1 3];[1 4];[2 3]; [2 4]; [3 4]];
display_configurations = [[-1 1 2 3];[-1 -1 4 5]; [-1 -1 -1 6]; [-1 -1 -1 -1]];
num_obs = size(Obs_1000_CMFs,3); %1000
E_Cell = cell(6, 11);
score_Cell = cell(6, 11);

for display_pattern_index = 1:6
    display_1 = display_patterns(display_pattern_index, 1);
    display_2 = display_patterns(display_pattern_index, 2);
    for color_index = 1:11
        display_spd_1 = display_spds{display_1};
        display_spd_2 = display_spds{display_2};
        E_set = [];
        for obs = 1:num_obs
            [X1,X2,Y1,Y2,Z1,Z2] = deal(0,0,0,0,0,0);
            for lamda = 380:1:780
                X1 = X1 + display_spd_1(color_index,lamda-379) * Obs_1000_CMFs_1nm(lamda-379,1,obs);
                X2 = X2 + display_spd_2(color_index,lamda-379) * Obs_1000_CMFs_1nm(lamda-379,1,obs);
                Y1 = Y1 + display_spd_1(color_index,lamda-379) * Obs_1000_CMFs_1nm(lamda-379,2,obs);
                Y2 = Y2 + display_spd_2(color_index,lamda-379) * Obs_1000_CMFs_1nm(lamda-379,2,obs);
                Z1 = Z1 + display_spd_1(color_index,lamda-379) * Obs_1000_CMFs_1nm(lamda-379,3,obs);
                Z2 = Z2 + display_spd_2(color_index,lamda-379) * Obs_1000_CMFs_1nm(lamda-379,3,obs);
            end
            XYZ1 = [X1*683,Y1*683,Z1*683];
            XYZ2 = [X2*683,Y2*683,Z2*683];

            color_rgb_1 = cm_xyz2rgb(XYZ1);
            color_rgb_normalized_1 = normalizeRGB(color_rgb_1);
            color_rgb_gamma_1 = real(gammaCorrection(color_rgb_normalized_1, 1 / 2.2));
            color_rgb_2 = cm_xyz2rgb(XYZ2);
            color_rgb_normalized_2 = normalizeRGB(color_rgb_2);
            color_rgb_gamma_2 = real(gammaCorrection(color_rgb_normalized_2, 1 / 2.2));

            E = CIE2000deltaE_XYZ(XYZ1,XYZ2,whitepoint( 'd65' )*100);
            E_set = [E_set,E];
        end
        E_Cell{display_pattern_index, color_index} = E_set;
    end   
end

function normalized_rgb = normalizeRGB(rgb)
    % Function to normalize RGB values by the white point
    normalized_rgb = rgb ./ (whitepoint( 'd65' ) * 100);
end

function gamma_corrected_rgb = gammaCorrection(rgb, gamma)
    % Function to apply gamma correction to RGB values
    gamma_corrected_rgb = rgb .^ gamma;
end

