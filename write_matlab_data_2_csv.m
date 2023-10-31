% Write Matlab Data to CSV (for Python read)
AllData = readtable('ObserverMetamerism\ObserverMetamerism\Data\YC_AllData.csv'); % 4884*5
AllData(4621:4752,:) = [];
MLE_data = readtable('mle_scaling\color_similarities_scaled.csv');
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
            E = CIE2000deltaE_XYZ(XYZ1,XYZ2,whitepoint( 'd65' )*100);
            E_set = [E_set,E];
        end
        E_Cell{display_pattern_index, color_index} = E_set;
    end   
end
for i = 1:exp_num
    display_1 = table2array(AllData(i,"Display1"));
    display_2 = table2array(AllData(i,"Display2"));
    color_index = table2array(AllData(i,"ColorIndex"));
    display_pattern_index = display_configurations(display_1,display_2);
    score = table2array(AllData(i,"MatchingScore"));
    score_Cell{display_pattern_index,color_index} = [score_Cell{display_pattern_index,color_index}, score];
end

E_mean = zeros(6, 11);
E_std = zeros(6, 11);
score_mean = zeros(6, 11);
score_std = zeros(6, 11);
MLE_mean = zeros(6, 11);
MLE_std = zeros(6, 11);
point_color = zeros(6, 11, 3);
point_color_index = zeros(6, 11);
point_display_pattern = zeros(6, 11);

for i = 1:6
    for j = 1:11
        % Check if all scores are zero
        if all(score_Cell{i, j} == 0)
            E_mean(i, j) = NaN;  % Set mean to NaN if all scores are zero
            E_std(i, j) = NaN;
            score_mean(i, j) = NaN;
            score_std(i, j) = NaN;
        else
            display_index_1 = display_patterns(i,1);
            display_index_2 = display_patterns(i,2);
            color_index = j;
            query_condition_id = sprintf('%d%d-%d', display_index_1, display_index_2, color_index);
            row_idx = strcmp(MLE_data.condition_id, query_condition_id);
            MLE_mean(i,j) = MLE_data.mle_mean_similarity(row_idx);
            MLE_std(i, j) = MLE_data.mle_std_similarity(row_idx);
            E_mean(i, j) = mean(E_Cell{i, j});
            E_std(i, j) = std(E_Cell{i, j});
            score_mean(i, j) = mean(score_Cell{i, j});
            score_std(i, j) = std(score_Cell{i, j}); % / sqrt(size(score_Cell{i, j},2));
        end
        point_color(i, j, :) = color_rgb_s(j, :);
        point_color_index(i, j) = j;
        point_display_pattern(i, j) = i;
    end
end


E_mean_flat = reshape(E_mean, [], 1);
E_std_flat = reshape(E_std, [], 1);
score_mean_flat = reshape(score_mean, [], 1);
score_std_flat = reshape(score_std, [], 1);
MLE_mean_flat = reshape(MLE_mean, [], 1);
MLE_std_flat = reshape(MLE_std, [], 1);
point_color_flat = reshape(point_color, [], 3);
point_color_index_flat = reshape(point_color_index, [], 1);
point_display_pattern_flat = reshape(point_display_pattern, [], 1);

valid_indices = ~isnan(E_mean_flat) & ~isnan(score_mean_flat);

results_table = table(E_mean_flat(valid_indices), E_std_flat(valid_indices), ...
                      score_mean_flat(valid_indices), score_std_flat(valid_indices), ...
                      MLE_mean_flat(valid_indices), MLE_std_flat(valid_indices), ...
                      point_color_flat(valid_indices, 1), point_color_flat(valid_indices, 2), point_color_flat(valid_indices, 3), ...
                      point_color_index_flat(valid_indices), point_display_pattern_flat(valid_indices), ...
                      'VariableNames', {'E_mean', 'E_std', 'Score_mean', 'Score_std', ...
                      'MLE_mean', 'MLE_std', 'Color_R', 'Color_G', 'Color_B', ...
                      'Color_Index', 'Display_Pattern'});

% 写入CSV文件
writetable(results_table, 'write_matlab_results_data.csv');