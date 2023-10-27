% The horizontal axis is the Standard Deviation of deltaE 2000 calculated by 
% the CMF of 1000 Simulated Observers. The vertical axis is the Standard Deviation 
% given by the corresponding 74 Observers. 
% All points with an average score of 0 are deleted.
AllData = readtable('ObserverMetamerism\ObserverMetamerism\Data\YC_AllData.csv'); % 4884*5
MLE_data = readtable('mle_scaling\color_similarities_scaled.csv');
AllData(4621:4752,:) = [];
ciexyz31_1 = table2array(readtable('ObserverMetamerism\ObserverMetamerism\Data\AuxData/ciexyz31_1.csv')); % 471*4
Obs_1000_CMF_struct = load("ObserverMetamerism\ObserverMetamerism\Data\AuxData/IndividualObs_2deg.mat");
Obs_1000_CMFs = Obs_1000_CMF_struct.xyz_CMFs;
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
            for lamda = 390:5:780
                X1 = X1 + display_spd_1(color_index,lamda-379) * Obs_1000_CMFs((lamda-385)/5,1,obs) * 5;
                X2 = X2 + display_spd_2(color_index,lamda-379) * Obs_1000_CMFs((lamda-385)/5,1,obs) * 5;
                Y1 = Y1 + display_spd_1(color_index,lamda-379) * Obs_1000_CMFs((lamda-385)/5,2,obs) * 5;
                Y2 = Y2 + display_spd_2(color_index,lamda-379) * Obs_1000_CMFs((lamda-385)/5,2,obs) * 5;
                Z1 = Z1 + display_spd_1(color_index,lamda-379) * Obs_1000_CMFs((lamda-385)/5,3,obs) * 5;
                Z2 = Z2 + display_spd_2(color_index,lamda-379) * Obs_1000_CMFs((lamda-385)/5,3,obs) * 5;
            end
            XYZ1 = [X1,Y1,Z1];
            XYZ2 = [X2,Y2,Z2];
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
MLE_mean = zeros(6, 11);
point_color = zeros(6, 11, 3);
point_display_pattern = zeros(6, 11);

for i = 1:6
    for j = 1:11
        % Check if all scores are zero
        if all(score_Cell{i, j} == 0)
            E_mean(i, j) = NaN;  % Set mean to NaN if all scores are zero
            MLE_mean(i, j) = NaN;
        else
            display_index_1 = display_patterns(i,1);
            display_index_2 = display_patterns(i,2);
            color_index = j;
            query_condition_id = sprintf('%d%d-%d', display_index_1, display_index_2, color_index);
            row_idx = strcmp(MLE_data.condition_id, query_condition_id);
            E_mean(i, j) = mean(E_Cell{i, j});
            MLE_mean(i, j) = MLE_data.mle_mean_similarity(row_idx); % / sqrt(size(score_Cell{i, j},2));
        end
        point_color(i, j, :) = color_rgb_s(j, :);
        point_display_pattern(i, j) = i;
    end
end

% Flatten the arrays for plotting
E_mean_flat = reshape(E_mean, [], 1);
MLE_mean_flat = reshape(MLE_mean, [], 1);
point_color_flat = reshape(point_color, [], 3);
point_display_pattern_flat = reshape(point_display_pattern, [], 1);

valid_indices = ~isnan(E_mean_flat);

% Plot the scatter plot with error bars
figure;
hold on;
E_mean_flat_valid = [];
MLE_mean_flat_valid = [];

for indice = 1:size(valid_indices,1)
    if (valid_indices(indice,1) == 0)
        continue
    else
        E_mean_flat_valid = [E_mean_flat_valid;E_mean_flat(indice)];
        MLE_mean_flat_valid = [MLE_mean_flat_valid;MLE_mean_flat(indice)];
        scatter(E_mean_flat(indice), MLE_mean_flat(indice), 'MarkerFaceColor', point_color_flat(indice, :), 'MarkerEdgeColor', point_color_flat(indice, :));
        text(E_mean_flat(indice), MLE_mean_flat(indice) - 0.01, sprintf("{%d, %d}", display_patterns(point_display_pattern_flat(indice),1), ...
            display_patterns(point_display_pattern_flat(indice),2)), 'HorizontalAlignment', 'center');
        
    end
end

X_axis_begin = 0.3;
Y_axis_begin = 4;
Y_axis_gap = 0.1;

correlation_coefficient_Pearson = corr(E_mean_flat_valid, MLE_mean_flat_valid, 'type','Pearson');
correlation_coefficient_Kendall = corr(E_mean_flat_valid, MLE_mean_flat_valid, 'type','Kendall');
correlation_coefficient_Spearman = corr(E_mean_flat_valid, MLE_mean_flat_valid, 'type','Spearman');
text(X_axis_begin, Y_axis_begin - 7 * Y_axis_gap, sprintf('Pearson Correlation: %.2f', correlation_coefficient_Pearson), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
text(X_axis_begin, Y_axis_begin - 8 * Y_axis_gap, sprintf('Kendall Correlation: %.2f', correlation_coefficient_Kendall), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
text(X_axis_begin, Y_axis_begin - 9 * Y_axis_gap, sprintf('Spearman Correlation: %.2f', correlation_coefficient_Spearman), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

for i = 1:11
    y = - (i - 9) * Y_axis_gap + Y_axis_begin;
    scatter(X_axis_begin + 0.03, y, 60, 'MarkerFaceColor', color_rgb_s(i,:), 'MarkerEdgeColor', color_rgb_s(i,:));
    text(X_axis_begin, y, ['Color', num2str(i)], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end

text(X_axis_begin, Y_axis_begin - 3 * Y_axis_gap, 'Display 1 - LG C2', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
text(X_axis_begin, Y_axis_begin - 4 * Y_axis_gap, 'Display 2 - Sony X310', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
text(X_axis_begin, Y_axis_begin - 5 * Y_axis_gap, 'Display 3 - Samsung Laser Projector', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
text(X_axis_begin, Y_axis_begin - 6 * Y_axis_gap, 'Display 4 - ASUS VG246', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

xlim([0, 0.35]);
% ylim([0, 1.5]);

xlabel('deltaE 2000 (Mean of 1000 values)');
ylabel('MLE Mean');
title('MLE Mean phi and 1000 CMFs deltaE 2000 Mean points');