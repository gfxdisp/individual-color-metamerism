% The horizontal axis is the Standard Deviation of deltaE 2000 calculated by 
% the CMF of 1000 Simulated Observers. The vertical axis is the Standard Deviation 
% given by the corresponding 74 Observers. 
% All points with an average score of 0 are deleted.
AllData = readtable('ObserverMetamerism\ObserverMetamerism\Data\YC_AllData.csv'); % 4884*5
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

E_means = zeros(6, 11);
E_std = zeros(6, 11);
score_means = zeros(6, 11);
score_std = zeros(6, 11);
point_color = zeros(6, 11, 3);
point_display_pattern = zeros(6, 11);

for i = 1:6
    for j = 1:11
        % Check if all scores are zero
        if all(score_Cell{i, j} == 0)
            E_means(i, j) = NaN;  % Set mean to NaN if all scores are zero
            E_std(i, j) = NaN;
            score_means(i, j) = NaN;
            score_std(i, j) = NaN;
        else
            E_means(i, j) = mean(E_Cell{i, j});
            E_std(i, j) = std(E_Cell{i, j});
            score_means(i, j) = mean(score_Cell{i, j});
            score_std(i, j) = std(score_Cell{i, j}); % / sqrt(size(score_Cell{i, j},2));
        end
        point_color(i, j, :) = color_rgb_s(j, :);
        point_display_pattern(i, j) = i;
    end
end

% Flatten the arrays for plotting
E_means_flat = reshape(E_means, [], 1);
E_std_flat = reshape(E_std, [], 1);
score_means_flat = reshape(score_means, [], 1);
score_std_flat = reshape(score_std, [], 1);
point_color_flat = reshape(point_color, [], 3);
point_display_pattern_flat = reshape(point_display_pattern, [], 1);

valid_indices = ~isnan(E_means_flat) & ~isnan(score_means_flat);

% Plot the scatter plot with error bars
figure;
hold on;
% errorbar(E_means_flat(valid_indices), score_means_flat(valid_indices), score_std_flat(valid_indices), 'bo', 'vertical');
% hold on;
% errorbar(E_means_flat(valid_indices), score_means_flat(valid_indices), E_std_flat(valid_indices), 'bo', 'horizontal');
% hold on;

for indice = 1:size(valid_indices,1)
    if (valid_indices(indice,1) == 0)
        continue
    else
        scatter3(E_means_flat(indice), E_std_flat(indice), score_std_flat(indice), 'MarkerFaceColor', point_color_flat(indice, :), 'MarkerEdgeColor', point_color_flat(indice, :));     
    end
end
grid on;
xlabel('deltaE 2000 (Mean of 1000 values)');
ylabel('deltaE 2000 (Standard Deviation of 1000 values)');
zlabel('Score (Standard Deviation of 72 values)');
title('72 Subject Scores and 1000 CMFs deltaE 2000 Standard deviation points');