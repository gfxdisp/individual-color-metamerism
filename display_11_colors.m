% 使用Rafal的代码进行sprad2color操作
% data_file = load("E:/About_Cambridge/All Research Projects/Color/ObserverMetamerism/ObserverMetamerism/Data/Spectra/Projector_Spectra");
% data = data_file.Projector_Spectra;
% data_file = load("E:/About_Cambridge/All Research Projects/Color/ObserverMetamerism/ObserverMetamerism/Data/Spectra/VG246_Spectra");
% data = data_file.VG246_Spectra;
data_file = load("E:/About_Cambridge/All Research Projects/Color/ObserverMetamerism/ObserverMetamerism/Data/Spectra/C2_Spectra");
data = data_file.C2_Spectra;
% data_file = load("E:/About_Cambridge/All Research Projects/Color/ObserverMetamerism/ObserverMetamerism/Data/Spectra/X310_Spectra");
% data = data_file.X310_Spectra;
[numColors, numWavelengths] = size(data);
wavelengths = 380:1:780;

sc = sprad2color('xyz1931');
all_colors = zeros(11, 3);


figure;
hold on;
for color_index = 1:numColors
    
    spectrum = data(color_index, :); %380nm ~ 780nm
    color_xyz = sc.get_color(wavelengths, spectrum); % spectrum to xyz
    color_rgb = cm_xyz2rgb(color_xyz);
    color_rgb_normalized = normalizeRGB(color_rgb);
    color_rgb_gamma = real(gammaCorrection(color_rgb_normalized, 1 / 2.2));

    subplot(numColors, 2, (color_index-1)*2 + 1);
    plot(wavelengths, spectrum, 'LineWidth', 2);
    xlabel('Wavelength (nm)');
    ylabel('Intensity');
    title(sprintf("Color index %d Spectra", color_index));
    grid on;

    subplot(numColors, 2, (color_index-1)*2 + 2); % 在第二个subplot中绘制纯色矩形和RGB值
    axis([0 1 0 1]); % 设置坐标范围为[0,1]
    R_display = color_rgb_gamma(1);
    G_display = color_rgb_gamma(2);
    B_display = color_rgb_gamma(3);
    all_colors(color_index, :) = [R_display*256, G_display*256, B_display*256];
    rectangle('Position', [0, 0, 1, 1], 'FaceColor', [R_display, G_display, B_display]);
    text(0.5, 0.5, sprintf('R:%.2f G:%.2f B:%.2f', R_display*256, G_display*256, B_display*256), 'Color', 'white', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

    title(sprintf("Color index %d Visualization", color_index));
    axis off; % 不显示坐标轴
end
hold off;
writematrix(all_colors, "E:/About_Cambridge/All Research Projects/Color/ObserverMetamerism/ObserverMetamerism/Data/all_11_colors_rgb.csv");
saveas(gcf, 'E:/About_Cambridge/All Research Projects/Color/ObserverMetamerism/ObserverMetamerism/Data/YC_Spectra_visualization/all_colors_visualization.png');

function normalized_rgb = normalizeRGB(rgb)
    % Function to normalize RGB values by the white point
    normalized_rgb = rgb ./ (whitepoint( 'd65' ) * 100);
end

function gamma_corrected_rgb = gammaCorrection(rgb, gamma)
    % Function to apply gamma correction to RGB values
    gamma_corrected_rgb = rgb .^ gamma;
end


