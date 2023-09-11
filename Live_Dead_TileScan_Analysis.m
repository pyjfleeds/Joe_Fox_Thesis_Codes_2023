%% Joe.F Tilescan Analysis and Quantification Nov 2022.

%% Resets the workspace
clear all
close all
clc

%% Reads in Live Stained Image

folder = uigetdir;
cd(folder);
I_dead = imread(ls('**/*ch02.tif')); % Reads in live stained image
I_live = imread(ls('**/*ch00.tif')); % Reads in dead stained image

%Creates Destination Folder For Saving Results
dest_dir = fullfile(folder,'Results_20230116_Thresholds_5');   %Change name of folder here
mkdir(dest_dir);

%% Thresholding Dead, FULL MERGED
%Generates blank matrix for manipulation
I_dead_thresh_merged = zeros(size(I_dead,1), size(I_dead,2)); 
% threshold value, in example was 128, 8 looks good for zoomed (dead)
t = 5;
% find values below
ind_below = (I_dead < t);
% find values above
ind_above = (I_dead >= t);

%Now edits the blank matrix, keeping the original for comparison
% set values below to black
I_dead_thresh_merged(ind_below) = 0;
% set values above to white
I_dead_thresh_merged(ind_above) = 255;

%Plot comparison
figure();
subplot(2,1,1);
imshow(I_dead);
title('Original - Dead');
subplot(2,1,2);
imshow(I_dead_thresh_merged);
title('Thresholded - Dead');
linkaxes([subplot(2,1,1), subplot(2,1,2)], 'xy');  %Used so that when one image is zoomed, both alter, for ease of comparison


exportgraphics(gcf,fullfile(dest_dir,'Thresholded_Dead_Merged.tif'),'Resolution',300);

%% Thresholding Live, FULL MERGED
%Generates blank matrix for manipulation
I_live_thresh_merged = zeros(size(I_live,1), size(I_live,2)); 
% threshold value, in example was 128, 8 looks good for zoomed (dead)
t = 5;
% find values below
ind_below = (I_live < t);
% find values above
ind_above = (I_live >= t);

%Now edits the blank matrix, keeping the original for comparison
% set values below to black
I_live_thresh_merged(ind_below) = 0;
% set values above to white
I_live_thresh_merged(ind_above) = 255;

%Plot comparison
figure();
subplot(2,1,1);
imshow(I_live);
title('Original - Live');
subplot(2,1,2);
imshow(I_live_thresh_merged);
title('Thresholded - Live');
linkaxes([subplot(2,1,1), subplot(2,1,2)], 'xy');  %Used so that when one image is zoomed, both alter, for ease of comparison

exportgraphics(gcf,fullfile(dest_dir,'Thresholded_Live_Merged.tif'),'Resolution',300);

%% Histogram on Thresholded Merged Images, Bins Vertically
clear Binned_Values_Live %Resets the Values array so bin can be used to maniupulate plot
clear Col_sums_Live
clear Binned_Values_Dead
clear Col_sums_Dead
clear x__bar_vals
clear Bin

%LIVE
Col_sums_Live = sum(I_live_thresh_merged(:,1:1:end),1); %sums all the columns vertically on the thresholded Image
% Combine the column sums in X(Bin) sized blocks
Bin = 50; %This defines how often the summation is made

counter = 1; %Used so values in loop do not overwrite
for i=1:Bin:Bin*floor((size(Col_sums_Live,2))/Bin); %Floor used to round down, get error because array size not divisible by 5
    Binned_Values_Live(counter,:) = sum(Col_sums_Live(:,i:i+(Bin-1)));
    counter = counter+1;
end

%Dead - Repeats the above loopCol_sums_Live
Col_sums_Dead = sum(I_dead_thresh_merged(:,1:1:end),1); %sums all the columns vertically on the thresholded Image
counter = 1; %Used so values in loop do not overwrite
for i=1:Bin:Bin*floor((size(Col_sums_Dead,2))/Bin); %Floor used to round down, get error because array size not divisible by 5
    Binned_Values_Dead(counter,:) = sum(Col_sums_Dead(:,i:i+(Bin-1)));
    counter = counter+1;
end

x_bar_vals = [Bin/2:Bin:size(Col_sums_Live,2)];

%Need this if statement as sometimes rounding gives size mismatch when try
%plot hist. Removes last value from x data so that X and Y can be plotted
if numel(x_bar_vals) ~= numel(Binned_Values_Live)
x_bar_vals(end)=[]; 
end

%Comparative subplot
figure(); 
subplot(2,1,1);
bar(x_bar_vals, Binned_Values_Live)
title('Thresholded+Binned, Live');
xlabel('X Distance (Pixels)');ylabel('Pixel Count');
subplot(2,1,2);
bar(x_bar_vals, Binned_Values_Dead)
title('Thresholded+Binned, Dead');
xlabel('X Distance (Pixels)');ylabel('Pixel Count');
hold off

exportgraphics(gcf,fullfile(dest_dir,'Live_Dead_Pixel_Count_Compared.tif'),'Resolution',300);

%Overlaid Plot
figure();
bar(x_bar_vals, Binned_Values_Live,'FaceColor', 'g')
xlabel('X Distance (Pixels)');ylabel('Pixel Count');
hold on
bar(x_bar_vals, Binned_Values_Dead,'FaceColor', 'r')
xlabel('X Distance (Pixels)');ylabel('Pixel Count');
legend('Live', 'Dead');
hold off

exportgraphics(gcf,fullfile(dest_dir,'Live_Dead_Pixel_Count_Overlay.tif'),'Resolution',300);

Binned_Pct_dead = (Binned_Values_Dead(:,1))./(Binned_Values_Dead(:,1) + Binned_Values_Live(:,1))*100; %.* gives elementwise maths
%% Quantification - Using Binned Values - Applied
Live_Intensity = sum(Col_sums_Live);
Dead_Intensity = sum(Col_sums_Dead);

Pct_dead = (sum(Col_sums_Dead))/(sum(Col_sums_Dead)+sum(Col_sums_Live))*100;
Pct_live = (sum(Col_sums_Live))/(sum(Col_sums_Dead)+sum(Col_sums_Live))*100;

% Comparative Quantification - Need to analyse central area of constant
% width, for image data to be comparable, for instance, one image may have
% analysed a longer region of the channel, and hence will have more live
% cells. Here we crop to only quantify a central region. 

%Cropped Data
MidPoint = ceil((numel(Col_sums_Live))/2);
Lower_Bound = MidPoint - 2250;
Upper_Bound = MidPoint + 2250;

Cropped_Col_sums_Live = Col_sums_Live(Lower_Bound:Upper_Bound);
Cropped_Col_sums_Dead = Col_sums_Dead(Lower_Bound:Upper_Bound);
%Values below offer inter-image comparison. Allows me to compare to DMEM
%control channels
Live_Intensity_Cropped = sum(Cropped_Col_sums_Live); 
Dead_Intensity_Cropped = sum(Cropped_Col_sums_Dead);

%% Fit Gaussian to Percentage Dead Plot
n   =  length(x_bar_vals);
ix = [round(0.1*n):round(0.95*n)];

[f,gof]         = fit(x_bar_vals(ix).',Binned_Pct_dead(ix),'gauss1','Display','iter');
fit_y_vals(1,:) = f(x_bar_vals);
d_fit           =  fit_y_vals(ix) - Binned_Pct_dead(ix);
peak_rmse       = sqrt(sum((d_fit(:).^2)/(n)));
% rmse_err = gof.rmse;
ci              = confint(f,0.95);
sigma           = f.c1;
FWHM            = 2*sqrt(2*log(2))*f.c1;
FWHM_err        = 2*sqrt(2*log(2))*0.5*(ci(6)-ci(5));
pk_height       = f.a1;
pk_x            = f.b1;

r_sq = gof.rsquare;
n = length(x_bar_vals);

figure();
bar(x_bar_vals, Binned_Pct_dead, 'FaceColor', 'r'); 
hold on;
plot(x_bar_vals,fit_y_vals, 'LineWidth',2, 'color', 'k');
xlabel('X Distance (Pixels)');ylabel('Dead Intensity (%)');

exportgraphics(gcf,fullfile(dest_dir,'Gauss_Fit_to_Pct_Dead.tif'),'Resolution',300);
%% Export to CSV

%Pixel Data
col_header= {'X Distance (Pixels)', 'Binned Values Live', 'Binned Values Dead', 'Binned Pct Dead','MatLab Gaussian Fit'};
data = table(x_bar_vals.',Binned_Values_Live,Binned_Values_Dead,Binned_Pct_dead,fit_y_vals.','VariableNames',col_header);
writetable(data,fullfile(dest_dir,'Pixel_Quantification_Data.csv'),'WriteVariableNames',true);

%Gaussian Data
col_header = {'X Peak (mm)','Y Peak (a.u.)','FWHM','FWHM Error','Peak RMSE','Bin Size','Percentage Live Intensity', 'Percentage Dead Intensity', 'Live Intensity', 'Dead Intensity','Cropped Live Intensity','Cropped Dead Intensity','Cropping Midpoint', 'Cropping Upper Bound','Cropping Lower Bound'};
summary_data = table(pk_x',pk_height,FWHM,FWHM_err,peak_rmse,Bin,Pct_live,Pct_dead, Live_Intensity,Dead_Intensity,Live_Intensity_Cropped,Dead_Intensity_Cropped,MidPoint,Upper_Bound,Lower_Bound,'VariableNames',col_header);
writetable(summary_data,fullfile(dest_dir,'Gaussian_Data.csv'),'WriteVariableNames',true);

%Saves workspace
save(fullfile(dest_dir,'Workspace'));

