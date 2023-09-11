%%%% Antibody Strip Analysis - Grayscale
%%%% Damien Batchelor
%%%% May 2020
%%%% University of Leeds
clear all
close all
clc
%%
global file scale  %Define file and scale as global variables.
[img_rgb,img_gray,file] = import_img; %Import image, create folder, and convert to grayscale(if RGB)
[red,green,blue,black,white,scale] = calibration(img_rgb);%Calibration for normalisation boxes.
[n,strip,name]=strip_detection(img_rgb); %User input for strip location.
[~,img_gray_norm] = gray_normalisation(img_gray,black,white); %Normalise grayscale image
[pk_data,gray_norm,gray_norm_std,dx_gray] = gray_analysis(img_gray_norm,n,strip);%Grayscale image analysis
grayscale_plot(gray_norm,dx_gray,pk_data,n,name); %Plot and write data for grayscale image.
grayscale_write(gray_norm,gray_norm_std,dx_gray,pk_data,n,name);
%Run these functions if the image is RGB
if size(img_rgb,3) == 3
    %Normalise image, run analysis and write and plot data.
    [norm_rgb] =rgb_normalisation(img_rgb,red,green,blue,black); %Calculate Normalisation Values for RGB
    [rgb_norm,rgb_norm_std,dx_rgb] = rgb_analysis(img_rgb,norm_rgb,n,strip,name);
    rgb_plot(rgb_norm,dx_rgb,n,name);
    rgb_write(rgb_norm,rgb_norm_std,dx_rgb,n,name);
end
%% List of Functions
% Function to import image, determing RGB or Gray, convert if needed.
% Make new folder for data to be saved.
function [img_rgb,img_gray,file] = import_img
[file, folder] = uigetfile('*.jpg'); % GUI for choosing file - saves file name and folder path
cd(folder); img_rgb  = imread(file); %img_rgb = imrotate(img_rgb,180);
if size(img_rgb,3) == 3
    img_gray = sum(img_rgb,3)/3;
elseif size(img_rgb,3) == 1
    img_gray = 255*im2double(img_rgb);
else
    return
end
mkdir([file(1:end-4)]); newFolder = fullfile(folder,file(1:end-4));
copyfile(file,newFolder); cd(newFolder);%Change directory to folder for saving
figure();imshow(img_rgb); %Show image and pop up instructions for drawing ROI calibration
end

%GUI Input for calibration points. Determine position of normalisation
%boxes and wait for user confirmation.
function [red,green,blue,black,white,scale] = calibration(img)
global file
figure();imshow(img); %Show image and pop up instructions for drawing ROI calibration
pos = ["Top Left" "Top Right" "Bottom Left" "Bottom Right"];
x = zeros(1,length(pos)); y=zeros(1,length(pos));
for j = 1:length(pos)
    str = "Click in the centre of the " +pos{j}+" +";
    msgfig = msgbox(str,'Success','modal');
    uiwait(msgfig)
    [x(j),y(j)] = ginput(1);
end
ds = [0.5*((x(2) - x(1))+(x(4)-x(3))) 0.5*((y(3) - y(1))+(y(4)-y(2)))]; %Distance between crosses for scaling (averaged)
scaling  = ds ./ [81 51.5]; scale = mean(scaling); %px/mm
r_h(1)   = 10 ; r_w(1) = 25 ;% (px) % [x y w h]; rectangle coords
red{1}   = [(x(1)+scale*(0.1*r_w(1)-9.2)) (y(1)+scale*(0.1*r_h(1)-9.2 -r_h(1))) (0.8*r_w(1)*scale) (0.7*r_h(1)*scale)];
green{1} = red{1}  + [r_w(1)*scale 0 0 0]; blue{1}  = green{1} + [r_w(1)*scale 0 0 0];
black{1} = blue{1} + [r_w(1)*scale 0 0 0];
white{1} = [red{1}(1) red{1}(2) 0.7*red{1}(3) red{1}(4)] + scale*[-20 0 0 0];
r_h(2)   = 17.5; r_w(2) = 10;
red{2}   = [(x(2)+scale*(9.5+0.2*r_w(2))) y(2)-scale*(9.2-0.1*r_h(2)) scale*0.7*r_w(2) scale*0.8*r_h(2)] ;
green{2} = red{2} + [0 r_h(2)*scale 0 0]; blue{2} = green{2} + [0 r_h(2)*scale 0 0];
black{2} = blue{2} + [0 r_h(2)*scale 0 0];
white{2} = [red{2}(1) white{1}(2) white{1}(3) white{1}(4)];
red{3}   = [(x(3)-scale*(9.2+0.9*r_w(2))) (y(3)-scale*(8.1-0.1*r_h(1))) (0.7*r_w(2)*scale) (0.8*r_h(2)*scale)];
green{3} = red{3}- [0 r_h(2)*scale 0 0]; blue{3}  = green{3} - [0 r_h(2)*scale 0 0];
black{3} = blue{3} - [0 r_h(2)*scale 0 0];
white{3} = [white{1}(1) black{3}(2) white{1}(3) white{1}(4)];
red{4}   = [(x(4)-scale*(15.4-0.1*r_w(1))) (y(4)+scale*(9.2+0.1*r_h(1))) (0.8*r_w(1)*scale) (0.7*r_h(1)*scale)];
green{4} = red{4} - [r_w(1)*scale 0 0 0]; blue{4}  = green{4} - [r_w(1)*scale 0 0 0];
black{4} = blue{4} - [r_w(1)*scale 0 0 0];
white{3} = [white{1}(1) black{4}(2) white{1}(3) white{1}(4)];
white{4} = [black{2}(1) black{4}(2) white{1}(3) white{1}(4)];
figure(); imshow(img);
for j=1:4
    images.roi.Rectangle(gca,'Position',red{j},'FaceAlpha',0,'Label', 'Red');
    images.roi.Rectangle(gca,'Position',green{j},'FaceAlpha',0,'Label','Green');
    images.roi.Rectangle(gca,'Position',blue{j},'FaceAlpha',0,'Label','Blue');
    images.roi.Rectangle(gca,'Position',black{j},'FaceAlpha',0,'Label','Black');
    images.roi.Rectangle(gca,'Position',white{j},'FaceAlpha',0,'Label','White');
end
answer = questdlg('Are the Boxes calibrated correctly?', ...
    'Calibration Check','Yes','No','No');
%         uiwait;
switch answer
    case 'Yes'
        
    case 'No'
        disp('Try again')
        error('Not calibrated correctly');
end
saveas(gcf,strcat(file,'_Calibration.tif'));
end

% Normalised RGB image using calibration data.
function [norm_rgb] = rgb_normalisation(img_rgb,red,green,blue,black)
global file 
[norm_r,norm_g,norm_b,norm_blk_r,norm_blk_g,norm_blk_b] = deal(cell(1,4));
for j = 1:4
    img_r = imcrop(img_rgb,red{j});
    norm_r{j}   = mean2(img_r(:,:,1)) - 0.5*(mean2(img_r(:,:,2))+mean2(img_r(:,:,3)));
    img_g = imcrop(img_rgb,green{j});
    norm_g{j} = mean2(img_g(:,:,2)) - 0.5*(mean2(img_g(:,:,1))+mean2(img_g(:,:,3)));
    img_b = imcrop(img_rgb,blue{j});
    norm_b{j}= mean2(img_b(:,:,3)) - 0.5*(mean2(img_b(:,:,1))+mean2(img_b(:,:,2)));
    img_blk = imcrop(img_rgb,black{j});
    norm_blk_r{j} = mean2(img_blk(:,:,1)) - 0.5*(mean2(img_blk(:,:,2))+mean2(img_blk(:,:,3)));
    norm_blk_g{j} = mean2(img_blk(:,:,2)) - 0.5*(mean2(img_blk(:,:,1))+mean2(img_blk(:,:,3)));
    norm_blk_b{j} = mean2(img_blk(:,:,3)) - 0.5*(mean2(img_blk(:,:,1))+mean2(img_blk(:,:,2)));
end
norm_rgb{1} = [ mean([norm_r{:}])  mean([norm_g{:}])  mean([norm_b{:}]) ];
norm_rgb{2} = [ mean([norm_blk_r{:}])  mean([norm_blk_g{:}])  mean([norm_blk_b{:}]) ];
figure();
X = categorical({'Top' 'Bottom' 'Left' 'Right'});
X = reordercats(X,{'Top' 'Bottom' 'Left' 'Right'});
subplot(1,3,1); bar(X,[norm_r{1} norm_r{3} norm_r{4} norm_r{2}],'r');
title('Normalised Red');lim = ylim; ylim([lim(1) 1.1*lim(2)]);
subplot(1,3,2); bar(X,[norm_g{1} norm_g{3} norm_g{4} norm_g{2}],'g');
title('Normalised Green');lim = ylim; ylim([lim(1) 1.1*lim(2)]);
subplot(1,3,3); bar(X,[norm_b{1} norm_b{3} norm_r{4} norm_b{2}],'b');
title('Normalised Blue');lim = ylim; ylim([lim(1) 1.1*lim(2)]);
savefig(strcat(file,'_CalibrationBoxes_RGB.fig'));
saveas(gcf,strcat(file,'_CalibrationBoxes_RGB.tif'));
end
%Normalise Grayscale image using calibration data.
function [norm_gray,img_gray_norm] = gray_normalisation(img_gray,black,white)
global file 
[norm_gray_blk, norm_gray_wht] = deal(cell(1,4)); 
for j = 1:4
    img_blk = imcrop(img_gray,black{j});
    norm_gray_blk{j} = mean2(img_blk);
    img_wht = imcrop(img_gray,white{j});
    norm_gray_wht{j} = mean2(img_wht);
end

norm_gray = [mean([norm_gray_wht{:}]) mean([norm_gray_blk{:}])];
img_gray_norm = 255*(img_gray - norm_gray(2)) /  ( norm_gray(1) - norm_gray(2));
figure('visible', 'off');; imshowpair(img_gray,img_gray_norm,'montage');
saveas(gcf,strcat(file,'_ImageMontage_Normalised.tif'));
figure();
X = categorical({'Top' 'Bottom' 'Left' 'Right'});
X = reordercats(X,{'Top' 'Bottom' 'Left' 'Right'});
subplot(1,2,1); bar(X,[norm_gray_blk{1} norm_gray_blk{3} norm_gray_blk{4} norm_gray_blk{2}],'k');
title('Mean Black');lim = ylim; ylim([lim(1) 1.1*lim(2)]);
subplot(1,2,2); bar(X,[norm_gray_wht{1} norm_gray_wht{3} norm_gray_wht{4} norm_gray_wht{2}],'w');
title('Mean White');lim = ylim; ylim([lim(1) 1.1*lim(2)]);
savefig(strcat(file,'_CalibrationBox_Gray.fig'));
saveas(gcf,strcat(file,'_CalibrationBox_Gray.tif'));
end
% GUI input for strip locations.
function [n,strip,name]=strip_detection(img)
    figure();imshow(img);
    n    = inputdlg('How many strips are there?'); n=str2double(n{1});%Number of strips = number of loops
    %Show image ready for rectangle input
    %Loop over image for number of stripes
    %Preallocation
    strip = cell(1,n); name = cell(1,n);
    for j=1:n
        msgfig = msgbox(sprintf('Draw a Rectangle on Strip %g', j),'Success','modal');
        uiwait(msgfig)
        strip{j}   =  drawrectangle('Label',sprintf('Strip %g',j),'Color',[1 1 1]);
        %User input for sample name
        name(j) = inputdlg(sprintf('Enter Sample Name for Strip %g', j),'Sample Name',[1 50],{' '},'on');
    end
end

%Analysis for RGB channels in RGB image only.
function [rgb_norm,rgb_norm_std,dx_rgb] = rgb_analysis(img_rgb,norm_rgb,n,strip,name)
global file scale
[rgb,rgb_norm,rgb_std,rgb_norm_std,dx_rgb] = deal(cell(1,n));
for j=1:n
    img_s   = imcrop(img_rgb,strip{j}.Position);
    sum_mean = squeeze(mean(img_s,1));
    s_mean = squeeze(mean(img_s,1));%Mean and StdDev across all horizontal profiles
    rgb{j}(:,1) = s_mean(:,1) - 0.5*(s_mean(:,2)+s_mean(:,3));
    rgb{j}(:,2) = s_mean(:,2) - 0.5*(s_mean(:,1)+s_mean(:,3));
    rgb{j}(:,3) = s_mean(:,3) - 0.5*(s_mean(:,1)+s_mean(:,2));
    rgb_norm{j} = (rgb{j} - norm_rgb{2})./(norm_rgb{1} - norm_rgb{2});
    s_stdev = squeeze(std(im2double(img_s),1));
    rgb_std{j} = s_stdev;
    rgb_norm_std{j} = s_stdev./(norm_rgb{1} - norm_rgb{2});
    dx_rgb{j} = (1/scale)*[1:length(rgb_norm{j})];
end
end
%Grayscale analysis for both RGB and Grayscale images.
function [pk_data,gray_norm,gray_norm_std,dx_gray] = gray_analysis(img_gray_norm,n,strip)
global file scale
[gray_norm,gray_norm_std,dx_gray,...
    pk_data] = deal(cell(1,n));
for j=1:n
    img_s   = imcrop(img_gray_norm,strip{j}.Position);
    sum_mean = squeeze(mean(img_s,1));
    gray_std_temp = squeeze(std(im2double(img_s),1));
    gray_norm{j} = 1-sum_mean/255;
    gray_norm_std{j} = gray_std_temp/255;
    dx_gray{j} = (1/scale)*[1:length(gray_norm{j})];
    pk_data{j} = peakfinder(dx_gray{j},gray_norm{j});
end
end
%Plot and write data for both RGB and Grayscale analysis.
function rgb_plot(rgb_norm,dx_rgb,n,name)
global file
figure()
for j=1:n
    subplot(2,round(n/2),j); %Subplot a 2 x n window
    %Plot each sample with corresponding title and r/g/b colours
    h=plot(dx_rgb{j},rgb_norm{j},'LineWidth',2);set(h,{'Color'},{'r';'g';'b'});
    title(name(j)); xlabel('Horizontal Distance (mm)');
    ylabel('Intensity (a.u.)');%ylim([0 3]);
end
%Save Figure as matlab .fig and .tif
savefig(strcat(file,'_Figure_RGB.fig'));
saveas(gcf,strcat(file,'_Figure_RGB.tif'));
%Save to separate excel spreadsheet for each strip.
%File name is 'imagename_samplename'
end
function grayscale_plot(gray_norm,dx_gray,pk_data,n,name)
global file
figure(); 
for j=1:n
    ax = subplot(2,round(n/2),j); %Subplot a 2 x n window
    plot(dx_gray{j},gray_norm{j},'LineWidth',2);
    hold on; plot(pk_data{j}(:,1),pk_data{j}(:,2),'o','MarkerSize',10); hold off
    xlabel('Horizontal Distance (mm)');
    ylabel('Normalised Intensity (a.u.)'); yL{j} = ylim;
end
ylim([ax(:)],[min([yL{:}]) max([yL{:}])]);
savefig(strcat(file,'_Profile_Grayscale.fig'));
saveas(gcf,strcat(file,'_Profile_Grayscale.tif'));
figure();
for j=1:n
    ax = subplot(2,round(n/2),j); %Subplot a 2 x n window
    bar(pk_data{j}(:,5));
    title(name(j)); xlabel('Peak Number');
    ylabel('Peak Area (a.u.)'); yL{j} = ylim;   
end
saveas(gcf,strcat(file,'_PeakArea_Gray.tif'));
savefig(strcat(file,'__PeakArea_Gray.fig'));
end
function grayscale_write(gray_norm,gray_norm_std,dx_gray,pk_data,n,name)
global file
col_header={'Horizontal Distance (px)','Normalised Intensity (a.u.)','StdDev'};
for k=1:n
    xlfilename = strcat(file,'_',name{k},'_Grayscale.xls');
    T = table(dx_gray{k}.',gray_norm{k}.',gray_norm_std{k}.','VariableNames',...
        col_header);
    writetable(T,xlfilename,'WriteVariableNames',true);

end
col_header={'Peak X (mm)','Peak Y (mm)','Peak Prominence (mm)','Peak Width(mm)', 'Peak Area (a.u.)'};
for k=1:n
    xlfilename = strcat(file,'_',name{k},'_Grayscale_PeakData.xls');
    T = table(pk_data{k}(:,1),pk_data{k}(:,2),pk_data{k}(:,3),...
        pk_data{k}(:,4), pk_data{k}(:,5),'VariableNames',...
        col_header);
    writetable(T,xlfilename,'WriteVariableNames',true);

end
end
function rgb_write(rgb_norm,rgb_norm_std,dx_rgb,n,name)
global file scale
col_header={'Horizontal Distance (mm)','Red','Red StdDev','Green',...
    'Green StdDev','Blue','Blue StdDev'};
for k=1:n
    xlfilename = strcat(file,'_',name{k},'_RGB_Profile.csv');
    
    T = table(dx_rgb{k}.',rgb_norm{k}(:,1),rgb_norm_std{k}(:,1),...
        rgb_norm{k}(:,2),rgb_norm_std{k}(:,2),...
        rgb_norm{k}(:,3),rgb_norm_std{k}(:,3),'VariableNames',col_header);
    writetable(T,xlfilename,'WriteVariableNames',true);
end
end
% Used to determine peaks in data.
function [pk_data] = peakfinder(x,y)
 [~,~,~,proms]=findpeaks(y,x,'Annotate','extents');
   [sortedX, ~] = sort(proms(:),'descend');
    [pk_y,pk_x,pk_width,pk_prom]=findpeaks(y,x,'MinPeakProminence',0.9*sortedX(2),'Annotate','extents');
    pk_area = pk_prom .* pk_width;
    pk_data = sortrows([pk_x.' pk_y.' pk_prom.' pk_width.' pk_area.']); 
end

