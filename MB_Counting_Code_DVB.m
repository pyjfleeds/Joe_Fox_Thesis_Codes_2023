clear all
close all
clc
%%
global bw_thresh circle_sens
% This function opens a GUI to select the folder in which the .jpg image
% files are located.
[myFolder,filePattern,theFiles] = findfolder();
% User input for required values.
[scaling,dilution,BinWidth,spacer_height]=calibration();
%GUI Thresholding for binarisation and circle detection sensitivity.
gui_threshold(theFiles,myFolder,scaling);
% Preallocation for variables.
[centers,radii_px,radii_um] = deal(cell(1,length(theFiles)));
concentration = zeros(1,length(theFiles));
% Loop over all the .jpg image files in the folder.
for k = 1:length(theFiles)
    fprintf(1, "Now reading %g of %g \n",[k,length(theFiles)]); 
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    img          = imread(fullFileName);
    %Convert image to grayscale and then binarize.
    gray_img     = grayscale_process(img);
    bw           = bw_process(gray_img,bw_thresh);
    %Analyse image to detection circles and remove false detections.
    [centers{k},radii_px{k}] = analysis_thresh(bw,circle_sens,scaling);
    radii_um{k} = scaling*radii_px{k};
    concentration(k)       = sum(~isnan(radii_um{k}));
    % Preview windows showing where the bubbles were detected in the first
    % 3 images.
    if k < 4
        figure(); %Create figure
        imshow(img,[]);%show together rgb image and fill image
        h = viscircles(centers{k},radii_px{k});  %display the detected circles
    end
end
%% Data Processing
[height, width, ~] = size(img); %find dimensions of image
V_img = height*width*scaling^2 * spacer_height/1e+12;
img_conc = concentration*dilution/V_img;
conc_avg = mean(img_conc);
conc_err = std([img_conc(:)]);
diam_all = 2*vertcat(radii_um{:});
diam_all = diam_all(~isnan(diam_all));
max_diam = max(diam_all);
diam_avg = mean(diam_all,'omitnan');
diam_std = std(diam_all,'omitnan');
N = length(diam_all);
number1  = length(find(diam_all < 8))*100/N; %Number of bubbles less than 8 or 10 microns
number2  = length(find(diam_all < 10))*100/N ;
nbins    = round(max_diam,0);
RSD      = (diam_std*100)/diam_avg;
DiamFile = fullfile(myFolder, 'DiameterList.txt'); %Save list of diameters as text file
save(DiamFile,'diam_all','-ascii');
%% Plot Bar Chart
%Define edges of bar chart
edges             = 0:BinWidth:ceil(max_diam);
Bins              = histc(diam_all,edges); %Get bin values based off the edges defined
BinConcentration  = (Bins * dilution)/(V_img*length(radii_um)); %Convert Counts into Concentration
EdgesNew          = edges-BinWidth/2; %Define the new edges to be centred
pos = get(0,'ScreenSize');
figure('position',[100 100 pos(3)/1.5 pos(4)/1.5]); 
b       = bar(EdgesNew,BinConcentration,'BarWidth',1,'EdgeColor','black'); %Plot the concentration bar chart
set(gca,'FontSize',12); 
xlabel('Diameter (\mum)'); %Give graph xlabel
ylabel('Concentration (/mL)');%y label
set(gca,'FontSize',14); %Set font size
xlim([0 10]); % X limits
set(gca,'FontSize',14) 
dim   = [.50 .2  .8 .6];%Defines position and height of the text box
str_g = sprintf('%0.3g +/- %0.3g',[conc_avg,conc_err]);
hold on
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
str = {(['Concentration = ' str_g ' bubbles/mL']),(['Mean = ' num2str(diam_avg,3) ' +/- ' num2str(diam_std,2) '\mum']), ...
    (['RSD = ' num2str(RSD,3) '%']), (['% of Bubbles < 10\mum =  ' num2str(number2,4)]), ...
   (['% of Bubbles < 8\mum =  ' num2str(number1,4)]), (['Max Diameter = ' num2str(max_diam,3) '\mum' ])...
    ...
   
 };
t       = annotation('textbox',dim,'String',str,'FitBoxToText','on');
t.FontSize = 10;
%Create legend stating number of bubbles
legend(['N =  ' num2str(N) ])
BinCentre         = edges+BinWidth/2;
%Define bin centres and save these and concentrations to csv file
col_header={'Bin Centre (um)','Concentration (/mL)'};
T = table(BinCentre.',BinConcentration,'VariableNames',col_header);
writetable(T,'Concentration.csv','WriteVariableNames',true);
 % Save histogram as .tif and as matlab figure
TifFile = fullfile(myFolder, 'Histogram.png');
FigFile = fullfile(myFolder, 'Histogram.fig');
saveas(gcf,TifFile);
saveas(gcf,FigFile);
savelog(myFolder,diam_avg,diam_std,conc_avg,conc_err,k,V_img,N);

function gui_threshold(theFiles,myFolder,scaling)
baseFileName = theFiles(1).name; % Define just the file name
fullFileName = fullfile(myFolder, baseFileName); %Create a full file address
img = imread(fullFileName);
gray = grayscale_process(img);
global bw_thresh circle_sens
binarythresh(gray);
uiwait;
bw = bw_process(gray,bw_thresh);
circle_thresh(bw,scaling);
uiwait;
end
function gray_img = grayscale_process(img)
%   gray       = rgb2gray(img); %convert to grayscale
  gray = img;
  filt       = medfilt2(gray,[2 2]); %Apply a grayscale contrast filter
  gray_img   = imgaussfilt(filt); %Apply a 2nd contrast filter
end
function bw = bw_process(gray,thresh)
    bw          = imbinarize(gray,'adaptive','ForegroundPolarity','dark','Sensitivity',thresh); %Convert to binary
    bw          = imcomplement(bw); %Invert the image eg. Black/white background/foreground to white/black
    bw          = imfill(bw,'holes'); %Fill holes in the image - eg. fill in the centre of the bubbles
end
function [scaling,dilution,BinWidth,spacer_height]=calibration()
%scaling value for microscope in microns per pixel
prompt     = {'Enter Image Scale (um/px):'}; %Create Prompt to enter and save value
dlg_title  = 'Input';
num_lines  = 1;
defaultans = {'0.19'};
scaling    = inputdlg(prompt,dlg_title,num_lines,defaultans); %Save user inputted value
scaling    = str2double(scaling(:)); %Convert string to double (number)
%%
prompt     = {'Enter Sample Dilution (e.g. 10)'};
dlg_title  = 'Input';
num_lines  = 1;
defaultans = {'10'};
dilution   = inputdlg(prompt,dlg_title,num_lines,defaultans);
dilution   = str2double(dilution(:));

prompt     = {'Enter Histogram Bin Width(um)'};
dlg_title  = 'Input';
num_lines  = 1;
defaultans = {'0.25'};
BinWidth   = inputdlg(prompt,dlg_title,num_lines,defaultans);
BinWidth   = str2double(BinWidth(:));

prompt     = {'Enter Height of Spacer/Chamber (um)'};
dlg_title  = 'Input';
num_lines  = 1;
defaultans = {'50'};
spacer_height   = inputdlg(prompt,dlg_title,num_lines,defaultans);
spacer_height   = str2double(spacer_height(:));
end
function bw_thresh = binarythresh(gray)
    figure('Visible','off');
    ax  = axes('Units','pixels');
    global binarysens
    binarysens = 0.5;
   figure(); set(gcf, 'Position', [100, 100, 1000, 500])
   
     %Create slider
     sld = uicontrol('Style', 'slider','Min',0,'Max',1,'Value',0.5,...
         'Position', [250 20 500 20],'Callback', @thresh1);
     btn = uicontrol('Style', 'pushbutton', 'String', num2str(binarysens),...
            'Position', [100 20 100 40]);
     bw_ini = bw_process(gray,0.5);
%    imshowpair(bw_ini,gray,'montage');
%       montage({bw_ini, gray},'DisplayRange',[]);
     subplot(1,2,1);imshow(bw_ini); subplot(1,2,2); imshow(gray,[]);
     
     
     

    function thresh1(source,event)

        btn = uicontrol('Style', 'pushbutton', 'String', num2str(binarysens),...
            'Position', [100 20 100 40],...
            'Callback', @pushbutton_callback);
        binarysens = source.Value;
        bw = bw_process(gray,binarysens);
%         imshowpair(bw,gray,'montage');
%         montage({bw, gray},'DisplayRange',[]);
        subplot(1,2,1);imshow(bw); subplot(1,2,2); imshow(gray,[]);

    end
     function pushbutton_callback(src,event)
         a = binarysens;
        assignin('base', 'bw_thresh', a)
        close
     end
   
         
 end
function circle_thresh(bw,scaling)
    ax = axes('Units','pixels');
    global sens
    sens = 0.8;
    figure();
    set(gcf, 'Position', [100, 100, 1000, 500])

   % Create slider
    sld = uicontrol('Style', 'slider', 'Min',0,'Max',1,'Value',0.8,...
        'Position', [250 20 500 20],'Callback', @thresh1); 
    btn = uicontrol('Style', 'pushbutton', 'String', num2str(sens),...
        'Position', [100 20 100 40]);   
    % Make figure visble after adding all components
    f.Visible = 'on';
    [centers,radii_px] = analysis_thresh(bw,0.8,scaling);
    imshow(bw)
    h = viscircles(centers,radii_px);
    
    function thresh1(source,event)
        sens = source.Value;
        [centers,radii_px] = analysis_thresh(bw,sens,scaling);
        imshow(bw)
         h = viscircles(centers,radii_px);
       btn = uicontrol('Style', 'pushbutton', 'String', num2str(sens),...
        'Position', [100 20 100 40],...
        'Callback', @pushbutton_callback);   
        
    end
     function pushbutton_callback(src,event)
         a = sens
        assignin('base', 'circle_sens', a)
        close
     end
 end
 function [centers,radii_px] = analysis_thresh(bw,sens,scaling)
 id  = 'images:imfindcircles:warnForSmallRadius';
 warning('off', id);
 id2 = 'images:imfindcircles:warnForLargeRadiusRange';
 warning('off',id2);
 [centers, radii_px , metric] = imfindcircles(bw,[2 round(20/scaling)],'ObjectPolarity','bright','Sensitivity',sens);
 diam_px   = 2*radii_px; %convert radius to diameter
 % Collate data centre and diameter of each circle/bubble into one big
 % array
 data      = cat(2,centers,diam_px);
 % This is some data filtering - loop throguh all the detected circles in
 % the image and if either the diameter is too small for the function or if it
 % isn't suitably circular, remove it. NaN =  Not a Number , basically a
 % blank space.
 for p = 1 : size(data,1)
     if data(p,3) <= 6 || metric(p) < 0.2
         data(p,:) = [NaN];
     end
 end
 %If any of the values above were changed to NaN, remove that circle from
 %the list by making it an empty vector eg. []
 data(~any(~isnan(data), 2),:) = [];
 
 % imfindcircles function struggles with smaller bubbles, so we use another
 % function regionprops to find smaller bubbles
 stats  = regionprops('table',bw,'Centroid',...
     'MajorAxisLength','MinorAxisLength','Eccentricity','Area','Perimeter');
 eccentricity     = stats.Eccentricity; %How Circular it is
 area_region      = stats.Area; %Area of circle
 Perimeter_region = stats.Perimeter; %Perimeter of circle
 circularity      = 4*pi*area_region ./ (Perimeter_region.^2); %Circularity
 centers_region   = stats.Centroid; %Centre of circle
 diam_px_region   = mean([stats.MajorAxisLength stats.MinorAxisLength],2); %Diameter of circle
 radii_px_region  = diam_px_region/2; %Radius of circle
 data_region      = cat(2,centers_region,diam_px_region); %Again, make big array of data
 for p=1:size(data_region,1)  %applying similar filtering as before however, getting rid of any big circles
     if data_region(p,3) > 6 || eccentricity(p) > 0.5
         data_region(p,:) = [NaN];
     end
 end
 %Same as before, removing circles from entry if they meet any of the
 %criteria
 data_region(~any(~isnan(data_region), 2),:) = [];
 % Now combine our two previous data vectors into one
 data = vertcat(data,data_region);
 data =  sortrows(data,3);
 
 %More Data Filtering - This section should remove any bubbles that
 %overlap with each other i.e, their centres are positioned within another
 %circle
 [distance,minDistance,minIndex] = deal(cell(1,size(data,1)));
 for p = 1:size(data,1)
     for pp = 1:size(data,1)
         distance{p}(pp,1)          = sqrt( (data(p,1) - data(pp,1)).^2 + (data(p,2) - data(pp,2)).^2 );
         distance{p}(distance{p} == 0)=NaN;
         [minDistance{p} ,minIndex{p}]= min(distance{p});
     end
 end
 % Create new variables for minDistance and minIndex
 data(:,4) = cell2mat(minDistance);
 data(:,5) = cell2mat(minIndex);
 
 for j = 1:size(data,1)
     % Remove the offending bubbles
     if data(j,4) <   0.5*data(data(j,5),3) || data(j,3) <= 4
         data(j,:) = [NaN];
    
     end
 end
 centers = [data(:,1) data(:,2)];
 radii_px = data(:,3)/2;
 end
 function [myFolder,filePattern,theFiles] = findfolder()
 myFolder=uigetdir; %open gui window to choose folder
 cd(myFolder); %Change directory to this folder, to load in and save files here automatically
 % % Check to see if the folder exists, if not then warn the user
 if ~isdir(myFolder)
     errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
     uiwait(warndlg(errorMessage));
     return;
 end
 % % Get a list of all files in the folder with the desired file name pattern.
 filePattern = fullfile(myFolder, '*.tif'); % Change to whatever pattern you need.
 theFiles    = dir(filePattern);
 end
 function savelog(myFolder,diam_avg,diam_std,conc_avg,conc_err,k,V_img,N)
global circle_sens bw_thresh
%Save a log file with meta data
LogFile = fullfile(myFolder,'log.txt');
fileID  = fopen(LogFile,'w');
fprintf(fileID,'%s %f ','Circle Threshold = ',circle_sens ); 
fprintf(fileID,'%s %f' ,' Binary Threshold = ', bw_thresh);
fprintf(fileID,'%s %f' , ' Mean Diameter = ', diam_avg);
fprintf(fileID,'%s %f' , ' Diameter Err = ', diam_std);
fprintf(fileID,'%s %e' , ' Concentration = ', conc_avg);
fprintf(fileID,'%s %e' , ' Concentration Err = ', conc_err);
fprintf(fileID,'%s %f' , ' Number of Images = ', k);
fprintf(fileID,'%s %e' , ' Image Volume (mL) = ', V_img);
fprintf(fileID,'%s %f' , ' Total Bubble Count = ', N);
fclose(fileID);
save('Workspace'); %Save workspace as a .mat file.
end
 