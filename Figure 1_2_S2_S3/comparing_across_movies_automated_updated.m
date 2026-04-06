% Comparing expression parameters output across movies - requires input from
% "spacer_analysis_code.m for each movie
close all
clear all

thr = 400; 
normMethod = "subtractMean"; % choose from "subtractMean", "subtractMean_divideStd", "divideMean", or "subtractMean_multiplyMean"  and make sure the appropriate threshold data has been generated from spacers_analysis_code_automated.mat


%% import filenames and gastrulation line details from excel file
XL_file = "C:\Users\Emilia Leyes Porello\OneDrive - PennO365\Penn Research\Enhancer-Promoter Spacers\Fly imaging\spacer analysis.xlsx";
sheet = "spacer analysis";
[~, txtData] = xlsread(XL_file,sheet);
[numData] = xlsread(XL_file,sheet);

rangeIDX = [13:15, 91:93];%[49:54, 64:66, 82:84];%[49:54, 55:57, 64:66];%[1:21, 25:27, 34:36, 43:45];%[1:3 10:15 19:21 25:27 34:36] ; %[10:15 19:21 25:27 34:36 43:45]; %[10:15]; %[1:21, 25:45];%[1:21];%[64:69];%[64:68];%[10 11 12 13 14 15]; %[1:24]; %all LacZ; %[28:33 37:42] %E1 & E2 %[10:15 25:27 34:36] t.yel vs. LacZ ; % enter embryo data to evaluate
emb_count = 3; % enter how many embryos evaluated per construct
nuc_diameter_est = 16; % approximate size of nuclei in NC14 (pixel count from imageJ -- approx. 4.75um)

% Automated construct name and count generation
names = cell(txtData(rangeIDX+1,5)); %names = names(~cellfun('isempty',names));
construct_count = length(names); 
% construct_names = names(1:construct_count)';
construct_names = cell(txtData(rangeIDX+1,3)); construct_names = construct_names(~cellfun('isempty',construct_names));

% Automated construct length generation
lengths = numData(rangeIDX,4); lengths = lengths(~isnan(lengths));
construct_lengths = lengths';

% Manual construct name and length generation - COMMENT OUT UNLESS NECESSARY
% construct_names = {'0kb', '1.3kb', '6.3kb', '7.6kb', '7.7kb', '8.8kb', '10kb', 'NoEnh', '7.7kb t.yel', '10kb t.yel'}; % manual input!!
% construct_lengths = [0 1.3 6.3 7.6 7.7 8.8 10 15 7.7 10]; % 15 as  placeholder for no enhancer

% start loop to iterate through all files
for k =rangeIDX %1:(size(txtData,1)-1)
    load_name = append(string(txtData(k+1,6)), "thr=", int2str(thr), '-', normMethod, "\analysis_output.mat");
    load(load_name)
    disp(append("loading data: ", string(txtData(k+1,5))))
    M2_list{k} = M2; % normalized ms2 trajectories (active and inactive) - for heatmap
    outputM_list{k} = outputM; % CHANGE
    output_lengths(k) = length(outputM);
    active_output_lengths(k) = length(active_outputM);
    active_rotated_list{k} = active_rotated; 
    frac_active_perFrame_list{k} = frac_active_perFrame;
    cumulative_active_count_list{k} = cumulative_active_count;
    total_active_count_list{k} = total_active_count;
    active_frac_perNuc_list{k} = active_frac_perNuc;
    active_frac_perNuc_sinceOnset_list{k} = active_frac_perNuc_sinceOnset;
    active_list{k} = active;
    active_outputM_list{k} = active_outputM; 
    gastrulation_line_list{k} = gastrulation_line;
    gastrulation_position_list{k} = gastrulation_position;
    top_width_list{k} = top_width;
    bottom_width_list{k} = bottom_width; 
    boundary_width_list{k} = boundary_width;
    onset_NC14prog_list{k} = onset_NC14prog;
    mean_amplitude_active_list{k} = mean_amplitude_active;
    max_amplitude_active_list{k} = max_amplitude_active;
    cumulative_outputM_perNuc_list{k} = cumulative_outputM_perNuc;
    cumulative_outputM_perNuc_active_list{k} = cumulative_outputM_perNuc_active;
    cumulative_outputM_perFrame_list{k} = cumulative_outputM_perFrame;
    inactive_to_active_count_perNuc_list{k} = inactive_to_active_count_perNuc;
    active_to_inactive_count_perNuc_list{k} = active_to_inactive_count_perNuc;
    burst_duration_list{k} = burst_duration;
    lineage_cy_list{k} = lineage_cy;
    lineage_cx_list{k} = lineage_cx;
    numChangesTo1_cell{k} = numChangesTo1;
    posFrames_cell{k} = posFrames;
    dist_posFrames_cell{k} = dist_posFrames_list;
end

max_nuc = max(output_lengths);
max_nuc_active = max(active_output_lengths);

% input names of each construct in order of loading
% names = {'10kb-emb01' '10kb-emb02' '10kb-emb03'...
%     '8.8kb-emb01' '8.8kb-emb02' '8.8kb-emb03'...
%     '7.7kb-emb01' '7.7kb-emb02' '7.7kb-emb03'...
%     '6.3kb-emb01' '6.3kb-emb02' '6.3kb-emb03'...
%     '1.3kb-emb01' '1.3kb-emb02' '1.3kb-emb03'};

% names = txtData(2:k+1,5);
% names = txtData(rangeIDX,5);
% names = cell(names);

% color_codes = [0, .6, .77; 0, .5, 0;.5, 0, 0; 0, 0, .5;.2, .2, .2];
color_codes = {'r', 'g', 'b', 'm', 'k', 'c', 'y'};

%% combining cells into matrices for relevant variables

mRNA_active_output_Matrix = nan(max_nuc_active,length(names));
for i = rangeIDX %1:length(names)
    mRNA_active_output_Matrix(1:length(active_outputM_list{i}),find(rangeIDX==i)) = active_outputM_list{i};
end

frac_time_active_Matrix = nan(max_nuc, length(names));
for i = rangeIDX %1:length(names)
    frac_time_active_Matrix(1:length(active_frac_perNuc_list{i}),find(rangeIDX==i)) = active_frac_perNuc_list{i};
end

frac_time_active_sinceOnset_Matrix = nan(max_nuc, length(names));
for i = rangeIDX %1:length(names)
    frac_time_active_sinceOnset_Matrix(1:length(active_frac_perNuc_sinceOnset_list{i}),find(rangeIDX==i)) = active_frac_perNuc_sinceOnset_list{i}; %#ok<*FNDSB> 
end

onset_NC14prog_Matrix = nan(max_nuc, length(names));
for i = rangeIDX %1:length(names)
    onset_NC14prog_Matrix(1:length(onset_NC14prog_list{i}),find(rangeIDX==i)) = onset_NC14prog_list{i};
end

onset_NC14prog_combinedPerConstruct = cell(1,construct_count);
for i = 1:length(names)
    pp = ceil(i/3);
    onset_NC14prog_combinedPerConstruct{pp} = [onset_NC14prog_combinedPerConstruct{pp}; onset_NC14prog_Matrix(:,i)];
    onset_NC14prog_combinedPerConstruct{pp} = onset_NC14prog_combinedPerConstruct{pp}(~isnan(onset_NC14prog_combinedPerConstruct{pp}));
end

mean_amplitude_active_Matrix = nan(max_nuc, length(names));
for i = rangeIDX %1:length(names)
    mean_amplitude_active_Matrix(1:length(mean_amplitude_active_list{i}),find(rangeIDX==i)) = mean_amplitude_active_list{i};
end

total_active_count_Matrix = nan(1, length(names));
for i = rangeIDX %1:length(names)
    total_active_count_Matrix(find(rangeIDX==i)) = total_active_count_list{i};
end

max_amplitude_active_Matrix = nan(max_nuc, length(names));
for i = rangeIDX %1:length(names)
    max_amplitude_active_Matrix(1:length(max_amplitude_active_list{i}),find(rangeIDX==i)) = max_amplitude_active_list{i};
end

for i = 1:length(names)
    pp = ceil(i/3);
    names_list(1:max_nuc_active,i) = string(construct_names{pp}); 
end

for i = 1:length(names)
    cy_active_list{i} = mean(lineage_cy_list{rangeIDX(i)}(:,active_list{rangeIDX(i)}));
    cy_active_normalized{i} = cy_active_list{i} - gastrulation_line_list{rangeIDX(i)};
end

mRNA_active_output_combinedPerConstruct = cell(1,construct_count);
for i = 1:length(names)
    pp = ceil(i/3);
    mRNA_active_output_combinedPerConstruct{pp} = [mRNA_active_output_combinedPerConstruct{pp}; mRNA_active_output_Matrix(:,i)];
    mRNA_active_output_combinedPerConstruct{pp} = mRNA_active_output_combinedPerConstruct{pp}(~isnan(mRNA_active_output_combinedPerConstruct{pp}));
end

% removing blanks from active nuclei list
active_list(cellfun('isempty', active_list))=[];

numChangesTo1_Matrix = nan(max_nuc, length(names));
for i = rangeIDX %1:length(names)
    numChangesTo1_Matrix(1:length(numChangesTo1_cell{i}),find(rangeIDX==i)) = numChangesTo1_cell{i};
end

posFrames_Matrix = nan(max(find_max_dimensions(posFrames_cell)), length(names));
for i = rangeIDX %1:length(names)
     posFrames_Matrix(1:length(posFrames_cell{i}),find(rangeIDX==i)) =  posFrames_cell{i};
end

dist_posFrames_Matrix = nan(max(find_max_dimensions(dist_posFrames_cell)), length(names));
for i = rangeIDX %1:length(names)
    dist_posFrames_cell{i} = dist_posFrames_cell{i}(~isnan(dist_posFrames_cell{i}));
    dist_posFrames_Matrix(1:length(dist_posFrames_cell{i}),find(rangeIDX==i)) = dist_posFrames_cell{i};
end

%% generate box plot of mRNA output for each video independently
figure;
boxplot(mRNA_active_output_Matrix)
set(gca,'XTick',1:size(onset_NC14prog_list,2),'XTickLabel',names)
ylabel('mRNA output (A.U)')
title(normMethod, 'Interpreter', 'none')

%% generate box plot of onset times for each video independently
figure;
boxplot(onset_NC14prog_Matrix)
set(gca,'XTick',1:size(onset_NC14prog_list,2),'XTickLabel',names)
ylabel('onset time (prog through NC14)')
title(normMethod, 'Interpreter', 'none')
 
%% generate box plot of active frac since onset for each video independently
figure;
boxplot(mean_amplitude_active_Matrix)
set(gca,'XTick',1:size(onset_NC14prog_list,2),'XTickLabel',names)
ylabel('mean amplitude of active frames')
title(normMethod, 'Interpreter', 'none')

%% generate box plot of mean amplitude for each video independently
figure;
boxplot(frac_time_active_sinceOnset_Matrix)
set(gca,'XTick',1:size(onset_NC14prog_list,2),'XTickLabel',names)
ylabel('frac time active since onset')
title(normMethod, 'Interpreter', 'none')

%% generate box plot of mRNA output for each construct
% figure;
% mRNA_active_output_combined = reshape(mRNA_active_output_Matrix,[],1);
% names_list_combined = reshape(names_list,[],1);
% boxplot(mRNA_active_output_combined, names_list_combined)
% set(gca,'XTick',1:length(construct_names),'XTickLabel',construct_names)
% ylabel('Total mRNA production (AU)')
% title(normMethod, 'Interpreter', 'none')
% 
% % combine mRNA output data by construct THEN remove top and bottomr 10% each construct
% mRNA_active_output_Matrix_innerPercent = nan(size(names_list,1)*3,length(construct_names));
% % remove top and bottom 10% of data points
% for i = 1:length(construct_names)
%     X = mRNA_active_output_combinedPerConstruct{i}; % X - list of mRNA output for data set i
%     percentile = prctile(X,[10 90]);
%     X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
%     mRNA_active_output_Matrix_innerPercent(1:length(X),i) = X;
% end
% 
% 
% % % remove 10% from each embryo data, before combining
% % mRNA_active_output_Matrix_innerPercent = nan(size(mRNA_active_output_Matrix));
% % % remove top and bottom 10% of data points
% % for i = 1:length(construct_names)
% %     X = mRNA_active_output_Matrix(:,i); % X - list of mRNA output for data set i
% %     percentile = prctile(X,[10 90]);
% %     X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
% %     mRNA_active_output_Matrix_innerPercent(:,i) = X;
% % end
% 
% % plot inner 80 percentile mRNA output
% figure;
% X_combined = reshape(mRNA_active_output_Matrix_innerPercent, [], 1);
% boxplot(X_combined, names_list_combined)
% set(gca,'XTick',1:length(construct_names),'XTickLabel',construct_names)
% ylabel('Total mRNA production (AU)')
% title(normMethod, 'Interpreter', 'none')


%% Stop auto-run code - manually run sections from here on

return


%% Generating mRNA output box plots using *Box Chart*

toPlot = [1 2 ]%1 2 4];%[1:8]; % construct idx to plot

% combine mRNA output data by construct THEN remove top and bottom 10% each construct
mRNA_active_output_Matrix_innerPercent = nan(size(names_list,1)*3,length(construct_names(toPlot)));
% remove top and bottom 10% of data points
count = 1; 
for i = toPlot%1:length(construct_names(toPlot))
    X = mRNA_active_output_combinedPerConstruct{i}; % X - list of mRNA output for data set i
    percentile = prctile(X,[10 90]);
    X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
    mRNA_active_output_Matrix_innerPercent(1:length(X),count) = X;
    count = count + 1;
end

xlabel_positions =[1 2 ];% [ 0    1.3   6.3    7.3    7.9    8.8   10 ]; %[1 2 ];% % manually set position of boxcharts

xStep = nan(size(mRNA_active_output_Matrix_innerPercent));
for i = 1:length(toPlot)
%     xStep(:,i) = construct_lengths(toPlot(i));
    xStep(:,i) = xlabel_positions(i);
end

%xStep_combined = reshape(xStep(:,toPlot([1 2 3 4])),[],1);
xStep_combined = reshape(xStep(:,[1 2 ]),[],1);
%output_combined = reshape(mRNA_active_output_Matrix_innerPercent(:,toPlot([1 2 3 4])), [], 1);
output_combined = reshape(mRNA_active_output_Matrix_innerPercent(:,[1 2 ]), [], 1);

figure;
swarmchart(xStep_combined, output_combined, 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on; boxchart(xStep_combined,output_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','b', 'BoxFaceColor','b', 'BoxWidth', 0.7);
xticks(xStep(1,1:length(toPlot)))
xticklabels(construct_names(toPlot))
% xticklabels(construct_lengths(toPlot([1 2 3 4 ])))
% xticklabels({'E1' 'E2' 'E1' 'E2'})
xtickangle(45)
%hold on; plot(xStep(1,[1 2 3 4 6]),nanmean(mRNA_active_output_Matrix_innerPercent(:,sort([1 2 3 4 6]))),'--o', 'Color', 'k', 'LineWidth',1.2)
ylabel('Mean mRNA production (AU)')
set(gca,'fontname','gillsans')
% pbaspect([1 1 1])


%% mRNA output Boxchart in different colors - run this section right after the previous one to overlay additional box charts (in diff colors) to the same plot

a1 = 2;
a2 = 3;
a3 = 4; 

hold on;
swarmchart(xStep(:, a1),mRNA_active_output_Matrix_innerPercent(:, a1), 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(xStep(:, a1),mRNA_active_output_Matrix_innerPercent(:, a1), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);
hold on;
swarmchart(xStep(:, a2),mRNA_active_output_Matrix_innerPercent(:, a2), 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(xStep(:, a2),mRNA_active_output_Matrix_innerPercent(:, a2), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);
hold on;
swarmchart(xStep(:, a3),mRNA_active_output_Matrix_innerPercent(:, a3), 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(xStep(:, a3),mRNA_active_output_Matrix_innerPercent(:, a3), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);


% xticks(xStep(1,[5 7]))
% xticklabels(construct_names([5 7]))
% xticklabels(construct_lengths([5 7]))
% hold on; plot(xStep(1,[5 7]),nanmean(mRNA_active_output_Matrix_innerPercent(:,[5 7])),'--o', 'Color', 'k', 'LineWidth',1.2)
% ylabel('Mean mRNA production (AU)')
% set(gca,'fontname','gillsans')
% xlim([7.52 7.78])

%% calculating active nuc # for E1/E2 data
% aa = mRNA_active_output_Matrix_innerPercent;
% for i = 1:4
% aa = a(:,i);
% B = aa(~isnan(aa));
% count(i) = length(B);
% end


%% Boxchart for onset times
% onset times - plotting all constructs combined (inner 80 percentile only)

toPlot = [1:2]; % construct idx to plot

% combine onset time data by construct THEN remove top and bottom 10% each construct
onset_NC14prog_perConstruct_innerPercent = nan(size(names_list,1)*3,length(construct_names(toPlot)));
% remove top and bottom 10% of data points
count = 1;
for i = toPlot%1:length(construct_names)
    X = onset_NC14prog_combinedPerConstruct{i}; % X - list of mRNA output for data set i
    percentile = prctile(X,[10 90]);
    X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
    onset_NC14prog_perConstruct_innerPercent(1:length(X),count) = X;
    count = count + 1;
end

xlabel_positions = [1 2  ];% [ 0    1.3   6.3    7.3    7.9    8.8   10 ];%[ 0    1.3   6.3    7.3    8.0    8.8   10   15]; % manually set position of boxcharts

xStep = nan(size(mRNA_active_output_Matrix_innerPercent));
for i = 1:length(construct_lengths(toPlot))
%     xStep(:,i) = construct_lengths(i);
    xStep(:,i) = xlabel_positions(i);
end

% xStep_combined = reshape(xStep(:,toPlot([1 ])),[],1);
% onset_combined = reshape(onset_NC14prog_perConstruct_innerPercent(:,toPlot([1 ])), [], 1);
xStep_combined = reshape(xStep(:,[1:2]),[],1);
onset_combined = reshape(onset_NC14prog_perConstruct_innerPercent(:,[1:2 ]), [], 1);

figure;

swarmchart(xStep_combined, onset_combined, 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on; boxchart(xStep_combined, onset_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','b', 'BoxWidth', 0.7);
xtickangle(45)
xticks(xStep(1,:))
xticklabels(construct_names(toPlot))
% xticklabels(construct_lengths(toPlot))
% xticklabels({'E1' 'E2' 'E1' 'E2'})
%hold on; plot(xStep(1,[1 2 3 4 6 ]),nanmean(onset_NC14prog_perConstruct_innerPercent(:,[1 2 3 4 6])),'--o', 'Color', 'k', 'LineWidth',1.2)
ylabel('Onset time (frac through NC14)')
set(gca,'fontname','gillsans')
ylim([0 1])
% pbaspect([1 1 1])

%% Onset Boxchart in different colors- run this section right after the previous one to overlay additional box charts (in diff colors) to the same plot
a1=2;
a2=3;
a3=4;

% figure;
hold on;
swarmchart(xStep(:, a1),onset_NC14prog_perConstruct_innerPercent(:, a1), 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(xStep(:, a1),onset_NC14prog_perConstruct_innerPercent(:, a1), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);
hold on;
swarmchart(xStep(:, a2),onset_NC14prog_perConstruct_innerPercent(:, a2), 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(xStep(:, a2),onset_NC14prog_perConstruct_innerPercent(:, a2), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);
hold on;
swarmchart(xStep(:, a3),onset_NC14prog_perConstruct_innerPercent(:, a3), 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(xStep(:, a3),onset_NC14prog_perConstruct_innerPercent(:, a3), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);


% xticks(xStep(1,[5 7]))
% xticklabels(construct_names([5 7]))
% % hold on; plot(xStep(1,[5 7]),nanmean(onset_NC14prog_perConstruct_innerPercent(:,[5 7])),'--o', 'Color', 'k', 'LineWidth',1.2)
% ylabel('Onset time (progress through NC14)')
% set(gca,'fontname','gillsans')
% xlim([7.55 7.75])

%% Comparing LacZ to trunc Yellow
% mRNA output
xStep = nan(size(mRNA_active_output_Matrix_innerPercent));
for i = 1:length(construct_lengths(toPlot))
    xStep(:,i) = construct_lengths(i);
%     xStep(:,i) = xlabel_positions(i);
end

xStep_manual = ones(length(xStep),2);
xStep_manual(:,2) = 2;
xStep_combined = reshape(xStep_manual,[],1);
output_combined = reshape(mRNA_active_output_Matrix_innerPercent(:,[1 3]), [], 1);


figure; 
swarmchart(xStep_combined,output_combined, 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on; 
boxchart(xStep_combined,output_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','b', 'BoxWidth', 0.5);
xtickangle(45)

% xticks([1 2])
% xticklabels(construct_names([1 3]))
% xtickangle(45)
% ylabel('mRNA output (A.U)')

xStep_manual = ones(length(xStep),2);
xStep_manual(:,1) = 4;
xStep_manual(:,2) = 5;
xStep_combined = reshape(xStep_manual,[],1);
output_combined = reshape(mRNA_active_output_Matrix_innerPercent(:,[2 4]), [], 1);

hold on;
swarmchart(xStep_combined,output_combined, 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(xStep_combined,output_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.5);
xticks([1 2 4 5])
xticklabels(construct_names([1 3 2 4]))
xtickangle(45)
ylabel('mRNA output (A.U)')
% legend({'LacZ (5 CTCF)' 't.Yellow (0 CTCF)' })


%onset time
xStep_manual = ones(length(xStep),2);
xStep_manual(:,2) = 2;
xStep_combined = reshape(xStep_manual,[],1);
output_combined = reshape(onset_NC14prog_perConstruct_innerPercent(:,[1 3]), [], 1);

figure;
swarmchart(xStep_combined,output_combined, 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(xStep_combined,output_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','b', 'BoxWidth', 0.5);
xtickangle(45)

xStep_manual = ones(length(xStep),2);
xStep_manual(:,1) = 4;
xStep_manual(:,2) = 5;
xStep_combined = reshape(xStep_manual,[],1);
output_combined = reshape(onset_NC14prog_perConstruct_innerPercent(:,[2 4]), [], 1);

hold on; 
swarmchart(xStep_combined,output_combined, 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(xStep_combined,output_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.5);
xticks([1 2 4 5])
xticklabels(construct_names([1 3 2 4]))
xtickangle(45)
ylabel('onset time (frac progress NC14)')
% legend({'LacZ (5 CTCF)' 't.Yellow (0 CTCF)' })


% fraction active after onset - have to run next section first
xStep_manual = ones(length(frac_time_active_sinceOnset_perConstruct),2);
xStep_manual(:,2) = 2;
xStep_combined = reshape(xStep_manual,[],1);
output_combined = reshape(frac_time_active_sinceOnset_perConstruct(:,[1 3]), [], 1);

figure;
swarmchart(xStep_combined,output_combined, 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(xStep_combined,output_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','b', 'BoxWidth', 0.5);
xtickangle(45)

xStep_manual = ones(length(frac_time_active_sinceOnset_perConstruct),2);
xStep_manual(:,1) = 4;
xStep_manual(:,2) = 5;
xStep_combined = reshape(xStep_manual,[],1);
output_combined = reshape(frac_time_active_sinceOnset_perConstruct(:,[2 4]), [], 1);

hold on; 
swarmchart(xStep_combined,output_combined, 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(xStep_combined,output_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.5);
xticks([1 2 4 5])
xticklabels(construct_names([1 3 2 4]))
xtickangle(45)
ylabel('Frac NC14 spent in ACTIVE state SINCE ONSET')
% legend({'LacZ (5 CTCF)' 't.Yellow (0 CTCF)' })

% Mean amplitude
xStep_manual = ones(length(mean_amplitude_active_perConstruct),2);
xStep_manual(:,2) = 2;
xStep_combined = reshape(xStep_manual,[],1);
output_combined = reshape(mean_amplitude_active_perConstruct(:,[1 3]), [], 1);

figure;
swarmchart(xStep_combined,output_combined, 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(xStep_combined,output_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','b', 'BoxWidth', 0.5);
xtickangle(45)

xStep_manual = ones(length(mean_amplitude_active_perConstruct),2);
xStep_manual(:,1) = 4;
xStep_manual(:,2) = 5;
xStep_combined = reshape(xStep_manual,[],1);
output_combined = reshape(mean_amplitude_active_perConstruct(:,[2 4]), [], 1);

hold on; 
swarmchart(xStep_combined,output_combined, 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(xStep_combined,output_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.5);
xticks([1 2 4 5])
xticklabels(construct_names([1 3 2 4]))
xtickangle(45)
ylabel('Mean Transcriptional Amplitude (AU)')
% legend({'LacZ (5 CTCF)' 't.Yellow (0 CTCF)' })

%% Box chart to compare time spent in active state (SINCE ONSET) per construct
clear frac_time_active_sinceOnset_perConstruct
clear frac_active_xGroup_combined
clear frac_active_ONSET_combined


xlabel_positions =  [1:2]; %[ 0    1.3   6.3    7.3    7.9    8.8   10 ];%construct_lengths;;%[ 0    1.3   6.3    7.3    8.0    8.8   10   15]; % manually set position of boxcharts

toEval = [1:2]; % list of construct indices to evaluate
count = 1;
for i = toEval
    pp = [1:3]+((i-1)*3);
    X = reshape(frac_time_active_sinceOnset_Matrix(:,pp),[],1);
    percentile = prctile(X,[10 90]); % can comment out this line and next to keep all data on plots
    X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
    frac_time_active_sinceOnset_perConstruct(1:3*(size(frac_time_active_sinceOnset_Matrix,1)),count) = X;
    frac_active_xGroup(1:1:3*(size(frac_time_active_sinceOnset_Matrix,1)),count) = xlabel_positions(count);%construct_lengths(count);
    count = count + 1;
end

toPlot = [1 2];

frac_active_xGroup_combined = reshape(frac_active_xGroup(:,toPlot([1:2])),[],1);
frac_active_ONSET_combined = reshape(frac_time_active_sinceOnset_perConstruct(:,toPlot([1:2])), [], 1);

figure; 
swarmchart(frac_active_xGroup_combined, frac_active_ONSET_combined, 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on; boxchart(frac_active_xGroup_combined, frac_active_ONSET_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','b', 'BoxWidth', 0.7);
xtickangle(45)
xticks(xStep(1,:))
% ticklabels(construct_names(sort([toPlot ])))
xticklabels(construct_names(toEval(toPlot)))
% xticklabels({'E1' 'E2' 'E1' 'E2'})
ylabel('NC14 fraction active (since onset)')
set(gca,'fontname','gillsans')


%% frac active since onset - Boxchart in different colors - run this section right after the previous one to overlay additional box charts (in diff colors) to the same plot

a1=2;
a2=3;
a3=4;

% figure;
hold on;
swarmchart(frac_active_xGroup(:, a1),frac_time_active_sinceOnset_perConstruct(:, a1), 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(frac_active_xGroup(:, a1),frac_time_active_sinceOnset_perConstruct(:, a1), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);
hold on;
swarmchart(frac_active_xGroup(:, a2),frac_time_active_sinceOnset_perConstruct(:, a2), 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(frac_active_xGroup(:, a2),frac_time_active_sinceOnset_perConstruct(:, a2), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);
hold on;
swarmchart(frac_active_xGroup(:, a3),frac_time_active_sinceOnset_perConstruct(:, a3), 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(frac_active_xGroup(:, a3),frac_time_active_sinceOnset_perConstruct(:, a3), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);
% 
% xticks(xStep(1,[5 7]))
% xticklabels(construct_names([5 7]))
% ylabel('Frac NC14 spent in ACTIVE state SINCE ONSET')
% set(gca,'fontname','gillsans')
% xlim([7.55 7.75])


%% mean amplitude - Box charts of mean amplitudes for active nuclei (only in their active frames)

% Uses same x spacing as above

% MEAN amplitudes - plotting all constructs combined (inner 80 percentile only)
clear mean_amplitude_active_perConstruct
clear meanAmp_xGroup_combined
clear meanAmp_combined

% xlabel_positions = [1 2] %construct_lengths;%[1 2 3 4];
count = 1;
toEval = [1:2]; % list of construct indices to evaluate
for i = toEval
    pp = [1:3]+((i-1)*3);
    X = reshape(mean_amplitude_active_Matrix(:,pp),[],1);
    percentile = prctile(X,[10 90]); % can comment out this line and next to keep all data on plots
    X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
    mean_amplitude_active_perConstruct(1:3*(size(mean_amplitude_active_Matrix,1)),count) = X;
    count = count + 1;
end

toPlot = [1 2 ];

meanAmp_xGroup_combined = reshape(frac_active_xGroup(:,toPlot([1:2 ])),[],1); % uses same matrix as onset time
meanAmp_combined = reshape(mean_amplitude_active_perConstruct(:,toPlot([1:2 ])), [], 1);

figure; 
swarmchart(meanAmp_xGroup_combined, meanAmp_combined, 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on; boxchart(meanAmp_xGroup_combined, meanAmp_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','b', 'BoxWidth', 0.7);
xtickangle(45)
xticks(xStep(1,:))
% xticks(xStep(1,sort([toPlot])))
% xticklabels(construct_names(sort([toPlot])))
xticklabels(construct_names(toEval(toPlot)))
% xticklabels({'E1' 'E2' 'E1' 'E2'})
ylabel('Mean Transcriptional Amplitude (AU)')
set(gca,'fontname','gillsans')

%% mean amplitude - Boxchart in different colors - run this section right after the previous one to overlay additional box charts (in diff colors) to the same plot

a1=2;
a2=3;
a3=4;

hold on;
swarmchart(frac_active_xGroup(:, a1),mean_amplitude_active_perConstruct(:, a1), 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(frac_active_xGroup(:, a1),mean_amplitude_active_perConstruct(:, a1), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);
hold on;
swarmchart(frac_active_xGroup(:, a2),mean_amplitude_active_perConstruct(:, a2), 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(frac_active_xGroup(:, a2),mean_amplitude_active_perConstruct(:, a2), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);
hold on;
swarmchart(frac_active_xGroup(:, a3),mean_amplitude_active_perConstruct(:, a3), 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on;
boxchart(frac_active_xGroup(:, a3),mean_amplitude_active_perConstruct(:, a3), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);


% xticks(xStep(1,[5 7]))
% xticklabels(construct_names([5 7]))
% ylabel('Mean Transcriptional Amplitude (AU)')
%set(gca,'fontname','gillsans')
%% Plotting bursting characterizations

% combine data per construct

numChangesTo1_combinedPerConstruct = nan(3*size(numChangesTo1_Matrix,1),ceil(length(names)/3));
count = 1;
for i = 1:3:length(names)
    combine = [numChangesTo1_Matrix(:,i); numChangesTo1_Matrix(:,i+1); numChangesTo1_Matrix(:,i+2)];
    numChangesTo1_combinedPerConstruct(:,count) = combine;
    count = count + 1; 
end

posFrames_combinedPerConstruct = nan(3*size(posFrames_Matrix,1),ceil(length(names)/3));
count = 1;
for i = 1:3:length(names)
    combine = [posFrames_Matrix(:,i); posFrames_Matrix(:,i+1); posFrames_Matrix(:,i+2)];
    posFrames_combinedPerConstruct(:,count) = combine;
    count = count + 1; 
end

dist_posFrames_combinedPerConstruct = nan(3*size(dist_posFrames_Matrix,1),ceil(length(names)/3));
count = 1;
for i = 1:3:length(names)
    combine = [dist_posFrames_Matrix(:,i); dist_posFrames_Matrix(:,i+1); dist_posFrames_Matrix(:,i+2)];
    dist_posFrames_combinedPerConstruct(:,count) = combine;
    count = count + 1; 
end

% % Remove zeros
% numChangesTo1_combinedPerConstruct(numChangesTo1_combinedPerConstruct==0) = NaN;
% posFrames_combinedPerConstruct(posFrames_combinedPerConstruct==0) = NaN;
% dist_posFrames_combinedPerConstruct(dist_posFrames_combinedPerConstruct==0) = NaN;

% test plot 
figure; boxplot(numChangesTo1_combinedPerConstruct)
title('numChangesTo1')
xticklabels(construct_names)

figure; boxplot(posFrames_combinedPerConstruct)
title('total positive frames')
xticklabels(construct_names)

figure; boxplot(dist_posFrames_combinedPerConstruct)
title('distribution of positive frames')
xticklabels(construct_names)

%PICK UP HERE!!! (1/17/25)
% plot # of bursting instances (# of pos slope instances, regardless of duration)
toPlot = [1:20];%[1:8]; % construct idx to plot

% combine mRNA output data by construct THEN remove top and bottomr 10% each construct
mRNA_active_output_Matrix_innerPercent = nan(size(names_list,1)*3,length(construct_names(toPlot)));
% remove top and bottom 10% of data points
count = 1; 
for i = toPlot%1:length(construct_names(toPlot))
    X = mRNA_active_output_combinedPerConstruct{i}; % X - list of mRNA output for data set i
    percentile = prctile(X,[10 90]);
    X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
    mRNA_active_output_Matrix_innerPercent(1:length(X),count) = X;
    count = count + 1;
end

xlabel_positions = toPlot; %[1 2 3 4];%[ 0    1.3   6.3    7.4    7.9    8.8   10   15]; % manually set position of boxcharts

xStep = nan(size(mRNA_active_output_Matrix_innerPercent));
for i = 1:length(toPlot)
%     xStep(:,i) = construct_lengths(toPlot(i));
    xStep(:,i) = xlabel_positions(i);
end

%xStep_combined = reshape(xStep(:,toPlot([1 2 3 4])),[],1);
xStep_combined = reshape(xStep(:,[1 2 3 4]),[],1);
%output_combined = reshape(mRNA_active_output_Matrix_innerPercent(:,toPlot([1 2 3 4])), [], 1);
output_combined = reshape(mRNA_active_output_Matrix_innerPercent(:,[1 2 3 4]), [], 1);

figure;
swarmchart(xStep_combined, output_combined, 5, 'o', 'MarkerEdgeColor', [.7 .7 .7])
hold on; boxchart(xStep_combined,output_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','b', 'BoxWidth', 0.7);
xticks(xStep(1,1:length(toPlot)))
xticklabels(construct_names(toPlot))
% xticklabels(construct_lengths(toPlot([1 2 3 4 ])))
% xticklabels({'E1' 'E2' 'E1' 'E2'})
xtickangle(45)
%hold on; plot(xStep(1,[1 2 3 4 6]),nanmean(mRNA_active_output_Matrix_innerPercent(:,sort([1 2 3 4 6]))),'--o', 'Color', 'k', 'LineWidth',1.2)
ylabel('Mean mRNA production (AU)')
set(gca,'fontname','gillsans')


%%  fraction active - Boxchart

clear frac_time_active_perConstruct
clear frac_active_xGroup_combined
clear frac_active_combined

count = 1;
toEval = [1:4]; % list of construct indices to evaluate
for i = toEval
    pp = [1:3]+((i-1)*3);
    X = reshape(frac_time_active_Matrix(:,pp),[],1);
    percentile = prctile(X,[10 90]); % can comment out this line and next to keep all data on plots
    X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
    frac_time_active_perConstruct(1:3*(size(frac_time_active_Matrix,1)),count) = X;
    count = count + 1;
end

toPlot = [1 2 3 4];

frac_active_xGroup_combined = reshape(frac_active_xGroup(:,toPlot),[],1);
frac_active_combined = reshape(frac_time_active_perConstruct(:,toPlot), [], 1);

figure; 
boxchart(frac_active_xGroup_combined, frac_active_combined, 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','b', 'BoxWidth', 0.7);
xtickangle(45)
xticks(xStep(1,sort([toPlot 5 7])))
xticklabels(construct_names(sort([toPlot 5 7])))
ylabel('Frac NC14 spent in ACTIVE state')
set(gca,'fontname','gillsans')

%% frac active - Boxchart in different colors - run this section right after the previous one to overlay additional box charts (in diff colors) to the same plot

% figure;
hold on;
boxchart(frac_active_xGroup(:, 5),frac_time_active_perConstruct(:, 5), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);
hold on;
boxchart(frac_active_xGroup(:, 7),frac_time_active_perConstruct(:, 7), 'Notch', 'on', 'MarkerStyle', '.', 'MarkerColor','k', 'BoxFaceColor','r', 'BoxWidth', 0.7);

% xticks(xStep(1,[5 7]))
% xticklabels(construct_names([5 7]))
% ylabel('Frac NC14 spent in ACTIVE state SINCE ONSET')
set(gca,'fontname','gillsans')

%% interpolation parameters of time series data
x_rawData = cell(1,length(names));
for i = rangeIDX%1:length(names)
    x_rawData{i} = linspace(0,1,length(cumulative_active_count_list{i}));
end
x_Int = linspace(0,1,180);

% cumulative active nuclei per frame across single constructs

%individual time course plots
figure;
count = 1; 
for i = rangeIDX% 1:length(names)
    pp = ceil(count/3);
    subplot(1, length(construct_names),pp); hold  on;
    y = interp1(x_rawData{i},cumulative_active_count_list{i}, x_Int);
   cumulative_active_count_list_interp{i} = y;
    plot(x_Int,y);
%     plot(linspace (0,1,length(frac_active_perFrame_list{i})),frac_active_perFrame_list{i})
    title(append( construct_names{pp}))
    legend('emb01', 'emb02', 'emb03', 'Location','northwest')
    xlabel('Frac Progress through NC14')
    ylabel('Cumulative # active nuclei')
    ylim([0 400])
    axis square
    count = count + 1; 
end

% time course plots with shaded error bars

% first we combine the data for each construct into cells
figure; sgtitle("Cumulative # active nuclei")
count = 1;
constructIDX = rangeIDX(1:3:end); % take every 3rd element
for i = constructIDX%1:length(construct_names) % # of constructs
%     i = ceil(count/3);
%     pp = (i-1)*3+1:(i-1)*3+3;
    pp = [i i+1 i+2];
    combine_cum_active_interp{count}(1,:) = cumulative_active_count_list_interp{pp(1)};
    combine_cum_active_interp{count}(2,:) = cumulative_active_count_list_interp{pp(2)};
    combine_cum_active_interp{count}(3,:) = cumulative_active_count_list_interp{pp(3)};
    mean_cum_active_interp{count} = mean(combine_cum_active_interp{count});
    stdev_cum_active_interp{count} = std(combine_cum_active_interp{count});
    up_error_cum_active_interp{count} =mean_cum_active_interp{count}+stdev_cum_active_interp{count};
    low_error_cum_active_interp{count} = mean_cum_active_interp{count}-stdev_cum_active_interp{count};

%     subplot(1, length(construct_names),count)
    h_line(count) = plot(x_Int,mean_cum_active_interp{count}, 'Color', string(color_codes(count)), 'LineWidth',1.5); hold on;
    patch([x_Int fliplr(x_Int)], smooth([low_error_cum_active_interp{count} fliplr(up_error_cum_active_interp{count})]), string(color_codes(count)), 'FaceAlpha',0.5, 'EdgeColor','none')
    %title(construct_names(count))
    xlabel("Progress through NC14")
    ylabel("Cumulative # active nuclei")
    ylim([0 500])
    %axis square
    alpha 0.3
    count = count + 1;
    hold on;
end
legend(h_line,construct_names)
hold off; 
%%
% Summary bar chart - not accounting for rotation
x = 1:length(constructIDX);
total_active = nan(1,length(x));
errhigh = nan(1,length(x));
errlow = nan(1,length(x));
for i = x
    total_active(i) = mean_cum_active_interp{i}(end);
    errhigh(i) = stdev_cum_active_interp{i}(end); % high and low error are same
    errlow(i) = stdev_cum_active_interp{i}(end);
end

up = [1 2 3 4 6];
down = [5 7];
figure; 
b = bar(x, total_active, 'FaceColor', 'b', 'BarWidth', 0.8);
b.FaceColor = 'flat';
b.CData(up,:) = ones(length(up),3).*[0 0 1];
b.CData(down,:) = ones(length(down),3).*[1 0 0];
b.CData(8,:) = [.7 .7 .7];
% b.CData(1,:) = 'r';
% hold on;
% bar(x(down), total_active(down), 'FaceColor', 'r', 'BarWidth', 0.8);
hold on;
er = errorbar(x,total_active,errlow,errhigh);    
er.Color = 'k';                            
er.LineStyle = 'none';  

xticklabels(construct_names)
ylabel('Total Active Nuclei')
title('Active Nuc# - NOT Accounting for Emb Rotation')


%%
% Summary bar chart - accounting for rotation
x = 1:length(constructIDX);
total_active_accountRot = nan(1,length(constructIDX));
errhigh_accountRot = nan(1,length(constructIDX));
errlow_accountRot = nan(1,length(constructIDX));
count = 1;
for i = constructIDX
    pp = [i i+1 i+2];
    total_active_accountRot(count) = mean(cell2mat(total_active_count_list(pp)));
    errhigh_accountRot(count) = std(cell2mat(total_active_count_list(pp)));
    errlow_accountRot(count) = errhigh_accountRot(count);
    count = count + 1;
end
up = [1 2 3 4 6];
down = [5 7];
figure; 
b = bar(x, total_active_accountRot, 'FaceColor', 'b', 'BarWidth', 0.8);
b.FaceColor = 'flat';
b.CData(up,:) = ones(length(up),3).*[0 0 1];
b.CData(down,:) = ones(length(down),3).*[1 0 0];
b.CData(8,:) = [.7 .7 .7];
% b.CData(1,:) = 'r';
% hold on;
% bar(x(down), total_active(down), 'FaceColor', 'r', 'BarWidth', 0.8);
hold on;
er = errorbar(x,total_active_accountRot,errlow_accountRot,errhigh_accountRot);    
er.Color = 'k';                            
er.LineStyle = 'none';  

xticklabels(construct_names)
ylabel('Total Active Nuclei')
title('Active Nuc# - Accounting for Emb Rotation')


%% Fraction active nuclei per frame
%individual plots
figure;
count = 1; 
for i = rangeIDX% 1:length(names)
    pp = ceil(count/3);
    subplot(1, length(construct_names),pp); hold  on;
    y = interp1(x_rawData{i},frac_active_perFrame_list{i}, x_Int);
    frac_active_perFrame_list_interp{i} = y;
    plot(x_Int,y);
%     plot(linspace (0,1,length(frac_active_perFrame_list{i})),frac_active_perFrame_list{i})
    title(append( construct_names{pp}))
    legend('emb01', 'emb02', 'emb03', 'Location','northwest')
    ylim([0 1])
    xlabel('Frac Progress through NC14')
    ylabel('Instantaneous frac active nuclei')
    axis square
    count = count + 1; 
end
%%
% shaded error bars
% color_codes = ['r' 'k'];
% first we combine the data for each construct into cells
figure; sgtitle("Fraction active nuclei per frame")
count = 1;
constructIDX = rangeIDX(1:3:end); % take every 3rd element
for i = constructIDX%1:length(construct_names) % # of constructs
%     i = ceil(count/3);
    pp = (i-1)*3+1:(i-1)*3+3;
    pp = [i i+1 i+2];
    combine_frac_active_interp{count}(1,:) = frac_active_perFrame_list_interp{pp(1)};
    combine_frac_active_interp{count}(2,:) = frac_active_perFrame_list_interp{pp(2)};
    combine_frac_active_interp{count}(3,:) = frac_active_perFrame_list_interp{pp(3)};
    mean_frac_active_interp{count} = mean(combine_frac_active_interp{count});
    stdev_frac_active_interp{count} = std(combine_frac_active_interp{count});
    up_error_frac_active_interp{count} =mean_frac_active_interp{count}+stdev_frac_active_interp{count};
    low_error_frac_active_interp{count} = mean_frac_active_interp{count}-stdev_frac_active_interp{count};

%     subplot(1, length(construct_names),count)
%     if count == 1; aa = 3; end
%     if count == 2; aa = 1; end
    h_line(count) = plot(x_Int,mean_frac_active_interp{count}, 'Color', string(color_codes(count)), 'LineWidth',1.5); hold on;
    patch([x_Int fliplr(x_Int)], smooth([low_error_frac_active_interp{count} fliplr(up_error_frac_active_interp{count})]), string(color_codes(count)), 'FaceAlpha',0.5, 'EdgeColor','none')
    %title(construct_names(count))
    xlabel("Progress through NC14")
    ylabel("Fraction active nuclei")
    ylim([0 1])
%     axis square
    alpha 0.3
    hold on;
    count = count + 1;
end
legend(h_line,construct_names);
hold off; 

%% generating heatmaps of MS2 signal over time (ACTIVE NUC ONLY)

% rangeIDX must be entered in ascending order at the very start of the code

% must run "interpolation" section before this

% sort data according to onset time
% plot using imagesc() and display "colorbar"

% interpolate x step and intensity (y) values to same length
x_rawData = cell(1,max(rangeIDX));
for i = rangeIDX%1:length(names)
    x_rawData{i} = linspace(0,1,size(M2_list{i},1));
end
x_Int = linspace(0,1,180);

M2_list_interp=cell(1,max(rangeIDX));
for i = rangeIDX
    for j = 1:size(M2_list{i},2)
    M2_list_interp{i}(:,j) = interp1(x_rawData{i},M2_list{i}(:,j), x_Int);
    end
end

count = 1;
midpoint = [];
reps_IDX = [1 2 3];
for kk = rangeIDX(1:3:end) % take every 3rd element
clear onFrame onset_frame_allNuc
% M2 = [M2_list_interp{kk}(:,active_list{reps_IDX(1)+(3*(count-1))}), M2_list_interp{kk+1}(:,active_list{reps_IDX(2)+(3*(count-1))}), M2_list_interp{kk+2}(:,active_list{reps_IDX(3)+(3*(count-1))})];
% use line below to plot individual embryo heat maps
M2 = [M2_list_interp{kk}(:,active_list{reps_IDX(1)+(3*(count-1))})];
onFrame = zeros(size(M2));
onFrame(M2>thr) = 1;
for i = 1:size(onFrame,2) %iterate through columns (i.e. nuclei)
    for j = 1:size(onFrame,1) %iterate through rows (i.e. frames)
        if onFrame(j,i) == 1
            onset_frame_allNuc(i) = j;
            break
        end
    end
end
[sortedOnset, sortedOrder] = sort(onset_frame_allNuc);
sortedM2 = M2(:,flip(sortedOrder));
rot270degSortedM2 = rot90(sortedM2,1);
x = linspace(0, 1, size(rot270degSortedM2,2));
y = [1: size(rot270degSortedM2,1)];
figure; imagesc(x,y,rot270degSortedM2,[400 4000]); colorbar; colormap turbo
ylabel('Nucleus number'); xlabel('Progress through NC14')
xlim([0 1])
title(construct_names(count)) 
set(gca,'YTick', [])

center = ceil(size(rot270degSortedM2,1)/2);
aboveThr = find(rot270degSortedM2(center,:)>thr); 
midpoint(count) = aboveThr(1);
numberActive(count) = size(rot270degSortedM2,1);

count = count + 1; 
end


%% UPDATED TO BOX CHART THROUGH HERE - all below is boxplot still

%% generate box plot of mRNA output for individual constructs - **input idx values**
% idx = [1:8]; % index of constructs to evaluate
% list = [(idx*3-2); (idx*3-2)+1; (idx*3-2)+2]';
% list = transpose(list);
% lin_list = reshape(list, [1, length(idx)*3]);
% 
% figure;
% mRNA_active_output_combined = reshape(mRNA_active_output_Matrix(:, lin_list),[],1);
% names_list_combined = reshape(names_list(:, lin_list),[],1);
% boxplot(mRNA_active_output_combined, names_list_combined)
% set(gca,'XTick',1:length(construct_names),'XTickLabel',construct_names(idx))
% ylabel('Mean mRNA production (AU)')
% 
% % combine mRNA output data by construct THEN remove top and bottom 10% each construct
% mRNA_active_output_Matrix_innerPercent = nan(size(names_list,1)*3, length(idx));
% % remove top and bottom 10% of data points
% count = 1;
% for i = idx
%     X = mRNA_active_output_combinedPerConstruct{i}; % X - list of mRNA output for data set i
%     percentile = prctile(X,[10 90]);
%     X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
%     mRNA_active_output_Matrix_innerPercent(1:length(X),count) = X;
%     count = count + 1;
% end
% 
% 
% % % remove 10% from each embryo data, before combining
% % mRNA_active_output_Matrix_innerPercent = nan(size(mRNA_active_output_Matrix));
% % % remove top and bottom 10% of data points
% % for i = 1:length(construct_names)
% %     X = mRNA_active_output_Matrix(:,i); % X - list of mRNA output for data set i
% %     percentile = prctile(X,[10 90]);
% %     X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
% %     mRNA_active_output_Matrix_innerPercent(:,i) = X;
% % end
% 
% % plot inner 80 percentile mRNA output
% figure;
% X_combined = reshape(mRNA_active_output_Matrix_innerPercent, [], 1);
% boxplot(X_combined, names_list_combined)
% set(gca,'XTick',1:length(idx),'XTickLabel',construct_names(idx))
% ylabel('Mean mRNA production (AU)')

%% generate box plot of activity onset time for each video independently
% figure;
% boxplot(onset_NC14prog_Matrix)
% set(gca,'XTick',1:size(onset_NC14prog_list,2),'XTickLabel',names)
% ylabel('Onset time (NC14 prog.)')
% title(normMethod, 'Interpreter', 'none')
 
%% Boxplot of onset times per construct
% % onset times - plotting all constructs combined (inner 80 percentile only)
% clear onset_NC14prog_perConstruct
% count = 1;
% toEval = [1:7]; % list of construct indices to evaluate
% for i = toEval
%     pp = [1:3]+((i-1)*3);
%     X = reshape(onset_NC14prog_Matrix(:,pp),[],1);
%     percentile = prctile(X,[10 90]); % can comment out this line and next to keep all data on plots
%     X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
%     onset_NC14prog_perConstruct(1:3*(size(onset_NC14prog_Matrix,1)),count) = X;
%     count = count + 1;
% end
% 
% figure; boxplot(onset_NC14prog_perConstruct,'Labels', construct_names(toEval));
% ylabel("Onset time (frac through NC14)")
% ylim([0 1])
% title(normMethod, 'Interpreter', 'none')

%% generate box plots of mean amplitudes for active nuclei (only in their active frames)

% boxplot of mean amplitudes per embryo
figure;
boxplot(mean_amplitude_active_Matrix)
set(gca,'XTick',1:size(mean_amplitude_active_list,2),'XTickLabel',names)
ylabel('Mean traj amplitude (active frames only)')

% MEAN amplitudes - plotting all constructs combined (inner 80 percentile only)
clear mean_amplitude_active_perConstruct
count = 1;
toEval = [1:7]; % list of construct indices to evaluate
for i = toEval
    pp = [1:3]+((i-1)*3);
    X = reshape(mean_amplitude_active_Matrix(:,pp),[],1);
    percentile = prctile(X,[10 90]); % can comment out this line and next to keep all data on plots
    X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
    mean_amplitude_active_perConstruct(1:3*(size(mean_amplitude_active_Matrix,1)),count) = X;
    count = count + 1;
end

figure; boxplot(mean_amplitude_active_perConstruct,'Labels', construct_names(toEval));
ylabel("Mean transcriptional burst amplitude (AU)")

%% generate box plots of max amplitudes for active nuclei (only in their active frames)

% boxplot of max amplitudes per embryo 
figure;
boxplot(max_amplitude_active_Matrix)
set(gca,'XTick',1:size(max_amplitude_active_list,2),'XTickLabel',names)
ylabel('Max traj amplitude (active frames only)')

% MAX amplitudes - plotting all constructs combined (inner 80 percentile only)
clear max_amplitude_active_perConstruct
count = 1;
toEval = [1:7]; % list of construct indices to evaluate
for i = toEval
    pp = [1:3]+((i-1)*3);
    X = reshape(max_amplitude_active_Matrix(:,pp),[],1);
    percentile = prctile(X,[10 90]); % can comment out this line and next to keep all data on plots
    X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
    max_amplitude_active_perConstruct(1:3*(size(max_amplitude_active_Matrix,1)),count) = X;
    count = count + 1;
end
figure; boxplot(max_amplitude_active_perConstruct,'Labels', construct_names(toEval));
ylabel("Max transcriptional burst amplitude (AU)")




%% Plotting as scatter plot on representative x-axis
% manual input for desired constructs to plot
x = [0 1.3 6.3 7.6 7.7 8.8 10];
%y_list = flip([2:7]);
y_list = 1:9;
clear y
clear err
% y = nan(1,length(y_list));
% err = nan(1,length(y_list));

count = 1; 
for i = y_list
    y(count) = median(mRNA_active_output_combinedPerConstruct{i});
    err(count) = sqrt(std(mRNA_active_output_combinedPerConstruct{i}));0
    count = count + 1;
end

enh_up = [1 2 3 4 6];
enh_down = [5 7];

figure;

errorbar(x(enh_up), y(enh_up), err(enh_up), 'Color','b','LineWidth',2, 'MarkerSize',1);
ylabel('mRNA output (A.U.)')
hold on;
errorbar(x(enh_down), y(enh_down), err(enh_down), 'Color','r','LineWidth',2, 'MarkerSize',1);
xlabel('E-P distance (kb)')
text(100, 800,"error = sqrt(stdev)")
ylim([0 1800])
legend({'enh-UP', 'enh-DONW'})
hold off;

%% Total expression boundary width 

% plotting max boundary width expression
max_width = nan(1,length(names));
max_width_combined = cell(1,length(names)/3);
for i = 1:length(names) 
    max_width(i) = (max(cy_active_list{i})-min(cy_active_list{i}))/nuc_diameter_est;
    pp = ceil(i/3);
    max_width_combined{pp}(end+1) = max_width(i);
end

for i = 1: length(names)/3
mean_max_width_combined(i) = mean(max_width_combined{i});
stdError_max_width_combined(i) = std(max_width_combined{i}) / sqrt(length(max_width_combined{i}));
end

figure;
bar([1:length(names)/3], mean_max_width_combined); 
hold on;
errlow = stdError_max_width_combined; 
errhigh = stdError_max_width_combined;
er = errorbar([1:length(names)/3], mean_max_width_combined, errlow, errhigh);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylim([0 30])
xticklabels(construct_names)
ylabel("DV expression domain width (# nuclei)")


%% Total Number of Active Nuclei

% plotting max boundary width expression
active_count = nan(1,length(names));
active_count_combined = cell(1,length(names)/3);
for i = 1:length(names) 
    active_count(i) = length(active_rotated_list{i});
    pp = ceil(i/3);
    active_count_combined{pp}(end+1) = active_count(i);
end

for i = 1: length(names)/3
mean_active_count_combined(i) = mean(active_count_combined{i});
stdError_active_count_combined(i) = std(active_count_combined{i}) / sqrt(length(active_count_combined{i}));
end

figure;
bar([1:length(names)/3], mean_active_count_combined); 
hold on;
errlow = stdError_active_count_combined; 
errhigh = stdError_active_count_combined;
er = errorbar([1:length(names)/3], mean_active_count_combined, errlow, errhigh);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
%ylim([0 30])
xticklabels(construct_names)
ylabel("# active nuclei")


%% UPDATED UP THROUGH HERE - 11/28/23

%%

% %% cumulative active nuclei across single constructs
figure;
count = 1;
for i = rangeIDX%1:length(names)
    pp = ceil(i/3);
    subplot(2, round(length(construct_names)/2),count); hold on;
    y = interp1(x_rawData{i},cumulative_active_count_list{i}, x_Int);
    cumulative_active_count_list_interp{i} = y;
    plot(x_Int,y);
    %plot(linspace(0,1,length(cumulative_active_count_list{i})),cumulative_active_count_list{i})%./max(cumulative_active_count_list{i}))
    title(append("Cumulative active nuc - ", construct_names{count}))
    legend('emb01', 'emb02', 'emb03', 'Location','southeast', 'FontSize', 8)
    ylim([0 400])
    axis square
    count = count + 1;
end

% combine triplicate data sets using shaded error bars

% first we combine the data for each construct into cells
% figure; sgtitle("Cumulative active nuclei")
% for i = 1:length(construct_names) % # of constructs
%     pp = (i-1)*3+1:(i-1)*3+3;
%     combine_cumulative_active_interp{i}(1,:) = cumulative_active_count_list_interp{pp(1)};
%     combine_cumulative_active_interp{i}(2,:) = cumulative_active_count_list_interp{pp(2)};
%     combine_cumulative_active_interp{i}(3,:) = cumulative_active_count_list_interp{pp(3)};
%     mean_cumulative_active_interp{i} = mean(combine_cumulative_active_interp{i});
%     stdev_cumulative_active_interp{i} = std(combine_cumulative_active_interp{i});
%     up_error_cumulative_active_interp{i} =mean_cumulative_active_interp{i}+stdev_cumulative_active_interp{i};
%     low_error_cumulative_active_interp{i} = mean_cumulative_active_interp{i}-stdev_cumulative_active_interp{i};
% 
%     subplot(1, length(construct_names),i)
%     %     if i == 2 || i == 4 || i == 5
%     plot(x_Int,mean_cumulative_active_interp{i}, 'Color', string(color_codes(i)), 'LineWidth',1.5); hold on;
%     patch([x_Int fliplr(x_Int)], smooth([low_error_cumulative_active_interp{i} fliplr(up_error_cumulative_active_interp{i})]), string(color_codes(i)), 'FaceAlpha',0.5, 'EdgeColor','none')
%     %title(construct_names(i))
%     xlabel("Progress through NC14")
%     ylabel("Cumulative active nuclei")
%     legend('10kb','','8.8kb', '', '7.7kb', '' , '6.3kb', '', '1.3kb', 'Location', 'northwest')
%     axis square
%     ylim([0 400])
%     alpha 0.3
% %     end
% end


%% fraction active nuclei per frame across single constructs

%individual plots
figure;
count = 1; 
for i = rangeIDX% 1:length(names)
    pp = ceil(count/3);
    subplot(1, length(construct_names),pp); hold  on;
    y = interp1(x_rawData{i},frac_active_perFrame_list{i}, x_Int);
    frac_active_perFrame_list_interp{i} = y;
    plot(x_Int,y);
%     plot(linspace (0,1,length(frac_active_perFrame_list{i})),frac_active_perFrame_list{i})
    title(append( construct_names{pp}))
    legend('emb01', 'emb02', 'emb03', 'Location','northwest')
    ylim([0 1])
    xlabel('Frac Progress through NC14')
    ylabel('Instantaneous frac active nuclei')
    axis square
    count = count + 1; 
end

% shaded error bars

% first we combine the data for each construct into cells
figure; sgtitle("Fraction active nuclei per frame")
count = 1;
constructIDX = rangeIDX(1:3:end); % take every 3rd element
for i = constructIDX%1:length(construct_names) % # of constructs
%     i = ceil(count/3);
%     pp = (i-1)*3+1:(i-1)*3+3;
    pp = [i i+1 i+2];
    combine_frac_active_interp{count}(1,:) = frac_active_perFrame_list_interp{pp(1)};
    combine_frac_active_interp{count}(2,:) = frac_active_perFrame_list_interp{pp(2)};
    combine_frac_active_interp{count}(3,:) = frac_active_perFrame_list_interp{pp(3)};
    mean_frac_active_interp{count} = mean(combine_frac_active_interp{count});
    stdev_frac_active_interp{count} = std(combine_frac_active_interp{count});
    up_error_frac_active_interp{count} =mean_frac_active_interp{count}+stdev_frac_active_interp{count};
    low_error_frac_active_interp{count} = mean_frac_active_interp{count}-stdev_frac_active_interp{count};

    subplot(1, length(construct_names),count)
    plot(x_Int,mean_frac_active_interp{count}, 'Color', string(color_codes(count)), 'LineWidth',1.5); hold on;
    patch([x_Int fliplr(x_Int)], smooth([low_error_frac_active_interp{count} fliplr(up_error_frac_active_interp{count})]), string(color_codes(count)), 'FaceAlpha',0.5, 'EdgeColor','none')
    title(construct_names(count))
    xlabel("Progress through NC14")
    ylabel("Fraction active nuclei")
    ylim([0 1])
    axis square
    alpha 0.3
    count = count + 1;
end
%% Compare total number of active nuclei for each video independently

% Scatter plot

% figure; scatter(1:length(total_active_count_list), cell2mat(total_active_count_list),'*')
% ylim([0 max(cell2mat(total_active_count_list))+10])
% set(gca,'xtick',[1:length(total_active_count_list)],'xticklabel',names)

% Bar chart
figure; bar(1:length(names), cell2mat(total_active_count_list))
ylim([0 max(cell2mat(total_active_count_list))+10])
set(gca,'xtick',[1:length(total_active_count_list)],'xticklabel',names)
ylabel('Total count of active nuclei')

%% Compare total number of active nuclei per construct - update code if we add more constructs to analyze

% generate table of values
array_total_active_count_list = cell2mat(total_active_count_list);
% if length(array_total_active_count_list) == 14
%     array_total_active_count_list = [array_total_active_count_list NaN];
% end
matrix_total_active_count_list = reshape(array_total_active_count_list,[3,construct_count]);


% generate boxplot 
figure;
boxplot(matrix_total_active_count_list)
set(gca,'XTick',1:construct_count,'XTickLabel',construct_names)
ylabel('Total count of active nuclei')


%% Compare time spent in active state per video
% boxplot
figure;
boxplot(frac_time_active_Matrix)
set(gca,'XTick',1:size(frac_time_active_Matrix,2),'XTickLabel',names)
ylabel('Fraction of NC14 spent in ACTIVE state')


%% Compare time spent in active state per construct
clear frac_time_active_perConstruct
count = 1;
toEval = [1:7]; % list of construct indices to evaluate
for i = toEval
    pp = [1:3]+((i-1)*3);
    X = reshape(frac_time_active_Matrix(:,pp),[],1);
    percentile = prctile(X,[10 90]); % can comment out this line and next to keep all data on plots
    X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
    frac_time_active_perConstruct(1:3*(size(frac_time_active_Matrix,1)),count) = X;
    count = count + 1;
end

figure; boxplot(frac_time_active_perConstruct,'Labels', construct_names(toEval));
ylabel("Frac NC14 spent in ACTIVE state")

%% Compare time spent in active state (SINCE ONSET) per video
% boxplot
figure;
boxplot(frac_time_active_sinceOnset_Matrix)
set(gca,'XTick',1:size(frac_time_active_sinceOnset_Matrix,2),'XTickLabel',names)
ylabel('Fraction of NC14 spent in ACTIVE state since onset')

%% Compare time spent in active state (SINCE ONSET) per construct
clear frac_time_active_sinceOnset_perConstruct
count = 1;
toEval = [1:7]; % list of construct indices to evaluate
for i = toEval
    pp = [1:3]+((i-1)*3);
    X = reshape(frac_time_active_sinceOnset_Matrix(:,pp),[],1);
    percentile = prctile(X,[10 90]); % can comment out this line and next to keep all data on plots
    X(X<percentile(1) | X>percentile(2)) = NaN; % insert NaN for any numbers outside the inner 80 percentile -- NaN excluded from box plots!
    frac_time_active_sinceOnset_perConstruct(1:3*(size(frac_time_active_sinceOnset_Matrix,1)),count) = X;
    count = count + 1;
end

figure; boxplot(frac_time_active_sinceOnset_perConstruct,'Labels', construct_names(toEval));
ylabel("Frac NC14 spent in ACTIVE state SINCE ONSET")
ylim([0 1.1])

%% Boundary width over time

% figure;
for i = 1:length(names)
    pp = ceil(i/3);
    subplot(2, (length(construct_names)/2),pp); hold on;
    plot(linspace(0,1,length(boundary_width_list{i})),boundary_width_list{i}./nuc_diameter_est)
    title(append('Boundary width - ', construct_names{pp}))
    ylabel('width (# nuclei)')
    xlabel('Progress through NC14')
    legend('emb01', 'emb02', 'emb03', 'Location','northwest')
    ylim([0 160])
    axis square
end



%%

% %% Cumulative mRNA output over time
% 
% figure;
% for i = 1:length(names)
%     pp = ceil(i/3);
%     subplot(1, length(construct_names),pp); hold on;
%     y = interp1(x_rawData{i},cumulative_outputM_perFrame_list{i}, x_Int);
%     cumulative_outputM_perFrame_list_interp{i} = y;
%     plot(x_Int,y);
% %     plot(linspace(0,1,length(cumulative_outputM_perFrame_list{i})),cumulative_outputM_perFrame_list{i})%./max(cumulative_outputM_perFrame_list{i})) % cumulative production normalized the total mRNA produced (i.e. fraction of total produced at given point)
%     title(append('Cumulative mRNA output over time - ', construct_names{pp}))
%     ylabel("Total mRNA production")
%     xlabel("Progress through NC14")
%     legend('emb01', 'emb02', 'emb03', 'Location','northwest')
%     ylim([0 max(nansum(mRNA_active_output_Matrix))])
% %     ylim([0 1])
%     axis square
% 
% end
% 
% % combine triplicate data into shaded error bars plot
% 
% % first we combine the data for each construct into cells
% figure; sgtitle("mRNA output")
% for i = 1:length(construct_names) % # of constructs
%     pp = (i-1)*3+1:(i-1)*3+3;
%     combine_cumulative_outputM_perFrame_interp{i}(1,:) = cumulative_outputM_perFrame_list_interp{pp(1)};
%     combine_cumulative_outputM_perFrame_interp{i}(2,:) = cumulative_outputM_perFrame_list_interp{pp(2)};
%     combine_cumulative_outputM_perFrame_interp{i}(3,:) = cumulative_outputM_perFrame_list_interp{pp(3)};
% 
%     mean_cumulative_outputM_interp{i} = nanmean(combine_cumulative_outputM_perFrame_interp{i});
%     stdev_cumulative_outputM_interp{i} = nanstd(combine_cumulative_outputM_perFrame_interp{i});
%     up_error_cumulative_outputM_interp{i} =mean_cumulative_outputM_interp{i}+stdev_cumulative_outputM_interp{i};
%     low_error_cumulative_outputM_interp{i} = mean_cumulative_outputM_interp{i}-stdev_cumulative_outputM_interp{i};
%     
%     if i == 2 || i == 4 || i == 5
% %     subplot(1, length(construct_names),i)
%     plot(x_Int,mean_cumulative_outputM_interp{i}, color_codes(i), 'LineWidth',1.5); hold on;
%     
%     patch([x_Int fliplr(x_Int)], smooth([low_error_cumulative_outputM_interp{i} fliplr(up_error_cumulative_outputM_interp{i})]), color_codes(i), 'FaceAlpha',0.5, 'EdgeColor','none')
%     %title(construct_names(i))
%     xlabel("Progress through NC14")
%     ylabel("Total mRNA production (AU)")
%     %     legend('10kb', '', '8.8kb', '' ,'7.7kb', '', '6.3kb', '', '1.3kb')
%     legend('8.8kb', '' , '6.3kb', '', '1.3kb')
%     ylim([0 max(nansum(mRNA_active_output_Matrix))])
%     xlim([0 1])
%     axis square
%     alpha 0.2
%     end
% end

% %% Cumulative mRNA output over time - binned into 10% progress through NC14
% 
% binning_intervals = linspace(0.1,1,10);
% binned_cumulative_mRNA = zeros(length(binning_intervals),length(names));
% for i = 1:length(names)
%     x = round(binning_intervals.*length(cumulative_outputM_perFrame_list{i}));
%     binned_cumulative_mRNA(1,i) = cumulative_outputM_perFrame_list{i}(x(1));
%     binned_cumulative_mRNA(2:length(x),i) = cumulative_outputM_perFrame_list{i}(x(2:length(x)))-cumulative_outputM_perFrame_list{i}(x(1:length(x)-1));
% %     binned_cumulative_mRNA(:,i) = binned_cumulative_mRNA(:,i)./max(cumulative_outputM_perFrame_list{i}); %normalize to fraction of total mRNA
% end
% 
% figure; sgtitle("Binned mRNA Output")
% for i = 1:length(construct_names)
%     pp = [1:3]+((i-1)*3);
%     subplot(1, length(construct_names),i); hold on;
%     bar(binned_cumulative_mRNA(:,pp));
%     set(gca,'xtick',[1:length(binning_intervals)],'xticklabel',append(string((binning_intervals-0.1).*100),'-',string(binning_intervals.*100)))
%     %title(append('Binned mRNA output - ', construct_names{i}))
%     title(construct_names{i})
%     ylabel("amount of mRNA produced")
%     xlabel("Binned % progress through NC14")
%     legend('emb01', 'emb02', 'emb03', 'Location','northwest')
%     ylim([0 max(max(binned_cumulative_mRNA))])
%     % ylim([0 0.25])
%     axis square
% end
% 
% %% Instances of active-to-inactive states (and vice versa)
% bursting_on = nan(max_nuc_active, length(names));
% bursting_off = nan(max_nuc_active, length(names));
% 
% for i = 1:length(names)
%     bursting_on(1:length(inactive_to_active_count_perNuc_list{i}),i) = inactive_to_active_count_perNuc_list{i};
%     bursting_off(1:length(active_to_inactive_count_perNuc_list{i}),i) = active_to_inactive_count_perNuc_list{i};
% end
% 
% figure;
% boxplot(bursting_on)
% set(gca,'XTick',1:size(mRNA_active_output_Matrix,2),'XTickLabel',names)
% ylabel('Instances of inactive to active - "bursting on"')
% ylim([0 max(max(bursting_on))+1])
% 
% figure;
% boxplot(bursting_off)
% set(gca,'XTick',1:size(mRNA_active_output_Matrix,2),'XTickLabel',names)
% ylabel('Instances of active to inactive')
% ylim([0 max(max(bursting_off))+1])

%% Plotting active burst durations
burst_duration_matrix = nan(10000, length(names));
for i = 1:length(names)
    burst_duration_matrix(1:length(burst_duration_list{i}),i) = burst_duration_list{i};
end

figure;
boxplot(burst_duration_matrix)
set(gca,'XTick',1:size(mRNA_active_output_Matrix,2),'XTickLabel',names)
ylabel('Duration of bursts (# of frames)')
% ylim([0 max(max(bursting_on))+1])


%% mRNA expression along dorsal-ventral axis

figure; sgtitle("mRNA output across DV axis")
count = 1;
plotIdx = [1:10];
for i = plotIdx
    pp = [1:3]+((i-1)*3);
    mRNAout_combined = [mRNA_active_output_Matrix(1:active_output_lengths(pp(1)),pp(1));...
        mRNA_active_output_Matrix(1:active_output_lengths(pp(2)),pp(2));...
        mRNA_active_output_Matrix(1:active_output_lengths(pp(3)),pp(3))]; 
    
    cy_active_combined = [cy_active_normalized{pp(1)}, cy_active_normalized{pp(2)}, cy_active_normalized{pp(3)}]';
    subplot(1,length(plotIdx),count) 
    scatter(cy_active_combined, mRNAout_combined, 'k.'); hold on; xline(0);
    title(construct_names(i))
    xlabel('DV position (pixels)')
    ylabel('mRNA output of active nuclei')
    xlim([-256 256])
    ylim([0 5000])
    axis square
    hold off;
    count = count + 1;
end

% %% mRNA output normalized to fraction of active frames
% % for each nucleus, divide mRNA output by fraction of movie since onset time
% 
% frac_since_Onset = nan(size(mRNA_active_output_Matrix));
% for i = 1:length(names)
%     onset = onset_NC14prog_Matrix(:,i);
%     onset = onset(~isnan(onset));
%     frac_since_Onset(1:length(onset),i) = 1 - onset; % fraction of movie since onset of activation 
% end
% 
% mRNA_output_norm_to_activeFrames = mRNA_active_output_Matrix./frac_since_Onset;
% 
% figure;
% boxplot(mRNA_output_norm_to_activeFrames,'symbol','')
% set(gca,'XTick',1:size(mRNA_output_norm_to_activeFrames,2),'XTickLabel',names)
% ylabel('Total mRNA production (AU)')
% title('mRNA output normalized to time after onset')
% ylim([0 8000])