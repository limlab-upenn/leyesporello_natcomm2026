% RNA FISH image analysis

% comparing 7.6kb enh-UP vs enh-DOWN

clear global

% select which construct to analyze; 1 = enhUP; 2 = enhDOWN
% construct_select = 2; 

for construct_select = 1:2
if construct_select == 1
    embryos_to_analyze =  ["02" "03" "04" "05" "06" "07" "08" "09" "11" "12" "13" "14" "15"]; % missing emb01 and emb10
    pathname = "C:\Users\Emilia Leyes Porello\OneDrive - PennO365\Penn Research\Enhancer-Promoter Spacers\Fly imaging\FISH images\7.6kb-enhUP-t.yel\TIF files\";
    construct_nickName = "7.6kb-enhUP-t.Yel"; 
elseif construct_select == 2
    embryos_to_analyze =  ["01" "02" "03" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15"]; % missing emb04
    pathname = "C:\Users\Emilia Leyes Porello\OneDrive - PennO365\Penn Research\Enhancer-Promoter Spacers\Fly imaging\FISH images\7.7kb-enhDOWN-t.yel\TIF files\";
    construct_nickName = "7.7kb-enhDOWN-t.Yel"; 
end
disp(append("Evaluating ", construct_nickName))

sensitivity_metrics = cell(1,length(embryos_to_analyze));

A = cell(1, length(embryos_to_analyze));
B = cell(1, length(embryos_to_analyze));
C = cell(1, length(embryos_to_analyze));
D = cell(1, length(embryos_to_analyze));

for j = 1:length(embryos_to_analyze)

% pathname = "C:\Users\Emilia Leyes Porello\OneDrive - PennO365\Penn Research\Enhancer-Promoter Spacers\Fly imaging\FISH images\7.6kb-enhUP-t.yel\TIF files\";
filename = append("MAX_emb", embryos_to_analyze(j), ".tif");

fname = append(pathname, filename);
% info = imfinfo(fname);
% num_images = numel(info);

export = [pathname 'segmentation\'];
export2 = pathname;

% Define each channel
his = 1; % Histone channel - DAPI - BLUE
snail_1 = 2; % sna-BIOT - 555nm - RED
yellow = 3; % yellow(yw)-FITC - 488nm - GREEN
snail_2 = 4; % sna-DIG - 657nm - MAGENTA

% Read in each channel
A{j} = imread(fname, his); % DAPI - blue
B{j} = imread(fname, snail_1); % snail - red
C{j} = imread(fname, yellow); % yellow - green
D{j} = imread(fname, snail_2); % snail - magenta

%% Calculating yellow probe singal metrics (can also do sensitivity analysis for selection of intensity threshold)

% Evaluate channel C (yellow probes) only

thr = 3000; % 1000:100:5000; % vary threshold or define a single one

C_filtered = cell(1,length(thr));
sensitivity_metrics{j} = nan(9,4); % col 1 - threshold value; col 2 - nonzero pixel count; col 3 - mean intensity of signals; col 4 - stdev signal intensity

for i = 1:length(thr)
    C_filtered{i} = C{j};
    C_filtered{i}(C_filtered{i}<thr(i))=0;
    nonZeroValues = C_filtered{i}; 
    nonZeroValues(nonZeroValues==0) = [];
    sensitivity_metrics{j}(i,1) = thr(i); % threshold value
    sensitivity_metrics{j}(i,2) = nnz(C_filtered{i}); % nonzero pixel count
    sensitivity_metrics{j}(i,3) = mean(double(nonZeroValues)); % mean of nonzero pixels
    sensitivity_metrics{j}(i,4) = std(double(nonZeroValues)); % median of nonzero pixels
    
end
end

% Ploting sensitivity metrics

% Plot nonzero pixel counts as function of threshold
figure;
subplot(1,3,1)
for i = 1:length(embryos_to_analyze)
    plot(sensitivity_metrics{i}(:,1), sensitivity_metrics{i}(:,2), '.-');
    hold on
end
xlabel('signal threshold value')
ylabel('nonzero pixel count')
%legend(embryos_to_analyze)
    
% Plot mean signal intensity as function of threshold
subplot(1,3,2)
for i = 1:length(embryos_to_analyze)
    plot(sensitivity_metrics{i}(:,1), sensitivity_metrics{i}(:,3), '.-');
    hold on
end
xlabel('signal threshold value')
ylabel('mean signal intensity')
%legend(embryos_to_analyze, 'Location', 'northwest')


% Plot stddev of signal intensities as function of threshold
subplot(1,3,3)
for i = 1:length(embryos_to_analyze)
    plot(sensitivity_metrics{i}(:,1), sensitivity_metrics{i}(:,4), '.-');
    hold on
end
xlabel('signal threshold value')
ylabel('stddev signal intensity')
legend(embryos_to_analyze)

sgtitle(construct_nickName)

% plot visual of emb 15 to qualitatively compare thresholds
% figure; 
% subplot(1,4,1); imshow(C{j}*10); title('original'); 
% C_filt = C{j}; 
% C_filt(C_filt<2000) = 0;
% subplot(1,4,2); imshow(C_filt*10); title('thr = 2000'); 
% C_filt(C_filt<3000) = 0;
% subplot(1,4,3); imshow(C_filt*10); title('thr = 3000');
% C_filt(C_filt<4000) = 0;
% subplot(1,4,4); imshow(C_filt*10); title('thr = 4000'); 
% 
% sgtitle(append(construct_nickName,' - emb', string(embryos_to_analyze(j))))

%% Identifying snail expression domain 

rpo_snailDomain = cell(1, length(embryos_to_analyze));

for j = 1:length(embryos_to_analyze)
    P1 = imgaussfilt(D{j},30); % gaussian filter/blur
    P2 = adapthisteq(P1, "NumTiles",[12 12]);
    bw = imbinarize(P2, graythresh(P2));
    bw2 = bwareaopen(bw, 50000);

    % to fill holes that contact the edges, we pad the relevant corners, one at a
    % time, call imfill, then remove the padding before repeating on the
    % other corners

    % to identify the relevant corners, we first need to identify the
    % orientation of the embryo (which diagonal is it on?) - we will use an
    % identity matrix and it's flip to do so

    identity = eye(length(D{j})); % create 1024 x 1024 identity matrix
    identify_flip = flip(identity); % flip the diagonal

%     test_1 = sum(sum(bw2.*identity)); % measure of diagonal alignment top left to bottom right (back slash)
%     test_2 = sum(sum(bw2.*identify_flip)); % measure of diagonal alignment bottom left to top right (forward slash)

    test_1 = sum(sum(bw2(1:200, 1:200)));
    test_2 = sum(sum(bw2(1:200, length(D{j})-200:length(D{j}))));
    
    if test_1 > test_2 % embryo in back slash orientation - pad top left and bottom right corners

        bw2(1,:) = 1; % pad top edge
        bw2(:,1) = 1; % pad left edge
        bw_a = imfill(bw2, 'holes');
        bw2(1,:) = 0; bw_a(1,:) = 0; % clear top edge
        bw2(:,1) = 0; bw_a(:,1) = 0; % clear left edge
 
        bw2(length(D{j}),:) = 1; % pad bottom edge
        bw2(:,length(D{j})) = 1; % pad right edge
        bw_b = imfill(bw2, 'holes'); 
        bw2(length(D{j}),:) = 0; bw_b(length(D{j}),:) = 0; % clear bottom edge
        bw2(:,length(D{j})) = 0; bw_b(:,length(D{j})) = 0; % clear right edge

        bw3 = bw_a | bw_b; 

    elseif test_2 > test_1 % embryo in forward slash orientation - pad bottom left and top right corners
               
        bw2(length(D{j}),:) = 1; % pad bottom edge
        bw2(:,1) = 1; % pad left edge
        bw_a = imfill(bw2, 'holes');
        bw2(length(D{j}),:) = 0; bw_a(length(D{j}),:) = 0; % clear bottom edge
        bw2(:,1) = 0; bw_a(:,1) = 0; % clear left edge
 
        bw2(1,:) = 1; % pad top edge
        bw2(:,length(D{j})) = 1; % pad right edge
        bw_b = imfill(bw2, 'holes'); 
        bw2(1,:) = 0; bw_b(1,:) = 0; % clear top edge
        bw2(:,length(D{j})) = 0; bw_b(:,length(D{j})) = 0; % clear right edge

        bw3 = bw_a | bw_b; 

    end

    rpo_snailDomain{j} = regionprops(bw3, 'area', 'PixelList'); % area of objects (select largest for domain area)

%     figure; subplot(1,5,1); imshow(D{j}*7); subplot(1,5,2); imshow(P1*7); subplot(1,5,3); imshow(bw); subplot(1,5,4); imshow(bw2); subplot(1,5,5); imshow(bw3); sgtitle(j)
end


%% Calculating signals normalized to snail domain area

C_filt_3000 = cell(1, length(embryos_to_analyze));
rpo_yellowSignals = cell(1, length(embryos_to_analyze));

for j = 1: length(embryos_to_analyze)
    C_filt_3000{j} = C{j};
    C_filt_3000{j}(C_filt_3000{j}<3000)=0;
    C_filt_3000_bw = imbinarize(C_filt_3000{j});
    rpo_yellowSignals{j} = regionprops(C_filt_3000_bw, 'area', 'PixelList');
end

% distribution of signal sizes (rpo area for each signal)
signal_size_list = [];
for i = 1:length(embryos_to_analyze) 
    signal_size_list = [signal_size_list, rpo_yellowSignals{i}.Area]; %#ok<AGROW> 
end

% normalizing signal count (ID'd by rpo) to snail expression domain
sig_count = [];
sig_count_normArea = [];
for i = 1:length(embryos_to_analyze) 
    sig_count(i) = length(rpo_yellowSignals{i}); % number of signals identified by rpo in emb i
    sig_count_normArea(i) = sig_count(i)/max(rpo_snailDomain{i}.Area); % number of signals / snail domain area calc'd by rpo
end

disp('start')
% calculating average intensity of each signal
sig_intensity_all = nan(max(sig_count),length(embryos_to_analyze));
sig_intensity_2pix = nan(max(sig_count),length(embryos_to_analyze));
sig_intensity_5pix = nan(max(sig_count),length(embryos_to_analyze));
sig_intensity_10pix = nan(max(sig_count),length(embryos_to_analyze));
for i = 1:length(embryos_to_analyze) % each emb
    for j = 1:length(rpo_yellowSignals{i}) % each signal in emb i
        pix = rpo_yellowSignals{i}(j,:).PixelList;
        for k = 1:length(pix) % each pixel in signal j
            x = pix(:,1);
            y = pix(:,2);
            pix_intensities = impixel(C{i},x,y);
            sig_intensity_all(j,i) = mean(pix_intensities(:,1)); % mean intensity of pixels identified as part of the signal
            if length(pix) > 2
                sig_intensity_2pix(j,i) = mean(pix_intensities(:,1)); % mean intensity of pixels identified as part of the signal
            elseif length(pix) > 5
                sig_intensity_5pix(j,i) = mean(pix_intensities(:,1)); % mean intensity of pixels identified as part of the signal
            elseif length(pix) > 10
                sig_intensity_10pix(j,i) = mean(pix_intensities(:,1)); % mean intensity of pixels identified as part of the signal
            end
        end
    end    
end
disp('end')

% calculating average intensity of yellow probe within snail expression domain
domain_intensity_mean = nan(1,length(embryos_to_analyze));
domain_intensity_distribution = [];
for i = 1:length(embryos_to_analyze) % each emb
    pix = rpo_snailDomain{i}(:).PixelList;
    x = pix(:,1);
    y = pix(:,2);
    pix_intensities = impixel(C{i},x,y);
    domain_intensity_mean(i) = mean(pix_intensities(:,1));
    domain_intensity_distribution = [domain_intensity_distribution pix_intensities(:,1)'];
end

% end of construct_select loop -- save or rename variables we want to compare across constructs here
construct_names{construct_select} = construct_nickName;
sig_count_byConstruct{construct_select} = sig_count;
sig_count_normArea_byConstruct{construct_select} = sig_count_normArea;
sig_intensity_all_byConstruct{construct_select} = sig_intensity_all;
sig_intensity_2pix_byConstruct{construct_select} = sig_intensity_2pix;
sig_intensity_5pix_byConstruct{construct_select} = sig_intensity_5pix;
sig_intensity_10pix_byConstruct{construct_select} = sig_intensity_10pix;
signal_size_list_byConstruct{construct_select} = signal_size_list;
domain_intensity_mean_byConstruct{construct_select} = domain_intensity_mean;
domain_intensity_distribution_byConstruct{construct_select} = domain_intensity_distribution; 
rpo_snailDomain_byConstruct{construct_select} = rpo_snailDomain; 

clearvars -except  sig_count_byConstruct sig_count_normArea_byConstruct...
    sig_intensity_all_byConstruct sig_intensity_2pix_byConstruct construct_names ...
    sig_intensity_5pix_byConstruct sig_intensity_10pix_byConstruct ...
    signal_size_list_byConstruct domain_intensity_mean_byConstruct ...
    domain_intensity_distribution_byConstruct rpo_snailDomain_byConstruct sig_intensity_5pix_byConstruct sig_intensity_10pix_byConstruct

end

%% Once all data is generated - plot/display it

% plot boxchart of signal count
plot_colors = ['b', 'r'];
figure;
for i = 1:length(construct_names)
    xGroup = ones(1,length(sig_count_byConstruct{i})).*i; % create xGroup list for boxchart input
    boxchart(xGroup, sig_count_byConstruct{i}, "BoxFaceColor", plot_colors(i));
    hold on;
end
xticks(1:length(construct_names))
xticklabels(construct_names)
ylabel('Number of Signals')
title('RNA FISH - yellow probes')

% plot boxchart of signal count normalized by area (muliplied by scalar 1000 since normalization by area yielded very tiny numbers)
plot_colors = ['b', 'r'];
figure;
for i = 1:length(construct_names)
    xGroup = ones(1,length(sig_count_normArea_byConstruct{i})).*i; % create xGroup list for boxchart input
    boxchart(xGroup, sig_count_normArea_byConstruct{i}.* 1000, 'MarkerColor','k', 'BoxFaceColor', plot_colors(i));
    hold on;
end
xticks(1:length(construct_names))
xticklabels(construct_names)
ylabel('Number of Signals (norm to snail domain area)')
title('RNA FISH - yellow probes')

% plot boxchart of average signal intensities
plot_colors = ['b', 'r'];
figure;
for i = 1:length(construct_names)
    sig_intensity_byConstruct_linear{i} = reshape(sig_intensity_2pix_byConstruct{i},1,[]); % only signals greater than 2 pixels
%     sig_intensity_byConstruct_linear{i} = reshape(sig_intensity_all_byConstruct{i},1,[]); % all sized signals
    xGroup = ones(1,length(sig_intensity_byConstruct_linear{i})).*i; % create xGroup list for boxchart input
    boxchart(xGroup, sig_intensity_byConstruct_linear{i}, "BoxFaceColor", plot_colors(i));
    hold on;
end
xticks(1:length(construct_names))
xticklabels(construct_names)
ylabel('signal intensity (A.U.)')
title('RNA FISH - yellow probes')


% plot boxchart of signal areas/size %%%%
plot_colors = ['b', 'r'];
figure;
for i = 1:length(construct_names)
    signal_size_list_byConstruct_linear{i} = reshape(signal_size_list_byConstruct{i},1,[]);
    xGroup = ones(1,length(signal_size_list_byConstruct_linear{i})).*i; % create xGroup list for boxchart input
    boxchart(xGroup, signal_size_list_byConstruct_linear{i}, "BoxFaceColor", plot_colors(i));
    hold on;
end
xticks(1:length(construct_names))
xticklabels(construct_names)
ylabel('signal area (pixel count)')
title('RNA FISH - yellow probes')


% plot boxchart of snail domain intensity (yellow probes)
plot_colors = ['b', 'r'];
figure;
for i = 1:length(construct_names)
    xGroup = ones(1,length(domain_intensity_mean_byConstruct{i})).*i; % create xGroup list for boxchart input
    boxchart(xGroup, domain_intensity_mean_byConstruct{i}, "BoxFaceColor", plot_colors(i));
    hold on;
end
xticks(1:length(construct_names))
xticklabels(construct_names)
ylabel('mean yellow probe intensity within snail domain')
title('RNA FISH - yellow probes')

%% histogram plotting (mean yellow intensity within snail domain)

x1 = sig_intensity_byConstruct_linear{1};%domain_intensity_mean_byConstruct{1};
figure;
histogram(x1,'Normalization','probability', 'FaceColor','b', 'EdgeColor','b')
[f_1,xi_1] = ksdensity(x1);
% hold on
% plot(xi_1,f_1, 'b')
xlim([0 10000])

x2 = sig_intensity_byConstruct_linear{2}; %domain_intensity_mean_byConstruct{2};
hold on;
histogram(x2,'Normalization','probability', 'FaceColor','r','EdgeColor','r')
[f_2,xi_2] = ksdensity(x2);
% hold on
% plot(xi_2,f_2, 'r')
xlim([0 10000])

xlabel('mean signal intensity (A.U.)')
ylabel('signal count (#)')

%%

x1 = domain_intensity_distribution_byConstruct{1}(domain_intensity_distribution_byConstruct{1}>3000);
figure;
histogram(x1,'Normalization','count', 'FaceColor','b', 'EdgeColor','b')
[f_1,xi_1] = ksdensity(x1);
% hold on
% plot(xi_1,f_1, 'b')
xlim([3000 30000])
% ylim([0 0.0001])

x2 = domain_intensity_distribution_byConstruct{2}(domain_intensity_distribution_byConstruct{2}>3000);
hold on;
histogram(x2,'Normalization','count', 'FaceColor','r','EdgeColor','r')
[f_2,xi_2] = ksdensity(x2);
% hold on
% plot(xi_2,f_2, 'r')
xlim([3000 30000])
% ylim([0 0.0001])

legend({'enhUP' 'enhDOWN'})