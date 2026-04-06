clear all
close all
tic
% %% load images and set parameters 
% 
% % automate images loading
% % import filenames from excel file
% XL_file = "C:\Users\Emilia Leyes Porello\OneDrive - PennO365\Penn Research\Enhancer-Promoter Spacers\Fly imaging\spacer analysis.xlsx";
% sheet = "spacer analysis";
% [~, txtData] = xlsread(XL_file,sheet);
% [numData] = xlsread(XL_file,sheet);
% 
% % start loop to iterate through all files
% for p = 1:29 %1:(size(txtData,1)-1)
%     clearvars -except p txtData
% 
%     fileName = string(txtData(p+1,1));
%     import = string(txtData(p+1,2));
%     disp(fileName);
%     
%     export2 = import;
%     load(append(export2, 'segmentation_lineage.mat'));
 

% % start here and delete "for"  loop to run single image manually
%export2 = 'C:\Users\Emilia Leyes Porello\OneDrive - PennO365\Penn Research\Enhancer-Promoter Spacers\Fly imaging\snaSE_7.6_evePr_truncYellow X 2color x yw - M2\Images and data for all replicates\MAX_snaSE_7.6_evePr_truncYellow_M2_emb03\';
export2 = 'C:\Users\Emilia Leyes Porello\OneDrive - PennO365\Penn Research\Enhancer-Promoter Spacers\Fly imaging\NewDet_E(3+7)-7.6-t.Yel-evePr x 2color - M5\Images and data for all replicates\MAX_E7-7.6-evePr_emb03\';

% load([export2 'trajectories.mat']);
load([export2 'segmentation_lineage.mat']);

% 
start = 1;
% T = start*0.3333 + 0.3333*(1:size(lineage_cx,1));


%% First, get the threshold for each time point 
for i=1:size(nuc_lineage,2) % nuclei index
    disp(i)
    for j=1:size(nuc_lineage,1)
        %  image with a single nucleus
        c_image = convex_image{j,nuc_lineage(j,i)};
%         c_image = imresize(c_image,1); % search area is 10% bigger
        
        % x, y distance translation
        r1 = regionprops(c_image,'PixelList','centroid');  % get pixellist and centroids for the nucleus
        
        % find the coordinates in the entire image (rather than within a
        % single nucleus)
        dx = lineage_cx(j,i) - r1.Centroid(1); dx = round(dx);
        dy = lineage_cy(j,i) - r1.Centroid(2); dy = round(dy);
        
        r1.PixelList(:,1) = r1.PixelList(:,1) + dx;
        r1.PixelList(:,2) = r1.PixelList(:,2) + dy;
        
        % correction for the coordinates (smaller/bigger than 1, 512)
        ddum = find(r1.PixelList(:,1) < 1); % too small
        r1.PixelList(ddum,1) = 1;
        ddum = find(r1.PixelList(:,2) < 1);
        r1.PixelList(ddum,2) = 1;
        ddum = find(r1.PixelList(:,1) > size(images_ms2,2));
        r1.PixelList(ddum,1) = size(images_ms2,2);
        ddum = find(r1.PixelList(:,2) > size(images_ms2,1));
        r1.PixelList(ddum,2) = size(images_ms2,1);
        
        % x and y coordinates of a given nucleus in the entire frame
        xx = r1.PixelList(:,1);
        yy = r1.PixelList(:,2);
        
        % remove the background, by subtracting the 50 percentile value
%         I_pp7 = images_pp7(:,:,j) - prctile(reshape(images_pp7(:,:,j),[1 size(images_ms2,1)*size(images_ms2,2)]),50);
        I_ms2 = images_ms2(:,:,j) - prctile(reshape(images_ms2(:,:,j),[1 size(images_ms2,1)*size(images_ms2,2)]),45);
%         I_ms2 = images_ms2(:,:,j);
        
%         P = zeros(size(images_ms2,1),size(images_ms2,2),'uint16');
%         O = zeros(size(images_ms2,1),size(images_ms2,2),'uint16');
%         O1 = zeros(size(images_ms2,1),size(images_ms2,2),'uint16');
        
        % extract pp7 and ms2 signals from all the pixel indices
        for k=1:size(xx,1)
%             Ps(k) = images_nuc(yy(k),xx(k),j);
            Ms(k) = images_ms2(yy(k),xx(k),j); % ms2 signal
%             Ps(k) = I_pp7(yy(k),xx(k));
%             Ms(k) = I_ms2(yy(k),xx(k));
%             O(yy(k),xx(k)) = I_ms2(yy(k),xx(k));
%             O1(yy(k),xx(k)) = images_ms2(yy(k),xx(k),j);
%             P(yy(k),xx(k),3) = 25500;
            
        end
        
%         Ps = double(Ps); Ms = double(Ms);
        
%         P_ratio(j,i) = max(Ps)/median(Ps); M_ratio(j,i) = max(Ms)/median(Ms);
        
%         O = imcrop(O,[min(xx) min(yy) max(xx)-min(xx)+1 max(yy)-min(yy)+1]);
%         O = imresize(O,2);
%         O1 = imcrop(O1,[min(xx) min(yy) max(xx)-min(xx)+1 max(yy)-min(yy)+1]);
%         O1 = imresize(O1,2);

%         
%         O_bw = im2bw(O,graythresh(O));
%         rpo = regionprops(O_bw,'centroid');
%         s = size(rpo,1);
%         
%         if s>2 || graythresh(O)<0.01% too many objects
%             O_bw = im2bw(O,0.09);
%             rpo = regionprops(O_bw,'centroid');
%             if size(rpo,1)>2
%                 O_bw = im2bw(O,0.1);
%                 rpo = regionprops(O_bw,'centroid');
%                 s = size(rpo,1);
%             else
%                 s = size(rpo,1);
%             end
%         end
            
        
%         P(:,:,2) = I_ms2*1.5;
%         P(:,:,1) = images_ms2(:,:,j)*2;
        
 
        Mm(j,i) = median(Ms);
        Mmax(j,i) = max(Ms);
        M_stddev(j,i) = std2(Ms); % standard deviation of MS2 signal

%         Pm(j,i) = median(Ps);
%         Pmax(j,i) = max(Ps);
%         M(j,i) = max(Ms); % maximum
        
        % edit - Feb2 average of top 5 pixels 
        Msort = sort(Ms,'descend');
        M(j,i) = mean(Msort(1:2));
%         Psort = sort(Ms,'descend');
%         Pj,i) = mean(Msort(1:2));
        
        
%           save(sprintf('%strajectories.mat',export2),'nuc_lineage','lineage_cx','lineage_cy',...
%         'M','T','Mm','Mmax');
    
%         figure(6); plot(T,M(:,i),'bo-');
%         figure(6); plot(M(:,i),'bo-');
        
        clear xx yy Ms Msort
        
%         w = waitforbuttonpress;
%         if w==1
%             continue
%         end
%         clf(1); clf(2); clf(3);
    end
    
    % correct for outlier
    tmp = find(M(:,i) == Mm(:,i));
    for k=1:length(tmp)
        if tmp(k)>1 && tmp(k)<size(M,1) && M(tmp(k)-1,i)>Mm(tmp(k)-1,i) && M(tmp(k)+1,i)>Mm(tmp(k)+1,i)
            M(tmp(k),i) = Mmax(tmp(k),i);
        end
    end
    
    tmp = find(M(:,i) > Mm(:,i));
    if isempty(tmp) == 0
        tmp1 = find(diff(tmp)==1);
        if isempty(tmp1) == 1 
            M(tmp,i) = Mm(tmp,i);
        else
            tmp2 = find(diff(tmp)>2);
            for k=1:length(tmp2)
                if tmp2(k) == length(tmp)-1
                    if (tmp(tmp2(k)+1)-tmp(tmp2(k)) ~=1)
                        M(tmp(tmp2(k)+1),i) = Mm(tmp(tmp2(k)+1),i);
                    end
                else
                    if (tmp(tmp2(k)+2)-tmp(tmp2(k)+1)~=1) && (tmp(tmp2(k)+1)-tmp(tmp2(k)) ~=1)
                        M(tmp(tmp2(k)+1),i) = Mm(tmp(tmp2(k)+1),i);
                    end
                end
            end
        end
    end
    
    
%     hold on 
%     plot(T,M(:,i),'b','linewidth',2); hold off;
%     w = waitforbuttonpress;
%     if w==1
%         continue
%     end
end
Mm = double(Mm); Mmax = double(Mmax);
      save(sprintf('%strajectories.mat',export2),'nuc_lineage','lineage_cx','lineage_cy',...
        'M','T','Mm','Mmax', 'M_stddev');
    
%end