% Perform cross-sectional transect sampling of Mg/Ca element ratios in
% natural carbonate crystals from Ocean Drilling Program Site 807
% measured via electron probe microanalysis at UC Davis

% Copyright 2019 Elizabeth Mitnick

% Lammers & Mitnick, Magnesian calcite solid solution thermodynamics
%       inferred from authigenic deep-sea carbonate
%       accessible @ https://doi.org/10.1016/j.gca.2019.01.006

clear
close all


maxTransectLength = 8; % transect length in pixels
numtransects = 8; % number of transects per particle
particleNum = 1; % particle ID
mapNo = 6;  % map ID
sf_quant = 111.312; % scaling factor


load microprobe_2017 % load structured data
path = Microprobe.Nov29.(strcat('area',num2str(mapNo))); % identify path
ca = path.ca; % Ca element map
mg = path.mg; % Mg element map
bse = path.bse; % Backscattered image element map
bse = im2bw(bse); % Make BSE map black and white
ratio = sf_quant.*mg./ca; % Calculate corrected Mg/Ca ratio
[L W] = size(ca); % determine size of map

% Background filtering
for i = 1:L
    for j = 1:W
        if bse(i,j) == 0
            ratio(i,j) = 0;
            mg(i,j) = 0;
            ca(i,j) = 0;
        end
        
        if ca(i,j) < 1000
            ratio(i,j) = 0;
            ca(i,j) = 0;mg(i,j) = 0;
        end
        
    end
end



sf2 = 0.32; % scaling factor
p2 = figure; % initialize figure
imshow(ratio.*sf2) % display scaled Mg/Ca map
colormap(p2,parula) % specify colormap

%%

summMatr = [nan nan];
numCenters = Microprobe.Nov29.(strcat('area',num2str(mapNo))).authCenters; % ID center position of particle
[p q] = size(numCenters);

% Starting at the center of the particle, record pixel Mg/Ca corrected
% value
for in3 = particleNum
    %coordinates of auth center point
    cx = numCenters(in3,1);%X-coord
    cy = numCenters(in3,2);%Y-coord
    
    transect = nan(maxTransectLength,2,numtransects);
    dx = 0.5;
    
    for in4 = 1:maxTransectLength
        in5=in4;
        iii=in4;
        %perpendicular
        
        if cy+in5 <= L
            transect(iii,1,1) = ratio(cy+in5,cx);
            perp2Bottom = cy+in5;
        end
        if cy-in5 >= 1
            transect(iii,1,2) = ratio(cy-in5,cx);
            perp2Top = cy-in5;
        end
        if cx+in4 <= W
            transect(iii,1,3) = ratio(cy,cx+in4);
            perp2Rt = cx + in4;
        end
        if cx-in4 >= 1
            transect(iii,1,4) = ratio(cy,cx-in4);
            perp2Lft = cx - in4;
        end
        
        %cross-lines
        if cx+in4 <= W && cy+in5 <= L
            transect(iii,1,5) = ratio(cy+in5,cx+in4);
            cross2Rt = cx+in4;cross2Bottom = cy+in5;
        end
        if cx+in4 <= W && cy-in5 >= 1
            transect(iii,1,6) = ratio(cy-in5,cx+in4);
            cross2Rt = cx+in4;cross2Top = cy-in5;
        end
        if cx-in4 >= 1 && cy+in5 <= L
            transect(iii,1,7) = ratio(cy+in5,cx-in4);
            cross2Lft = cx-in4;cross2Bottom = cy+in5;
        end
        if cx-in4 >= 1 && cy-in5 >= 1
            transect(iii,1,8) = ratio(cy-in5,cx-in4);
            cross2Lft = cx-in4;cross2Top = cy-in5;
        end
        
        
        
        
        
    end
    
    %transect = position vs val for all transects (=numtransects)
    for ind5 = 1:numtransects %for all transects
        c=0; %start counter at zero
        for ind6 = 1:maxTransectLength %for each pixel length along transect
            
            %if transect bridges background, all pixels from
            %background and on are removed from analysis
            if transect(ind6,1,ind5) == 0 %if pixel value = 0 (background)
                if ind6 == maxTransectLength %if at the end of the transect,
                    continue %then quit this for-loop
                end
                transect(ind6+1,1,ind5) = 0;%if not, set the next pixel value to 0
            end
            
            %if pixel is non-zero (i.e, not background), start counting
            if transect(ind6,1,ind5) > 0
                c = c+1; %number of continous non-zero pixels
            end
            r = c/2; %length (in microns) of transect with background removed
            lArr = 0:0.5:r-0.5;
            lArr = sort(lArr,'descend');
            transect(1:c,2,ind5) = lArr';
            
        end
        
    end
    
    
    
    %Draw transect lines
    line([cx cx],[cy perp2Bottom],'LineWidth',1,'Color',rgb('OrangeRed'));
    %center down
    hold on
    line([cx cx],[cy perp2Top],'LineWidth',1,'Color',rgb('OrangeRed'));
    %center up
    hold on
    line([cx perp2Rt],[cy cy],'LineWidth',1,'Color',rgb('OrangeRed'));
    %center right
    hold on
    line([cx perp2Lft],[cy cy],'LineWidth',1,'Color',rgb('OrangeRed'));
    %center left
    hold on
    
    line([cx cross2Rt],[cy cross2Bottom],'LineWidth',1,'Color',rgb('OrangeRed'));
    hold on
    line([cx cross2Rt],[cy cross2Top],'LineWidth',1,'Color',rgb('OrangeRed'));
    hold on
    line([cx cross2Lft],[cy cross2Bottom],'LineWidth',1,'Color',rgb('OrangeRed'));
    hold on
    line([cx cross2Lft],[cy cross2Top],'LineWidth',1,'Color',rgb('OrangeRed'));
    hold on
    
end

%Plot Mg/Ca molar ratios as function of distance from particle exterior
figure
set(gca,'fontsize',14)
for i = 1:numtransects
    distFrEdge = transect(:,2,i);
    MgCa = transect(:,1,i);
    plot(distFrEdge,MgCa,'k')
    hold on
    summMatr = [summMatr; transect(:,:,i)];
end
ylim([0 6])
set(gca,'fontsize',14)
xlabel('Distance from crystal exterior (\mum)','fontsize',16)
ylabel(' Mg/Ca Molar Ratio (mmol/mol)','fontsize',16)
title('Mg/Ca Ratio of Calcite Crystal Measured via EPMA')
