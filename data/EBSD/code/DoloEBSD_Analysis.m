%% Import Script for EBSD Data
%
% This script was automatically created by the import wizard. You should
% run the whoole script or parts of it in order to import your data. There
% is no problem in making any changes to this script.

close all
clear
clc

saveFlag = 1;
reduceF = 1;
arLimit = 1.5;


exp = {'MD6_258_6_250nm'};
% exp = {'MD6_258_6_250nm_cropDonut','MD6_258_6_250nm'};
% exp = {'SB1_ONC_250nm_crop','SB1_ONC_250nm'};
% exp = {'SB1_CON_250nm_cropCircle','SB1_CON_250nm'};

for i = 1:size(exp,1)

    close all
    ename = exp{i,1}

    if contains(ename,'ONC')
        poly = [236.0469 59.0279; 277.4257 91.8187; 307.8744 111.3370;
                339.8844 108.9948; 420.3000 154.2774; 456.9944 178.4801;
                472.6091 194.0948; 481.9779 205.8058; 549.9017 240.1581;
                568.6393 247.9654; 606.1145 280.7562; 670.1347 338.5305;
                738.8392 182.3838; 465.5825 0; 302.4092 0];
    else
        poly = [];
    end

    %% Specify Crystal and Specimen Symmetries
    
    % crystal symmetry
    CS = {... 
      'notIndexed',...
      crystalSymmetry('3', [4.8 4.8 16], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Dolomite', 'color', [0.53 0.81 0.98])};
    
    % plotting convention
    setMTEXpref('xAxisDirection','east');
    setMTEXpref('zAxisDirection','intoPlane');
    
    %% Specify File Names
    
    % path to files
    if contains(ename,'MD6_258_6_250nm')
        pname = 'G:\Shared drives\RIDL\_SEM-EBSD\_Data\_Collaborators\Kristin Bergmann\MD6_258_6_250nm\Map 1';
        fname = [pname,'\',ename,'.ctf'];
    elseif contains(ename,'SB1_CON_250nm')  
        pname = 'G:\Shared drives\RIDL\_SEM-EBSD\_Data\_Collaborators\Kristin Bergmann\SB1_CON_250nm';
        fname = [pname,'\',ename,'.ctf'];
    elseif contains(ename,'SB1_ONC_250nm')  
        pname = 'G:\Shared drives\RIDL\_SEM-EBSD\_Data\_Collaborators\Kristin Bergmann\SB1_ONC_250nm';
        fname = [pname,'\',ename,'.ctf'];
    end    
    
    %% Import the Data
    
    % create an EBSD variable containing the data
    ebsd = EBSD.load(fname,CS,'interface','ctf',...
      'convertEuler2SpatialReferenceFrame');
    
    rawebsd = ebsd;
    ebsd = reduce(ebsd,reduceF);
    
    % Load the uncropped map
    if size(exp,2) > 1
    ebsdWhole = EBSD.load([pname,'\',exp{i,2},'.ctf'],CS,'interface','ctf',...
      'convertEuler2SpatialReferenceFrame');
    ebsdWhole = reduce(ebsdWhole,reduceF);
    end

    if isempty(poly)
        poly = [0 0 range(ebsd.x) range(ebsd.y)];
    end

    %% Calculate grains
    
    disp(' ')
    disp('  Constructing grains...')
    
    % Calculate grains using a 10-degree misorientation threshold
    [grains,ebsd.grainId] = calcGrains(ebsd,'angle',10*degree,'boundary','tight');
    disp('  ...Grain calculation complete!')
    
    % Remove wild spikes
    grains(grains.grainSize == 1).phase = 0;
    % And nullify corresponding pixels
    ebsd(grains(grains.grainSize == 1)).phase = 0;
    
    % Smooth grain boundaries
    grains = smooth(grains,3);
    
    
    %% Plot IPF-x map
    
    % Calculate IPF-x colors
    ipfKey = ipfTSLKey(ebsd);
    ipfKey.inversePoleFigureDirection = vector3d.Y;
    IPFy = ipfKey.orientation2color(ebsd('indexed').orientations);
    
    figure
    if size(exp,2) > 1
        plot(ebsdWhole,'facecolor','k')
            hold on
        plot(ebsdWhole('indexed'),ebsdWhole('indexed').bc)
            colormap gray
            hold on
    end
    plot(ebsd('indexed'),IPFy)
        set(findall(gca,'type','hgtransform'),'visible','off')
        set(gca,'color','k')
        legend off
    if saveFlag
        export_fig([pname,'\Images\',fname(length(pname)+2:end-4),'_IPF-y.png'],'-m5','-r500')
    end

        return

        
    figure
    plot(ipfKey)
        set(findall(gcf,'type','text'),'visible','off')
    if saveFlag
        export_fig([pname,'\Images\',fname(length(pname)+2:end-4),'_IPFkey.png'],'-m3','-r300')
    end
    
    if size(exp,2) > 1
        figure
        plot(ebsdWhole('indexed'),ipfKey.orientation2color(ebsdWhole('indexed').orientations))
        set(findall(gca,'type','hgtransform'),'visible','off')
        set(gca,'color','k')
        legend off
        if saveFlag
            export_fig([pname,'\Images\',fname(length(pname)+2:end-4),'_whole_IPF-y.png'],'-m5','-r500')
        end
    end


    %% Calculate reference orientation deviation (GROD) for each pixel, where the reference orientation is the grain long axis

    ebsdtmp = ebsd('indexed');
    grainstmp = grains('indexed');

    condition = grainstmp.aspectRatio > arLimit & grainstmp.grainSize >= 10;
    grainstmp = grainstmp(condition);
    ebsdtmp = ebsd(ismember(ebsd.grainId,grainstmp.id));
    
    disp('  Calculating angle between long axis and c-axis for each pixel...')
    
    % Find the grain that each pixel belongs to
    [~,idx] = ismember(ebsdtmp.grainId,grains.id);
    
    % Get corresponding long axis for each pixel
    ebsdtmp.prop.longAx = nan(length(ebsdtmp),1);
    tmp = grains.longAxis;
    ebsdtmp.longAx = tmp(idx);
    
    % Find vector corresponding to c-axis of each pixel
    ebsdtmp.prop.cVec = nan(length(ebsdtmp),1);
    ebsdtmp.cVec = ebsdtmp.orientations * Miller(0,0,0,1,ebsd('Dolomite').CS);
    
    % Calculate angle between c-axis and long axis for each grain
    ebsdtmp.prop.theta = nan(length(ebsdtmp),1);
    ebsdtmp.theta = angle(ebsdtmp.longAx,ebsdtmp.cVec)/degree;
    
    % Plot on a map (only for elongate grains)
    figure
    if size(exp,2) > 1
        plot(ebsdWhole,'facecolor','k')
            hold on
        plot(ebsdWhole('indexed'),ebsdWhole('indexed').bc)
            colormap gray
            freezeColors
            hold on
    end
    plot(ebsdtmp,ebsdtmp.theta)
        mtexColorMap LaboTeX
        clim([0 90])
        set(gca,'color','k')
        legend off
        % mtexColorbar('title','c-axis to elongation angle (\circ)','location','southoutside')
        hold on
        set(findall(gca,'type','hgtransform'),'visible','off')
    if saveFlag
        export_fig([pname,'\Images\',fname(length(pname)+2:end-4),'_angleMap_areaWeighted.png'],'-m5','-r500')
    end
    
    
    %% Plot angle between c-axis and long axis as a histogram

    ebsdtmp = ebsdtmp(~inpolygon(ebsdtmp,poly));
    grainstmp = grainstmp(~inpolygon(grainstmp,poly));

    % For each pixel (i.e., area-weighted)
    figure
    histogram(ebsdtmp.theta,0:2:90,'normalization','probability','edgecolor','none')
        xlabel('Misorientation angle (\circ)')
        ylabel('%')
        ylim([0 ceil(50*max(get(gca,'ylim')))/50])
        set(gca,'color','none')
    if saveFlag
        saveas(gcf,[pname,'\Images\',fname(length(pname)+2:end-4),'_angleHist_areaWeighted.svg'])
    end
    
        
    %% Now find the crystal direction parallel to the long axis of each grain
    
    % Find out which crystal direction lies parallel to the long axis of each grain
    longAxMillerPx = inv(ebsdtmp.orientations) .* ebsdtmp.longAx;
    
    odf_Miller = calcDensity(longAxMillerPx,'halfwidth',5*degree);
    
    figure
    plot(odf_Miller,'fundamentalRegion',...
        'hemisphere','upper','projection','eangle')
        colormap parula
        % mtexColorbar
        set(findall(gcf,'type','text'),'visible','off')
        set(gca,'color','none')
        hold on
    plot(odf_Miller,'fundamentalRegion',...
        'hemisphere','upper','projection','eangle',...
        'contour',0:max(get(gca,'clim')),'linecolor','k')
        set(gca,'clim',[0 3])
    if saveFlag
        export_fig([pname,'\Images\',fname(length(pname)+2:end-4),'growthAxisIPF_areaWeighted.png'],'-m3','-r300')
    end
    
    %% Plot c-axis pole figure

    % One point per pixel (i.e., area-weighted)
    figure
    plotPDF(ebsdtmp.orientations,Miller(0,0,0,1,ebsd('Do').CS),'antipodal',...
        'hemisphere','lower','projection','earea',...
        'smooth','halfwidth',7.5*degree)
        colormap parula
        % mtexColorbar
        hold on
    plotPDF(ebsdtmp.orientations,Miller(0,0,0,1,ebsd('Do').CS),'antipodal',...
        'hemisphere','lower','projection','earea',...
        'smooth','halfwidth',7.5*degree,'contour',0:max(get(gca,'clim')),'linecolor','k')
        set(findall(gcf,'type','text'),'visible','off')
        set(gca,'color','none')
        clim([0 2])
    if saveFlag
        export_fig([pname,'\Images\',fname(length(pname)+2:end-4),'_PF_allPoints.png'],'-m3','-r300')
    end

    % One point per grain (i.e., not area-weighted)
    figure
    plotPDF(grainstmp.meanOrientation,Miller(0,0,0,1,ebsd('Do').CS),'antipodal',...
        'hemisphere','lower','projection','earea',...
        'smooth','halfwidth',7.5*degree)
        colormap parula
        % mtexColorbar
        hold on
    plotPDF(grainstmp.meanOrientation,Miller(0,0,0,1,ebsd('Do').CS),'antipodal',...
        'hemisphere','lower','projection','earea',...
        'smooth','halfwidth',7.5*degree,'contour',0:max(get(gca,'clim')),'linecolor','k')
        set(findall(gcf,'type','text'),'visible','off')
        set(gca,'color','none')
        clim([0 2])
    if saveFlag
        export_fig([pname,'\Images\',fname(length(pname)+2:end-4),'_PF_onePerGrain.png'],'-m3','-r300')
    end

    N_grains(1) = length(grainstmp);
    
      

    %% NOW DO THE SAME THING FOR EACH GRAIN RATHER THAN FOR EACH PIXEL
    % - THIS EFFECTIVELY REMOVES THE AREA WEIGHTING
    
    
    
    %% Calculate angle between long axis and average [0001]-axis orientation for each grain
    
    disp('  Calculating angle between long axis and mean c-axis for each grain...')

    % For each grain, derive a vector parallel to the c-axis
    % (using the mean grain orientation)
    grainstmp.prop.cVec = nan(length(grainstmp),1);
    grainstmp.cVec = grainstmp.meanOrientation * Miller(0,0,0,1,ebsd('Dolomite').CS);
    
    % Calculate angle between c-axis and long axis for each grain
    grainstmp.prop.theta = angle(grainstmp.longAxis,grainstmp.cVec)/degree;
    
        
    %% Plot angle between c-axis and long axis as a histogram
    
    figure
    histogram(grainstmp.theta,0:2:90,'normalization','probability','edgecolor','none')
        xlabel('Misorientation angle (\circ)')
        ylabel('%')
        ylim([0 ceil(50*max(get(gca,'ylim')))/50])
        set(findall(gca,'type','hgtransform'),'visible','off')
    if saveFlag
        saveas(gcf,[pname,'\Images\',fname(length(pname)+2:end-4),'_angleHist.svg'])
    end
    
    
    %% Now find the crystal direction parallel to the long axis of each grain
    
    % Find out which crystal direction lies parallel to the long axis of each grain
    longAxMiller = inv(grainstmp.meanOrientation) .* grainstmp.longAxis;
    
    odf_Miller = calcDensity(longAxMiller,'halfwidth',5*degree);
    
    figure
    plot(odf_Miller,'fundamentalRegion',...
        'hemisphere','upper','projection','eangle')
        colormap parula
        % mtexColorbar
        set(findall(gcf,'type','text'),'visible','off')
        set(gca,'color','none')
        hold on
    plot(odf_Miller,'fundamentalRegion',...
        'hemisphere','upper','projection','eangle',...
        'contour',0:max(get(gca,'clim')),'linecolor','k')
        clim([0 3])
    if saveFlag
        export_fig([pname,'\Images\',fname(length(pname)+2:end-4),'growthAxisIPF.png'],'-m3','-r300')
    end
  
    N_grains(2) = length(grainstmp);

end



















return
    
    
    %% Calculate misorientation between each pixel and the mean orientation of the grain it belongs to
    
    ebsd.prop.mis2mean = calcGROD(ebsd,grains);
    
    figure
    plot(ebsd,ebsd.mis2mean.angle/degree,'colorrange',[0 40])
        hold on
    plot(grains.boundary,'linewidth',0.1)
        % mtexColorbar
    
    
    %% Find misorientation between each pixel and the long axis of the grain it belongs to
    
    ebsd.prop.GROD = nan(length(ebsd),1);
    ebsd('indexed').prop.GROD = angle(longAxMillerPx,longAxMiller(idx),'antipodal')/degree;
    
    figure
    plot(ebsd,ebsd.GROD,'colorrange',[0 40])
        % mtexColorbar('title','Misorientation angle relative to elongation/growth axis (\circ)')
        hold on
    plot(grains.boundary,'linewidth',0.1)

%%

misoAxis = axis(ebsd('indexed').orientations, grains(idx).meanOrientation);
thetaSigned = sign(misoAxis.xyz(:,3)) .* angle(ebsd('indexed').orientations, grains(idx).meanOrientation)/degree;

figure
plot(ebsd('indexed'),thetaSigned,'colorrange',[-30 30])
    mtexColorMap red2Blue
    % mtexColorbar('title','Misorientation angle relative to elongation/growth axis (\circ)')
    hold on
plot(grains.boundary,'linewidth',0.1)
    

%%
close all

% Get grain and pixel information for example grain
gtmp = grains(grains.id==4735);
ebsdtmp = ebsd(gtmp);

% Find EBSD orientation parallel to grain long axis
longAxTmpPx = inv(ebsdtmp.orientations) .* repmat(gtmp.longAxis,length(ebsdtmp),1);
    longAxTmpPx = project2FundamentalRegion(longAxTmpPx);

% Find mean grain orientation along long axis
longAxTmp = inv(gtmp.meanOrientation) .* gtmp.longAxis;
    longAxTmp = project2FundamentalRegion(longAxTmp);

misoAxis = axis(ebsdtmp.orientations, gtmp.meanOrientation);
thetaSigned = sign(misoAxis.xyz(:,3)) .* angle(ebsdtmp.orientations, gtmp.meanOrientation)/degree;

figure
plotPDF(ebsdtmp.orientations,thetaSigned,longAxTmp,'antipodal',...
    'hemisphere','lower','projection','earea','colorrange',[-20 20])
mtexColorMap red2Blue
    hold on
plotPDF(gtmp.meanOrientation,longAxTmp,'antipodal',...
    'hemisphere','lower','projection','earea')

rot = rotation.map(ebsdtmp.orientations*longAxTmp,gtmp.meanOrientation*longAxTmp);
dot(rot.axis,longAxTmp)
