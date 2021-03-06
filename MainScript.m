%% Script for processing edge vs fret cross correlation dynamics
%% Parameters for cell edge parametrization
nFretWindows=60;                                    % Number of windows to use for FRET vs edge correlation measurements
edgeOversamplingParam=5;                            % How many times more points should the edge of the cell be tracked at
nPointsParam=nFretWindows*edgeOversamplingParam;    % Number of points to track on the cell edge
pdSmoothing=3;                                      % Used with imclose to make the selection of points for tracking less dependent on noise or wrinkles in the cell edge
edgeDepthDist=10;                                   % Number of pixels deep for the windows for computing FRET values.
interiorDepthDist=50;                               % Number of pixels deep for the analyzing the cell interior near the front (distance from center of the leading edge)
%%
mainFolder='C:\Test';
dataDir=[mainFolder filesep 'Data'];
%% Align and Process images, get Coordinates and Measurements, and save the result
clear allCoors; clear subdir; clear allTraj;
clear maskBGinitial; clear maskFinal; clear coors;
subdir=getSubdirectories(mainFolder);

allCoors=cell(length(subdir),1);
for folderNum=1:length(subdir)
    
    folder=[mainFolder filesep subdir{folderNum}];
    fprintf('%2i %s\n',i,subdir{folderNum});
    filenames=getFilenames(folder);
    FRETfile=filenames(boolRegExp(filenames,'FRET'));
    CFPfile=filenames(boolRegExp(filenames,'CFP'));
    clear maskBGinitial; clear maskFinal; clear coors;
    % -- image processing --
    for j=1:length(filenames)
        im(:,:,1)=imread([folder filesep FRETfile{j}]);
        im(:,:,2)=imread([folder filesep CFPfile{j}]);
        imSum=sum(double(im),3);
        [maskBGinitial{j},objectRectangles]=fretGetInitialObjectsAndBGmask(imSum);
        imForSeg0=imageSubtractBackgroundWithObjectMaskingAndPadding(imSum,80,60,maskBGinitial{j});
        imFRET=imageSubtractBackgroundWithObjectMaskingAndPadding(im(:,:,1),80,60,maskBGinitial{j});
        imCFP=imageSubtractBackgroundWithObjectMaskingAndPadding(im(:,:,2),80,60,maskBGinitial{j});
        [maskFinal{j},coors{j}]=fretGetCellMasks(imForSeg0,objectRectangles,2000,imCFP,imFRET);
        %coors[1,2 position; 3 area; 4:7 boundingBox (make the smallest rectangle containing cells); 8,9 centroid; 10 (try this one rather than 11) mean (using sum) intensity of ratio image; 11 median intensity of ratio image; 12 mean intensity of CFP; 13 mean intensity of FRET; 14 ??; 15 16
    end
    save([folder filesep 'masks and coors.mat'],'maskBGinitial','maskFinal','coors');
    allCoors{i}=coors;
end
% Track cells
for i=1:length(allCoors)
    allTraj{i}=ultTrackAnnSearch(allCoors{i},'pairrule','fwdbckmtch','maxdisp',100,'verbose',false);
end
save([mainFolder filesep mainFolder{sub} filesep 'allCoors and allTraj.mat'],'allCoors','allTraj');
%% Compute and Save Cell Edge Dynamics Data
clear edgeData;
cellCount=0;
for snum=1:length(subdir)
    subfolder=[mainFolder filesep subdir{snum}];
    fprintf('%s\n',subfolder);
    filenames=getFilenames(subfolder,'.tif$');
    FRETfile=filenames(boolRegExp(filenames,'FRET'));
    CFPfile=filenames(boolRegExp(filenames,'CFP'));
    try
        timeData=fretReadTimeDataFromMetadata(subfolder,[],numBefore);
    catch
        fprintf('%-8s %2i,%2i - could not read time data.\n',gtpases{gnum},enum,snum);
        timeData=(0:5:5*(length(filenames)-1))';
    end
    load([subfolder filesep 'masks and coors.mat'],'maskBGinitial','maskFinal','coors');
    fprintf('%i trajectories.\n',length(traj));
    for trajNum=1:length(traj)
        if size(traj{trajNum},1)>=20
            cellCount=cellCount+1;
            frameCount=0;
            clear edgeCoors; clear edgeCoorsSmoothed; clear windowCoors; clear protvals; clear fretvals; clear indFrontCenter; clear windowCoors
            startFrame=traj{trajNum}(1,end);
            endFrame=traj{trajNum}(end,end);
            for imageNum=startFrame:endFrame
                frameCount=frameCount+1;
                
                % Read images and compute FRET values
                imFRET=imread([folder filesep FRETfile{j}]);
                imCFP=imread([folder filesep CFPfile{j}]);
                imFRET=imageSubtractBackgroundWithObjectMaskingAndPadding(im(:,:,1),80,60,maskBGinitial{imageNum});
                imCFP=imageSubtractBackgroundWithObjectMaskingAndPadding(im(:,:,2),80,60,maskBGinitial{imageNum});
                
                % Identify the correct cell and make a one cell mask
                objects=regionprops(maskFinal{imageNum},'PixelIdxList','PixelList','Centroid','BoundingBox');
                thisTraj=traj{trajNum};
                cellCent=round(thisTraj(find(thisTraj(:,end)==imageNum,1),1:2));
                cellNum=find(round(arrayfun(@(x) x.Centroid(1),objects))==cellCent(1) & round(arrayfun(@(x) x.Centroid(2),objects))==cellCent(2));
                thisMask=false(size(maskFinal{imageNum}));
                thisMask(objects(cellNum).PixelIdxList)=true;
                
                % Parametrize cell edge and compute protrusion values
                if frameCount==1
                    [edgeCoors{frameCount},edgeCoorsSmoothed{frameCount}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing));
                else
                    [edgeCoors{frameCount},edgeCoorsSmoothed{frameCount}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing),edgeCoors{frameCount-1});
                    windowCoors{frameCount}=edgeCoors{frameCount}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);
                    protvals(:,frameCount-1)=vect(computeProtrusionValues(edgeCoors{frameCount-1},edgeCoorsSmoothed{frameCount-1},edgeCoors{frameCount}));
                end
                windowCoors{frameCount}=edgeCoors{frameCount}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);
                
                % Compute edge regions and local fret ratio values
                labelMask=getWindowLabelMap(thisMask,windowCoors{frameCount},edgeDepthDist);
                for k=1:size(windowCoors{frameCount},1)
                    fretvals(k,frameCount)=sum(imFRET(labelMask==k))/sum(imCFP(labelMask==k));
                end
            end
            protvalsWindow=zeros(size(fretvals)-[0 1]);
            for k=1:edgeOversamplingParam
                protvalsWindow=protvalsWindow+protvals(k:edgeOversamplingParam:end,:);
            end
            % Normalize protrusion data by dividing by the edgeOversampling factor and dividing by the time between images.
            protvalsWindow=(protvalsWindow/edgeOversamplingParam)./repmat(diff(timeData(startFrame:endFrame))',[nFretWindows 1]);
            
            % Fit a bleaching parameter
            bleachCorrectionLinear=fitline(1:size(fretvals,2),nanmean(fretvals,1));
            bleachCorrectionExp=fitline(1:size(fretvals,2),log(nanmean(fretvals,1)));
            
            % Build an edgeData structure
            edgeData(cellCount).edgeCoors=edgeCoors;
            edgeData(cellCount).edgeCoorsSmoothed=edgeCoorsSmoothed;
            edgeData(cellCount).windowCoors=windowCoors;
            edgeData(cellCount).cellArea=traj{trajNum}(:,3);
            edgeData(cellCount).protvalsAll=protvals;
            edgeData(cellCount).protvalsWindow=protvalsWindow;
            edgeData(cellCount).fretvals=fretvals;
            edgeData(cellCount).rootDir=mainFolder;
            edgeData(cellCount).folder=subfolder;
            edgeData(cellCount).trajNum=trajNum;
            edgeData(cellCount).startFrame=startFrame;
            edgeData(cellCount).endFrame=endFrame;
            edgeData(cellCount).bleachCorrectionLinear=bleachCorrectionLinear;
            edgeData(cellCount).bleachCorrectionExp=bleachCorrectionExp;
            bleachCorrAllValues.m(cellCount,1)=edgeData(cellCount).bleachCorrectionLinear.m;
            bleachCorrAllValues.b(cellCount,1)=edgeData(cellCount).bleachCorrectionLinear.b;
        end
    end
    % Compare bleaching correction values
    bleachCorrAllValues.rel=bleachCorrAllValues.m./bleachCorrAllValues.b;
end
save([dataDir filesep 'Cell edge data.mat'],'edgeData','bleachCorrAllValues');  % Save the data for this experimental group