function [bgMask,objectRectangles] = fretGetInitialObjectsAndBGmask(imSum)
%fretGetInitialObjectsAndBG does an initial pass object detection and
%background subtraction to identify objects for finer processing.

%Do an initial segmentation of background and foreground
imNorm=mat2gray(imSum,[prctile(imSum(:),5) prctile(imSum(:),99.5)]);
imNorm_s=imadjust(mat2gray(imfilter(imNorm, fspecial('disk'),'symmetric')));
threshInitialSeg=graythresh(imNorm_s(imNorm_s>0 & imNorm_s<1));
bgMask0=imNorm_s<(threshInitialSeg);%
bgMask=imerode(imdilate(bgMask0,strel('disk',5)),strel('disk',20));
%figure(1); imagesc(bgMask)
fgMask=imerode(imdilate(~bgMask0,strel('disk',5)),strel('disk',10));

%Identify objects and bounding rectangles
objects=regionprops(fgMask,'BoundingBox');
% figure;
% imagesc(imSum); hold on;
for i=1:length(objects)
    width=fretGetObjectRectanglePaddingSize;
    objectRectangles{i}=objects(i).BoundingBox+[-1*width -1*width 2*width 2*width];
%     rectangle('Position',objectRectangles{i});
end

end

