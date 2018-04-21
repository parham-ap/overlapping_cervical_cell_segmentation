function [nuclei, contourArea, contourSize] = SegmentNuclei(I)
%%
    if (size(I, 3) == 3)
        I = rgb2gray(I);
    end
    minAcceptableBoundDiff = 15;
    
    cells = struct('MinSize', 110, 'MinMean', 60, 'MaxMean', 150, 'MinSolidity', 0.9);
    
    lowN = floor(cells.MinMean / 10) * 10;
    highN = ceil(cells.MaxMean / 10) * 10;
    
    
    I = wiener2(I, [5 5]);
    nuclei = zeros(size(I));

    allPixels = length(I(:));

    for thresh = lowN: 10: highN
        binaryImage = I <= thresh;
        if sum(binaryImage(:)) > allPixels / 5
            break
        end

        blobs = bwlabel(binaryImage);
        regProp = regionprops(blobs, 'Area', 'Solidity', 'PixelIdxList');
        addTheseRegions = true(length(regProp), 1);
        removeHighSTDTooConcaveTooSmallTooLargeBlobs = ([regProp.Area] < cells.MinSize) | ...
                                                        ([regProp.Solidity] < cells.MinSolidity);
        
        addTheseRegions(removeHighSTDTooConcaveTooSmallTooLargeBlobs) = false;

        pixelsAlreadyInNuclei = find(blobs ~= 0 & nuclei ~= 0);
        blobsAlreadyInNuclei = unique(blobs(pixelsAlreadyInNuclei));
        
        
        nuclei = bwlabel(nuclei);
        nucRegProp = regionprops(nuclei, 'Solidity', 'PixelIdxList');
        if (~isempty(blobsAlreadyInNuclei))
            for j = blobsAlreadyInNuclei'
                intersectWithThese = unique(nuclei(blobs == j));
                if (regProp(j).Solidity < max([nucRegProp(intersectWithThese(intersectWithThese > 0)).Solidity]))
                    addTheseRegions(j) = false;
                end
            end
        end
        
        nuclei(cat(1, regProp(addTheseRegions).PixelIdxList)) = 1;
        nuclei = logical(nuclei);
        
    end
    nuclei = bwareaopen(imfill(nuclei, 'holes'), floor(cells.MinSize), 8);
    
    dilatedSeg = imdilate(nuclei, strel('disk', 1));
    [regionsLabel, ~] = bwlabel(nuclei);
    [dilatedRegionsLabel, numOfDilatedRegions] = bwlabel(dilatedSeg);
    for l = 1: numOfDilatedRegions
        if (length(unique(regionsLabel(dilatedRegionsLabel == l))) >= 3)
            dilatedSeg(dilatedRegionsLabel == l) = nuclei(dilatedRegionsLabel == l);
        end
    end
    nuclei = dilatedSeg;
    
    [labels, totalLabels] = bwlabel(nuclei);
    meanDiff = zeros(totalLabels, 1);
    for i = 1: totalLabels
        dilatedNucleus = imdilate(labels == i, strel('disk', 3, 0));
        meanDiff(i) = mean(I(dilatedNucleus & (labels == 0))) - mean(I(labels == i));
    end
    regProp = regionprops(labels, 'PixelIdxList');
    removeNotHighMeanDiff = meanDiff < minAcceptableBoundDiff;
    nuclei(cat(1, regProp(removeNotHighMeanDiff).PixelIdxList)) = 0;
    
    [L, num] = bwlabel(nuclei);
    [contourArea, contourSize] = deal(cell(num, 1));
    for i = 1: num
        contourArea{i} = L == i;
        contourSize{i} = length(find(contourArea{i}));
    end