function cytoplasms = SegmentCytoplasms(I, loadedVolumeImages, contourSize, contourArea, allClumps, alpha, beta)
%%
% This section approximates the cytoplasm boundaries.
gridWidth = 8;

includeDetectedNucleus = zeros(length(contourArea), 1);
for i = 1: size(includeDetectedNucleus, 1)
    for j = 1: length(allClumps)
        if (nnz(contourArea{i}(:, :) & allClumps{j}) >= .7 * contourSize{i})
            includeDetectedNucleus(i) = j;
            break
        end
    end
end

nucleusGridSquareRatio = 0.1;
squareArea = gridWidth * gridWidth;

cytoplasms = cell(size(contourSize));
canvas = false(size(allClumps{1}));
for c = 1: length(allClumps)
    nucleiInClump = find(includeDetectedNucleus == c);
    numOfNuclei = length(nucleiInClump);
    if (isempty(nucleiInClump))
        continue
    end

    sClump = regionprops(allClumps{c}, 'BoundingBox');

    W1 = max(floor(sClump(1).BoundingBox(2) - gridWidth), 1);
    W2 = min(ceil(sClump(1).BoundingBox(2) + sClump(1).BoundingBox(4) + gridWidth), size(canvas, 1));
    H1 = max(floor(sClump(1).BoundingBox(1) - gridWidth), 1);
    H2 = min(ceil(sClump(1).BoundingBox(1) + sClump(1).BoundingBox(3) + gridWidth), size(canvas, 2));
    if (W1 == 1)
        W2 = W2 - mod(W2 - W1 + 1, gridWidth);
    else
        W1 = W1 + mod(W2 - W1 + 1, gridWidth);
    end
    if (H1 == 1)
        H2 = H2 - mod(H2 - H1 + 1, gridWidth);
    else
        H1 = H1 + mod(H2 - H1 + 1, gridWidth);
    end

    gridSquareNucleiEffect = zeros(floor((W2 - W1 + 1) / gridWidth), floor((H2 - H1 + 1) / gridWidth), numOfNuclei);

    if (~isempty(loadedVolumeImages))
        gridSquareFocusMeasure = zeros(length(loadedVolumeImages), floor((W2 - W1 + 1) / gridWidth), floor((H2 - H1 + 1) / gridWidth));
        tempVector = zeros(length(loadedVolumeImages), 1);
        for i = 1: (W2 - W1 + 1) / gridWidth
            for j = 1: (H2 - H1 + 1) / gridWidth
                for k = 1: length(loadedVolumeImages)
                    tempGridSquare = loadedVolumeImages{k}((i - 1) * gridWidth + W1: i * gridWidth + W1 - 1, (j - 1) * gridWidth + H1: j * gridWidth + H1 - 1);

                    tempVector(k) = std(double(tempGridSquare(:)));
                end
                if (any(tempVector))
                    tempVector = (tempVector - min(tempVector(:))) / (max(tempVector(:)) - min(tempVector(:)));
                end
                gridSquareFocusMeasure(:, i, j) = tempVector;
            end
        end
    end        

    nucleusGridSquares = cell(numOfNuclei, 1);
    nucleusCenterCoord = cell(numOfNuclei, 1);
    nucleusCenterCoordOrig = cell(numOfNuclei, 1);
    for n = 1: numOfNuclei
        nucleusGridSquares{n} = cell(0);
        nucleusInClumpArea = contourArea{nucleiInClump(n)}(W1: W2, H1: H2);
        nucleusProp = regionprops(nucleusInClumpArea, 'Centroid');
        nucleusCenterCoordOrig{n} = round(nucleusProp.Centroid);
        nucleusCenterCoord{n} = ceil(nucleusProp.Centroid / gridWidth);
        for i = 1: (W2 - W1 + 1) / gridWidth
            for j = 1: (H2 - H1 + 1) / gridWidth
                if (nnz(nucleusInClumpArea((i - 1) * gridWidth + 1: i * gridWidth, (j - 1) * gridWidth + 1: j * gridWidth)) > nucleusGridSquareRatio * squareArea)
                    nucleusGridSquares{n}{end + 1} = [i, j];
                end
            end
        end

        if (isempty(nucleusGridSquares{n}))
            return
        end

        for i = 1: (W2 - W1 + 1) / gridWidth
            for j = 1: (H2 - H1 + 1) / gridWidth
                for k = 1: length(nucleusGridSquares{n})
                    if (~isempty(loadedVolumeImages))
                        gridSquareNucleiEffect(i, j, n) = gridSquareNucleiEffect(i, j, n) + ...
                            exp(-(sqrt(sum((gridSquareFocusMeasure(:, i, j) - gridSquareFocusMeasure(:, nucleusGridSquares{n}{k}(1), nucleusGridSquares{n}{k}(2))) .^ 2))^2 / (2 * alpha^2))) * ...
                            exp(-(norm([i - nucleusGridSquares{n}{k}(1), j - nucleusGridSquares{n}{k}(2)])^2 / (2 * alpha^2)));
                    else
                        gridSquareNucleiEffect(i, j, n) = gridSquareNucleiEffect(i, j, n) + ...
                            exp(-(norm([i - nucleusGridSquares{n}{k}(1), j - nucleusGridSquares{n}{k}(2)])^2 / (2 * alpha^2)));
                    end
                end
                gridSquareNucleiEffect(i, j, n) = gridSquareNucleiEffect(i, j, n) / length(nucleusGridSquares{n});
            end
        end
    end

    sumOfNucleiEffect = sum(gridSquareNucleiEffect, 3);
    shrinkedClump = imresize(allClumps{c}(W1: W2, H1: H2), 1 / gridWidth);
    for n = 1: numOfNuclei
        thisNucleusEffect = beta * gridSquareNucleiEffect(:, :, n) - sumOfNucleiEffect;

        cellBlobs = imfill(thisNucleusEffect > 0.0, 'holes');
        cellBlobs = imopen(cellBlobs & shrinkedClump, strel('diamond', 1));
        
        % Removing unreachable grid squares (pixels in cellBlobs)...
        needToContinue = true;
        while (needToContinue)
            needToContinue = false;
            [pixX, pixY] = find(cellBlobs);
            for i = 1: length(pixX);
                [inBetweenX, inBetweenY] = bresenham(pixX(i), pixY(i), nucleusCenterCoord{n}(2), nucleusCenterCoord{n}(1));
                inBetweenIndexes = sub2ind(size(cellBlobs), inBetweenX, inBetweenY);
                if ~all(cellBlobs(inBetweenIndexes))
                    cellBlobs(pixX(i), pixY(i)) = false;
                    needToContinue = true;
                end
            end
        end

        cytoplasms{nucleiInClump(n)} = canvas;
        cytoplasms{nucleiInClump(n)}(W1: W2, H1: H2) = allClumps{c}(W1: W2, H1: H2) & imresize(cellBlobs, gridWidth);
        [labeledCellBlobs, ~] = bwlabel(cytoplasms{nucleiInClump(n)}(W1: W2, H1: H2), 4);
        cytoplasms{nucleiInClump(n)}(W1: W2, H1: H2) = labeledCellBlobs == labeledCellBlobs(nucleusCenterCoordOrig{n}(2), nucleusCenterCoordOrig{n}(1));
        if (~labeledCellBlobs(nucleusCenterCoordOrig{n}(2), nucleusCenterCoordOrig{n}(1)))
            disp('ERROR!')
        end

    end
end

cytoplasms = {cytoplasms{includeDetectedNucleus ~= 0}}';

%%
% This section shrinks the cytoplasm (refer to the article).
% It is supposed that nucleus center is at most 150 pixels far horizontally
% or vertically from the cell boundary and we discretize it to 0.5 units.
[discUnits, angleSin] = meshgrid((0: 1: 300), sin((0: 2 * pi / 360: 2 * pi - 0.00001)));
xInc = round(discUnits .* angleSin);
[~, angleCos] = meshgrid((0: 1: 300), cos((0: 2 * pi / 360: 2 * pi - 0.00001)));
yInc = round(discUnits .* angleCos);

contourArea = contourArea(includeDetectedNucleus ~= 0);
includeDetectedNucleus = includeDetectedNucleus(includeDetectedNucleus ~= 0);

cytoIsChanging = true(length(cytoplasms), 1);
for X = 1: 20
    if (~any(cytoIsChanging))
        break;
    end

    if (X == 1)

        for i = 1: length(contourArea)
            expanded = imdilate(contourArea{i}, strel('disk', 7, 0));
            expanded(contourArea{i}) = false;
            I(imdilate(contourArea{i}, strel('disk', 1, 0))) = mean(mean(I(expanded)));
        end

        meanNucleiPix = regionprops(logical(sum(reshape([contourArea{:}], size(I, 1), size(I, 2), []), 3)), I, 'MeanIntensity');
        I = max(I, min([meanNucleiPix.MeanIntensity]));

        I = wiener2(I, [5 5]);
        intensityImage = im2double(I);
        foreground = any(reshape([allClumps{:}], size(intensityImage, 1), size(intensityImage, 2), []), 3);
        intensityImage(~foreground) = 1;
    end

    a = 10;
    c = .5;

    for i = 1: length(cytoplasms)
        if (~cytoIsChanging(i))
            continue
        end

        cytoplasmProp = regionprops(cytoplasms{i}, 'BoundingBox', 'PixelIdxList', 'PixelList');
        nucleusCent = regionprops(contourArea{i}, 'Centroid');
        nucX = round(nucleusCent.Centroid(2));
        nucY = round(nucleusCent.Centroid(1));

        x1 = ceil(cytoplasmProp.BoundingBox(2));
        x2 = x1 + ceil(cytoplasmProp.BoundingBox(4)) - 1;
        y1 = ceil(cytoplasmProp.BoundingBox(1));
        y2 = y1 + ceil(cytoplasmProp.BoundingBox(3)) - 1;
        allX = nucX - xInc;
        allY = nucY + yInc;

        inImage = (allX <= size(intensityImage, 1)) & (allX >= 1) & (allY <= size(intensityImage, 2)) & (allY >= 1);
        inCytoplasm = (allX <= x2) & (allX >= x1) & (allY <= y2) & (allY >= y1);
        inRangeInd = sub2ind(size(intensityImage), allX(inCytoplasm), allY(inCytoplasm));
        inCytoplasm(inCytoplasm) = ismember(inRangeInd, cytoplasmProp.PixelIdxList);

        allBoundaryPoints = zeros(360, 2);
        previousBoundaryPoints = zeros(360, 2);
        includeBoundaryPoints = true(360, 1);
        boundaryPointChanged = true(360, 1);

        for r = 1: 360
            if (X == 1)
                firstPixOutCyto = find(~inCytoplasm(r, :), 1);
                lastPixInCyto = firstPixOutCyto - 1;
                rayLength = min(2 * lastPixInCyto, size(allX, 2));
                firstPixOutImg = find(~inImage(r, :), 1);
                if (~isempty(firstPixOutImg))
                    rayLength = min(rayLength, firstPixOutImg - 1);
                end
                rayWeightVector = 1 - abs((1: rayLength) - lastPixInCyto) / lastPixInCyto;
                rayWeightVector = sigmf(rayWeightVector, [a c]);
            else
                firstPixOutCyto = find(~inCytoplasm(r, :), 1);
                lastPixInCyto = firstPixOutCyto - 1;
                rayLength = lastPixInCyto;
                firstPixOutImg = find(~inImage(r, :), 1);
                if (~isempty(firstPixOutImg))
                    rayLength = min(rayLength, firstPixOutImg - 1);
                end
                rayWeightVector = (1: rayLength) / rayLength;
                rayWeightVector = sigmf(rayWeightVector, [a c]);
            end

            rayX = allX(r, 1: rayLength);
            rayY = allY(r, 1: rayLength);
            previousBoundaryPoints(r, 1) = allX(r, lastPixInCyto);
            previousBoundaryPoints(r, 2) = allY(r, lastPixInCyto);
            indexes = sub2ind(size(intensityImage), rayX, rayY);
            rayPixels = intensityImage(indexes);

            tmp = [0, 0, rayPixels, 1, 1];
            diff = zeros(length(rayPixels), 1);
            for xR = 3: length(tmp) - 2
                diff(xR - 2) = sum(tmp(xR - 1)) - sum(tmp(xR + 1));
            end
            [~, pos] = min(rayWeightVector' .* diff);

            if (rayLength - pos <= 5)
                boundaryPointChanged(r) = false;
            end

            [bX, bY] = ind2sub(size(intensityImage), indexes(pos));
            allBoundaryPoints(r, :) = [bX, bY];
        end
        allBoundaryPoints = allBoundaryPoints(includeBoundaryPoints, :);

        fittingPolynomialOrder = 3;
        fittingWindowWidth = 2 * floor(size(allBoundaryPoints, 1) / 8) + 1;

        fittingMargin = floor(fittingWindowWidth / 2);
        smoothX = sgolayfilt([allBoundaryPoints(end - fittingMargin + 1: end, 2); ...
            allBoundaryPoints(:, 2); allBoundaryPoints(1: fittingMargin, 2)], ...
            fittingPolynomialOrder, fittingWindowWidth);
        smoothY = sgolayfilt([allBoundaryPoints(end - fittingMargin + 1: end, 1); ...
            allBoundaryPoints(:, 1); allBoundaryPoints(1: fittingMargin, 1)], ...
            fittingPolynomialOrder, fittingWindowWidth);
        smoothX = smoothX(fittingMargin + 1: end - fittingMargin);
        smoothY = smoothY(fittingMargin + 1: end - fittingMargin);

        distSmoothBoundary = [smoothY, smoothX] - allBoundaryPoints;
        distValue = sqrt(sum(distSmoothBoundary .^ 2, 2));

        newBoundaryPoints = allBoundaryPoints(distValue <= 10, :);
        fittingWindowWidth = 2 * floor(size(newBoundaryPoints, 1) / 16) + 1;
        fittingMargin = floor(fittingWindowWidth / 2);
        smoothX = sgolayfilt([newBoundaryPoints(end - fittingMargin + 1: end, 2); ...
            newBoundaryPoints(:, 2); newBoundaryPoints(1: fittingMargin, 2)], ...
            fittingPolynomialOrder, fittingWindowWidth);
        smoothY = sgolayfilt([newBoundaryPoints(end - fittingMargin + 1: end, 1); ...
            newBoundaryPoints(:, 1); newBoundaryPoints(1: fittingMargin, 1)], ...
            fittingPolynomialOrder, fittingWindowWidth);
        smoothX = smoothX(fittingMargin + 1: end - fittingMargin);
        smoothY = smoothY(fittingMargin + 1: end - fittingMargin);

        newCyto = poly2mask(smoothX, smoothY, size(intensityImage, 1), size(intensityImage, 2));

        addX = round(smoothX);
        addY = round(smoothY);
        addX = min(addX, size(newCyto, 2));
        addY = min(addY, size(newCyto, 1));
        addX = max(addX, 1);
        addY = max(addY, 1);
        newCyto(sub2ind(size(newCyto), addY, addX)) = true;

        newCyto = newCyto & allClumps{includeDetectedNucleus(i)};
        newCyto = imreconstruct(contourArea{i}, newCyto) | contourArea{i};

        newAndPreDiff = (newCyto & ~cytoplasms{i}) | (~newCyto & cytoplasms{i});

        if ((nnz(newAndPreDiff) / nnz(cytoplasms{i}) < 0.01))
            cytoIsChanging(i) = false;
        end

        cytoplasms{i} = newCyto;
    end
end
