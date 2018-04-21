function allClumps = SegmentClumps(I)
%%    
    q = 0.06;
    
    W = wiener2(I, [5 5]);
    i1 = 0;
    i2 = 0;
    while ((max(i1, i2) < 150) || (gObj.mu(1) == gObj.mu(2)))
        gObj = gmdistribution.fit(double(W(:)), 2);
        i1 = icdf('Normal', q, gObj.mu(1), sqrt(gObj.Sigma(1)));
        i2 = icdf('Normal', q, gObj.mu(2), sqrt(gObj.Sigma(2)));
    end
    
	minMeanIntensity = max(icdf('Normal', .0001, gObj.mu(1), sqrt(gObj.Sigma(1))), ...
        icdf('Normal', .0001, gObj.mu(2), sqrt(gObj.Sigma(2))));
    
    clump = I <= max(i1, i2);
    clump = bwareaopen(imopen(clump, strel('disk', 5)), 2000);
    
    clumpProps = regionprops(clump, W, 'MeanIntensity', 'PixelIdxList');
    clump(vertcat(clumpProps([clumpProps.MeanIntensity] > minMeanIntensity).PixelIdxList)) = 0;
    
    [L, num] = bwlabel(clump);
    allClumps = cell(num, 1);
    for i = 1: num
        allClumps{i} = ~bwareaopen(~imdilate(L == i, strel('disk', 1)), 20);
    end
    
%     if (isempty(allClumps))
%         wait  = 5;
%     end    