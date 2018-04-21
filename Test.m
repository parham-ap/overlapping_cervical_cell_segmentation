dataDirectory = 'G:\Google Drive\{Datasets}\{challenges}\Overlapping Cervical Cytology Image Segmentation Challenge - ISBI 2015\Testing';

testCases = dir(fullfile(dataDirectory, 'EDF', '*.png'));

fprintf('\nReading cytoplasm ground truth... ');
testResults = cell(length(testCases), 1);
load(fullfile(dataDirectory, 'AnnotationTest.mat'))
fprintf('Done.\n');

%%
% Segmenting nuclei and cell clumps...
tic
fprintf('\nSegmenting nuclei and cell clumps... \n');
imageInfo = cell(length(testCases), 1);
parfor i = 1: length(testCases)
    fprintf('\t%s: ', testCases(i).name(1: end - 4));
    imageName = fullfile(dataDirectory, 'EDF', testCases(i).name);
    imageInfo{i}.loadedEDF = imread(imageName);
    [imageInfo{i}.nuclei, imageInfo{i}.contourArea, imageInfo{i}.contourSize] = SegmentNuclei(imageInfo{i}.loadedEDF);
    fprintf('nuclei segmented, ');
    imageInfo{i}.allClumps = SegmentClumps(imageInfo{i}.loadedEDF);
    fprintf('cell clump segmented.\n');
end
fprintf('Done.\n');
toc


allVolumeImages = [];
%%
% Loading stack images...
% 
% If there are no volume images available, simply comment the code in this
% section. The method for cytoplasm segmentation will work without the
% volume/stack images as well, although it was originally designed to use
% volume images and may not perform as well without them (refer to the
% article).
fprintf('\nLoading stack images... ');
allVolumeImages = cell(length(testCases), 1);
parfor i = 1: length(testCases)
    imageName = fullfile(dataDirectory, 'EDF', testCases(i).name);
    volumeImageList = dir(fullfile(dataDirectory, [imageName(end - 11: end - 4), '_stack'], '*.png'));
    loadedvolumeImages = cell(length(volumeImageList), 1);
    for j = 1: length(volumeImageList)
        loadedvolumeImages{j} = imread(fullfile(dataDirectory, [imageName(end - 11: end - 4), '_stack'], volumeImageList(j).name));
    end
    allVolumeImages{i} = loadedvolumeImages;
end
fprintf('Done.\n');

%%
% Segmenting cytoplasms...
tic
alpha = 1.75;
beta = 20;
fprintf('\nSegmenting overlapping cytoplasms (alpha = %0.2f, beta = %0.2f)...\n', alpha, beta);

parfor i = 1: length(testCases)
    fprintf('\t%s: ', testCases(i).name(1: end - 4));
    cytoplasms = SegmentCytoplasms(imageInfo{i}.loadedEDF, allVolumeImages{i}, imageInfo{i}.contourSize, imageInfo{i}.contourArea, imageInfo{i}.allClumps, alpha, beta);
    testResults{i} = cytoplasms;
    fprintf('cytoplasms segmented.\n');
end
toc

%%
% Evaluating segmentation...
[DSC, FNRo, TPRp, FDRo, stdDSC, stdFNRo, stdTPRp, stdFDRo] = EvaluateSegmentation(AnnotationTest, testResults);
fprintf(['\nDSC  ', char(177), ' std\tFNRo ', char(177), ' std\tTPRp ', char(177), ' std\tFDRo ', char(177), ' std\n%.3f', char(177), '%.3f\t%.3f', char(177), '%.3f\t%.3f', char(177), '%.3f\t%.3f', char(177), '%.3f\n'], ...
    DSC, stdDSC, FNRo, stdFNRo, TPRp, stdTPRp, FDRo, stdFDRo);
