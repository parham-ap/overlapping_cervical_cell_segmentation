# Cervical Nuecleus Detection and Segmentation, and Overlapping Cytoplasm Segmentation

This is the code of the method in *A Framework for Nucleus and Overlapping Cytoplasm Segmentation in Cervical Cytology Extended Depth of Field and Volume Images* published in Computerized Medical Imaging and Graphics (CMIG) in 2017. The is the final version of the method in *A New Approach to Detect and Segment Overlapping Cells in Multi-Layer Cervical Cell Volume Images* published in International Symposium on Biomedical Imaging (ISBI) 2016.

## Some Notes
* The codes don't include comprehensive comments (except for the EvaluationSegmentation.m). To run the codes on the ISBI 2015 Overlapping Cervical Cytology Image Segmentation Challenge run Test.m after modifying the the path to the dataset directory (first line). You will need to have the ground truth for the test dataset in the challenge to evaluate the segmentation. Test.m uses the same folder structure that was used by the challenge organizers.
* Unfortunately, the original evaluation code provided during the challenges had a bug that reports wrong results. The publication in CMIG explains this and uses this evaluation code to evaluate the segmentation performance. EvaluateSegmentation.zip has the evaluation code along with an example. Moreover, the code reports False Discovery Rate at object level instead of False Positive Rate at pixel level (refer to the article).
* The code can be used also on datasets that do not have volume (stack) images. To use the code on those datasets you need to pass an empty arraty as the second parameter of SegmentCytoplasms function. The nucleus detection/segmentation and cell clump segmentation methods do not require volume images.
* If you use parts of this code, please consider citing the works published in ISBI 2016 and CMIG 2017:

```
@inproceedings{phoulady2016segmentation, 
	author={{Ahmady Phoulady}, Hady and Goldgof, Dmitry B. and Hall, Lawrence O. and Mouton, Peter R.}, 
	booktitle={2016 IEEE 13th International Symposium on Biomedical Imaging (ISBI)}, 
	title={A new approach to detect and segment overlapping cells in multi-layer cervical cell volume images}, 
	year={2016}, 
	pages={201--204}, 
	keywords={Approximation algorithms;Cervical cancer;Gaussian mixture model;Image segmentation;Testing;Training data;cytoplasm segmentation;multi-layer cytology;nucleus detection;overlapping cell segmentation}, 
	doi={10.1109/ISBI.2016.7493244}, 
	month={April}
}
```

```
@article{phoulady2017framework,
	title={A framework for nucleus and overlapping cytoplasm segmentation in cervical cytology extended depth of field and volume images},
	author={{Ahmady Phoulady}, Hady and Goldgof, Dmitry and Hall, Lawrence O and Mouton, Peter R},
	journal={Computerized Medical Imaging and Graphics},
	volume={59},
	pages={38--49},
	year={2017},
	publisher={Elsevier}
}
```