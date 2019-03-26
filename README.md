This is a GitHub share folder for the analysis in 

[**Integration of multi-omics datasets enables molecular classification of COPD**](<https://erj.ersjournals.com/content/51/5/1701930.long>)

Chuan-Xing Li, Craig E. Wheelock, C. Magnus Sköld, Åsa M. Wheelock
European Respiratory Journal Jan 2018, 1701930; DOI: 10.1183/13993003.01930-2017

The respiratory is under construction and will continually update until finished.
If you have any questions or suggestions, please contact Chuan-xing Li (lichuanxing@gmail.com).

This repo is a generalized version of the multi-omics analysis, which is mainly focused on the investigation of a number of omics, sample size, and their influence in multi-omics based prediction. Both supervised and unsupervised analysis were used with their corresponding parameters, including permutation number in Leave-one-out cross-validation(LOOCV), number of clusters in unsupervised Spectral clustering. The default parameters for [SNFtool](<https://github.com/maxconway/SNFtool>) is K, alpha, and t, and we search the best parameter combinations within the suggested range by the developer. 

The main code is ***./r/analysis.R***

Here is the workflow for the project:
![alt text](https://github.com/clisweden/ERJ_SNF/blob/master/doc/design_2019-03-26.png)