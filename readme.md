This repository contains the code for a submitted paper for review. The code is in Matlab. The steps for reproducing the figures and results in the paper are given below. First, install Matlab and open it. Then, change the current directory of Matlab to 'FTR-main', i.e., the folder containing the code ('run_algo.m' and all the subdirectories).

# 1. Preprocessing

First, download the data files from ADNI/OASIS/NACC into specific 'data' folders. Then, change the current directory to the specific directory (ADNI/OASIS/NACC) under './preprocess' and run the pipeline code. The details for each cohort is given below.

## 1.1. ADNI

To generate preprocessed ADNI data, first create a 'data' subfolder under './preprocess/ADNI' and download the following data files into './preprocess/ADNI/data' from [ADNI](https://ida.loni.usc.edu/).

*UCSFFSX6_08_17_22_16Aug2023.csv*  
*UCSFFSX51_11_08_19_16Aug2023.csv*  
*UCSFFSX_11_02_15_24Sep2023.csv*  
*UCBERKELEYAV45_8mm_02_17_23_16Aug2023.csv*  
*ADNIMERGE_16Aug2023.csv*  
*UPENNBIOMK_MASTER_FINAL_04Oct2023.csv*  

Then, change the current directory to './preprocess/ADNI' (e.g. `cd ./preprocess/ADNI`).

Finally, execute

```
run_ADNI_pipeline;
```

The resulting preprocessed files will be saved in 'FTR-main/input/ADNI'.

## 1.2. OASIS

To generate preprocessed OASIS data, first create a 'data' subfolder under './preprocess/OASIS' and download the following data files into './preprocess/OASIS/data' from [OASIS]([OASIS Brains](https://sites.wustl.edu/oasisbrains/home/oasis-3/)).

*OASIS3_amyloid_centiloid.csv*  
*OASIS3_demographics.csv*  
*OASIS3_Freesurfer_output.csv*   
*OASIS3_UDSa1_participant_demo.csv*  
*OASIS3_UDSb4_cdr.csv*  

Then, change the current directory to './preprocess/OASIS'.

Finally, execute 

```
run_OASIS_pipeline;
```

The resulting preprocessed files will be saved in 'FTR-main/input/OASIS3'.

## 1.3. NACC

To generate preprocessed NACC data, first create a 'data' subfolder under './preprocess/NACC' and download the following data files into './preprocess/NACC/data' from [NACC (SCAN)](https://naccdata.org/data-collection/forms-documentation/biomarker-imaging#researchers%20data%20dictionary%20%20scan%20mri).

*investigator_ftldlbd_nacc69.csv*  
*UCBERKELEY_AMYLOID_MRIFREE_GAAIN_04Apr2025.csv*  
*UCDMRISBM_04Apr2025.csv*  

Then, change the current directory to './preprocess/NACC'.

Finally, execute

```
run_NACC_pipeline;
```

The resulting preprocessed files will be saved in 'FTR-main/input/NACC'.

# 2. Run the algorithm or show the results

Once the data are preprocessed, we can run the algorithm on them or analyze the results. Note that the trajectory and subtype/stage results have been saved. Hence, these results can be directly visualized by loading the preprocessed datasets. 

To either run the algorithm or visualize the results, we need to first change the current directory to 'FTR-main', and then add the following paths in the command window.

```
addpath(genpath("FTR_code/"));
addpath(genpath("utils/"));
addpath(genpath('analysis'));
```

## 2.1. Show the results

To reproduce the figures from the saved results, first create a 'figures/all' folder under FTR-main, then use the functions in the 'analysis' folder with the specific number and letter in the postfix. For example, 

1. To see the trajectories from FTR on all the datasets, execute

```
analyze_trajectory_S12;
```

The results will be saved in './figures/all/S12_*'.

2. To see that the subtypes have distinct clinical characteristics, execute

```
compare_subtype_demo_T1;
```

The table will be shown in the command window.

3. To see that the subtypes have different progression rates of MMSE in the ADNI discovery set, execute

```
compare_subtype_progression_3A;
```

The figure will be saved as '3A.jpg' in './figures/all'.

4. To see that the stages are associated with MMSE in the ADNI discovery set, execute

```
analyze_stages_3C;
```

The figure will be saved as '3C.jpg' in './figures/all'.

## 2.2. Run the algorithm

We can also run the algorithm to regenerate the trajectory and subtype/stage results (it may take hours on the ADNI 82-region dataset). Or, we can just see how FTR recovers trajectory functions in the demo.

1. To see that FTR can recover a sigmoid function or an arbitrary monotonically increasing function using just the CDF, execute the following code in the command window for the sigmoid function
   
   ```
   run_demo_single(1, 'sigmoid');
   ```
   
   and for the arbitrary monotonically increasing function.
   
   ```
   run_demo_single(1, 'arbitrary');
   ```

2. To run FTR on ADNI/OASIS/NACC (it may take hours on the full dimensional dataset), execute the following code for ADNI
   
   ```
   run_algo('ADNI_FSX_HS', 'FTR_MCEM');
   ```
   
   or OASIS
   
   ```
   run_algo('OASIS3_ROD1_HS', 'FTR_MCEM');
   ```
   
   or NACC.
   
   ```
   run_algo('NACC_HS', 'FTR_MCEM');
   ```

The postfix 'HS' represents the dimensionality of the dataset. The postfixes can be 'HS', 'HM', 'LS', 'LM', which correspond to 82 (64), 41 (32), 26 (14), 13 (7) brain regions for ADNI/OASIS (NACC) respectively.

# Contact

For any questions, please contact Yuan Zhou at zhouyuanzxcv@gmail.com for information.
