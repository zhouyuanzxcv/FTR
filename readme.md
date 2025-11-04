The code is in Matlab. First, install Matlab and open it. Then, change the current directory of Matlab to 'FTR', i.e., the folder containing the code ('run_algo.m' and all the subdirectories.)

# Preprocessing

First, download the data files from ADNI/OASIS/NACC into specific folders. Then, change the current directory to the specific directory (ADNI/OASIS/NACC) under './preprocess' and the run the code.

## ADNI

To generate preprocessed ADNI data, first download the following data files into './preprocess/ADNI/data' from ADNI

UCSFFSX6_08_17_22_16Aug2023.csv
UCSFFSX51_11_08_19_16Aug2023.csv
UCSFFSX_11_02_15_24Sep2023.csv
UCBERKELEYAV45_8mm_02_17_23_16Aug2023.csv
ADNIMERGE_16Aug2023.csv
UPENNBIOMK_MASTER_FINAL_04Oct2023.csv

Then, change the current directory to 'FTR/preprocess/ADNI' (e.g. `cd ./preprocess/ADNI`).

Finally, execute

```run_ADNI_pipeline;```

The resulting preprocessed files will be saved in 'FTR/input/ADNI'.

## OASIS

To generate preprocessed OASIS data, first download the following data files into './preprocess/OASIS/data' from OASIS

OASIS3_amyloid_centiloid.csv
OASIS3_demographics.csv
OASIS3_Freesurfer_output.csv
OASIS3_PUP.csv
OASIS3_UDSa1_participant_demo.csv
OASIS3_UDSb4_cdr.csv

Then, change the current directory to 'FTR/preprocess/OASIS'.

Finally, execute 

```run_OASIS_pipeline;```

# Run the algorithm or show the results

To run the algorithm or analyze the results, change the current directory to 'FTR'. Then, add the following paths in the command window

```
addpath(genpath("FTR_code/"));
addpath(genpath("utils/"));
addpath(genpath('analysis'));
addpath(genpath('competing_methods'));
addpath('figures');
```

## Show the results

The trajectory and subtype/stage results have been saved. To reproduce the figures, use the functions in the 'analysis' folder with the specific number and letter in the postfix. For example, 

1. To see the trajectories from FTR on all the datasets, execute

```analyze_trajectory_S12;```

The results will be saved in './figures/all/S12_*'.

2. To see that the subtypes have distinct clinical characteristics, execute

```compare_subtype_demo_T1;```

3. To see that the subtypes have different progression rates of MMSE in the ADNI discovery set, execute

```compare_subtype_progression_3A;```

4. To see that the stages are associated with MMSE in the ADNI discovery set, execute

```analyze_stages_3C;```

## Run the algorithm

1. To see that FTR can recover a sigmoid function or an arbitrary monotonically increasing function using just the CDF, execute the following code in the command window:
   
   ```
   run_demo_single(1, 'sigmoid');
   run_demo_single(1, 'arbitrary');
   ```
2. To run FTR on ADNI or OASIS (it will take several minutes), execute:
   
   ```
   run_algo('ADNI_FSX_HS', 'FTR_MCEM');
   run_algo('OASIS3_ROD1_HS', 'FTR_MCEM');
   ```
