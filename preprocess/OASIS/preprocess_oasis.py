import pandas as pd
import warnings
import numpy as np
import scipy
from typing import *
import os

version = ['HM', 'HS', 'LM', 'LS']


def sort_except_first(lst):
    first_element = lst[0]
    sorted_elements = sorted(lst[1:], key=str.lower)
    return [first_element] + sorted_elements


df_mapping = pd.read_csv('./FSX_biomarker_mapping.csv')

biomarker_ver_map = {element: [item.replace('FSX-', '') if item.startswith('FSX-') else item for item in
                               sort_except_first(df_mapping[element].unique().tolist())] for element in version}

map_LM_HM = {element: [item.replace('FSX-', '') for item in
                       df_mapping[df_mapping['LM'] == 'FSX-' + element]['HM'].unique().tolist()] for element in
             biomarker_ver_map['LM'][1:]}

map_LS_HS = {element: [item.replace('FSX-', '') for item in
                       df_mapping[df_mapping['LS'] == 'FSX-' + element]['HS'].unique().tolist()] for element in
             biomarker_ver_map['LS'][1:]}


def delete_subject_discrete(df, column_name, indicators):
    df_delete = df[~df[column_name].isin(indicators)]
    subjects = df_delete["OASISID"].unique()
    df_deleted = df[~df["OASISID"].isin(subjects)]
    return df_deleted


def choose_subject_discrete(df, column_name, indicators):
    df_chosen = df[df[column_name].isin(indicators)]
    subjects = df_chosen["OASISID"].unique()
    df_chosen = df[df["OASISID"].isin(subjects)]
    return df_chosen

def delete_duplicate_row(df, columns_to_compare):
    for i in range(len(df) - 1):
        if (df.loc[i, columns_to_compare] == df.loc[i + 1, columns_to_compare]).all():
            df.loc[i, columns_to_compare] = np.nan

    return df

def get_volume_data(subjects, output_postfix):
    def strip_prefix(x):
        return int(x[4:])

    def get_years(x):
        return int(x.split("d")[1]) / 365
    
    def get_fldstrength(x):
        if x == 'freesurfer-Linux-centos6_x86_64-stable-pub-v5.3.0-HCP-patch':
            return 3
        else:
            return 1.5

    df_volume = pd.read_csv("./data/OASIS3_Freesurfer_output.csv")
    df_volume = df_volume[df_volume["Subject"].isin(subjects)]
    df_volume["RID"] = df_volume["Subject"].apply(strip_prefix)
    df_volume["years"] = df_volume["MR_session"].apply(get_years)
    df_volume["ICV"] = df_volume["IntraCranialVol"]
    df_volume['FLDSTRENG'] = df_volume['version'].apply(lambda x: 3 if x == 'freesurfer-Linux-centos6_x86_64-stable-pub-v5.3.0-HCP-patch' else 1.5) 

    biomarkers = df_mapping['BIOMARKER'].unique().tolist()
    df_biomarker = df_volume[["RID", "years", "FLDSTRENG"]]
    for region in biomarkers:
        if region == 'ICV':
            df_biomarker['ICV'] = df_volume["ICV"]
        elif region.split('-', 1)[1][0].isupper():
            df_biomarker[region] = df_volume[f"{region}_volume"]
        else:
            if region.split('-', 1)[0] == 'Left':
                df_biomarker[region] = df_volume[f"lh_{region.split('-', 1)[1]}_volume"]
            else:
                df_biomarker[region] = df_volume[f"rh_{region.split('-', 1)[1]}_volume"]

    for ver in version:
        df_ver = df_biomarker[["RID", "years", "FLDSTRENG"]]
        biomarkers_ver = biomarker_ver_map[ver]
        for region in biomarkers_ver:
            if region != 'ICV':
                source_col = df_mapping[df_mapping[ver] == 'FSX-' + region]['BIOMARKER'].unique().tolist()
            else:
                source_col = df_mapping[df_mapping[ver] == region]['BIOMARKER'].unique().tolist()
            df_ver.loc[:, region] = df_biomarker[source_col].sum(axis=1)
        df_ver = pd.concat([df_ver[df_ver.columns[:3]],
                            df_ver[df_ver.columns[3:]].reindex(sorted(df_ver.columns[3:], key=str.lower), axis=1)],
                           axis=1)

        df_ver.to_csv(f"./Result/volume_{output_postfix}_{ver}.csv", index=0)


# def get_abeta_form():
#     df_abeta = pd.read_csv('./data/OASIS3_PUP.csv')
#     df_abeta['OASISID'] = df_abeta['PUP_PUPTIMECOURSEDATA ID'].apply(lambda x: x[:8])
#     df_abeta['RID'] = df_abeta['PUP_PUPTIMECOURSEDATA ID'].apply(lambda x: int(x[4:8]))
#     df_abeta['years'] = df_abeta['PUP_PUPTIMECOURSEDATA ID'].apply(lambda x: int(x.split("d")[1]) / 365)

#     df_abeta_PIB = df_abeta[(df_abeta['tracer'] == 'PIB') & (df_abeta['PET TC QC Status'] == 'Passed')]

#     df_abeta_name = pd.read_csv('abeta_biomarker.tsv', sep='\t')

#     df_result = df_abeta_PIB[['RID', 'years']]

#     df_result[df_abeta_name['BIOMARKER'].tolist()] = df_abeta_PIB[
#         ['PET_fSUVR_TOT_' + item for item in df_abeta_name['ABETACOL'].tolist()]]
#     df_result[['Left-' + item for item in df_abeta_name['BIOMARKER'].tolist()]] = df_abeta_PIB[
#         ['PET_fSUVR_L_' + item for item in df_abeta_name['ABETACOL'].tolist()]]
#     df_result[['Right-' + item for item in df_abeta_name['BIOMARKER'].tolist()]] = df_abeta_PIB[
#         ['PET_fSUVR_R_' + item for item in df_abeta_name['ABETACOL'].tolist()]]
    
#     os.makedirs('./abeta', exist_ok=True)

#     df_result.to_csv('./abeta/abeta_all.csv', index=False)

def get_input_form(output_path, group, hemisphere_status='HM', isextra=False, match=False, val=False):
    def get_years(x):
        return int(x.split("d")[1]) / 365

    def strip_prefix(x):
        return int(x[4:])

    def count_number_4(n):
        str_n = str(n)
        count = str_n.count('4')
        return count

    def cdr2dia(x):
        if isinstance(x, float) and x >= 1:
            return 1
        else:
            return x

    def weighted_sum(weight_lst, item_lst):
        return sum([weight_lst[i] * item_lst[i] for i in range(len(weight_lst))]) / sum(weight_lst)
    
    def match_demo(df, normalize=True, n_matches=2, with_replacement=True):

        from scipy.spatial.distance import cdist
        from sklearn.preprocessing import StandardScaler

        # demo_columns = ['AGE', 'PTGENDER', 'PTEDUCAT', 'APOE', 'CDRTOT', 'MMSE']
        demo_columns = ['AGE', 'PTGENDER']

        if not all(col in df.columns for col in demo_columns):
            raise ValueError("Some demo_columns are not present in the DataFrame.")

        df_group1 = df[df['group'] == 1]
        df_group0 = df[df['group'] == 0]

        if df_group1[demo_columns].isnull().any().any() or df_group0[demo_columns].isnull().any().any():
            df_group1 = df_group1.fillna(df_group1[demo_columns].mean())
            df_group0 = df_group0.fillna(df_group0[demo_columns].mean())
        
        data_group1 = df_group1[demo_columns].values
        data_group0 = df_group0[demo_columns].values

        if data_group1.size == 0 or data_group0.size == 0:
            return pd.DataFrame(columns=df.columns)

        if normalize:
            scaler = StandardScaler()
            all_data = df[demo_columns].values
            scaler.fit(all_data)
            data_group1 = scaler.transform(data_group1)
            data_group0 = scaler.transform(data_group0)

        if with_replacement:
            distances = cdist(data_group1, data_group0, 'euclidean')
            min_distance_indices = np.argsort(distances, axis=1)[:, :n_matches]
            min_distance_indices = np.unique(min_distance_indices.ravel())
            matched_df = df_group0.iloc[min_distance_indices]
        else:
            distances = cdist(data_group1, data_group0, 'euclidean')
            matched_indices = []
            remaining_group0_indices = list(range(len(df_group0)))

            for i in range(len(df_group1)):
                if len(remaining_group0_indices) < n_matches:
                    print(f"Warning: Not enough group0 samples remaining for {n_matches} matches.")
                    break
                    
                current_distances = distances[i, remaining_group0_indices]
                closest_indices = np.argsort(current_distances)[:n_matches]
                selected_indices = [remaining_group0_indices[idx] for idx in closest_indices]
                matched_indices.extend(selected_indices)
                remaining_group0_indices = [idx for idx in remaining_group0_indices if idx not in selected_indices]
            
                matched_df = df_group0.iloc[matched_indices]

        df_control = matched_df.reset_index(drop=True)
        df_case = df[df['group'] == 1]

        df_result = pd.concat([df_control, df_case], axis=0, ignore_index=True).sort_values(by=['RID', 'years'])

        return df_result

    NC_volume_source = f"./Result/volume_NC_{hemisphere_status}.csv"
    AD_volume_source = f"./Result/volume_{group}_{hemisphere_status}.csv"
    if val:
        AD_val_volume_source = f"./Result/volume_{group}_val_{hemisphere_status}.csv"

    df_ad = pd.read_csv(AD_volume_source)
    df_ad["group"] = 1
    df_nc = pd.read_csv(NC_volume_source)
    df_nc["group"] = 0

    if val:
        df_ad_val = pd.read_csv(AD_val_volume_source)
        df_ad_val["group"] = 2
        df_ad = pd.concat([df_ad, df_ad_val], axis=0, ignore_index=True)
        

    df_result = pd.concat([df_ad, df_nc], axis=0, ignore_index=True)
    df_result = df_result.sort_values(by=['RID', 'years'])

    df_cdr = pd.read_csv('./data/OASIS3_UDSb4_cdr.csv')
    df_cdr['years'] = df_cdr['OASIS_session_label'].apply(get_years)
    df_cdr['RID'] = df_cdr["OASISID"].apply(strip_prefix)

    df_cdr_first = df_cdr[['RID', 'dx1_code', 'CDRTOT', 'years']].groupby('RID').first().reset_index()

    subject_dxbl_nan = df_cdr_first[(df_cdr_first['years'] >= 0.5) | (df_cdr_first['years'] <= -0.5)]['RID'].unique().tolist()
    df_cdr_first['diagnosis_bl'] = df_cdr_first['CDRTOT'].apply(cdr2dia)
    df_cdr_first.loc[df_cdr_first['RID'].isin(subject_dxbl_nan), 'diagnosis_bl'] = np.nan

    df_result = pd.merge(df_result, df_cdr_first[['RID', 'diagnosis_bl']], on='RID')

    df_result_RID_list = []
    for RID in df_result['RID'].unique().tolist():
        df_result_RID = df_result[df_result['RID'] == RID]
        df_cdr_RID = df_cdr[df_cdr['RID'] == RID][['RID', 'years', 'CDRTOT', 'MMSE']]
        df_result_RID = pd.merge_asof(df_result_RID, df_cdr_RID, on="years", by='RID', tolerance=0.5, direction='nearest')
        df_result_RID_list.append(df_result_RID)
    df_result = pd.concat(df_result_RID_list)

    df_result['diagnosis'] = df_result['CDRTOT'].apply(cdr2dia)

    df_A1 = pd.read_csv("./data/OASIS3_UDSa1_participant_demo.csv")
    df_A1["OASISID"] = df_A1["OASISID"].apply(strip_prefix)
    df_demo = pd.read_csv("./data/OASIS3_demographics.csv")
    df_demo["OASISID"] = df_demo["OASISID"].apply(strip_prefix)

    df_gender = df_demo[["OASISID", "GENDER"]]
    d = {key: value for key, value in zip(df_gender['OASISID'], df_gender["GENDER"] % 2)}
    df_result["PTGENDER"] = df_result["RID"].map(d)

    df_age = df_demo[["OASISID", "AgeatEntry"]]
    d = {key: value for key, value in zip(df_age['OASISID'], df_age["AgeatEntry"])}
    df_result["AGE_baseline"] = df_result["RID"].map(d)
    df_result["AGE"] = df_result["AGE_baseline"] + df_result["years"]

    df_educ = df_demo[["OASISID", "EDUC"]]
    d = {key: value for key, value in zip(df_educ['OASISID'], df_educ["EDUC"])}
    df_result["PTEDUCAT"] = df_result["RID"].map(d)

    df_apoe = df_demo[["OASISID", "APOE"]]
    df_apoe["APOE"] = df_apoe["APOE"].apply(count_number_4)
    d = {key: value for key, value in zip(df_apoe['OASISID'], df_apoe["APOE"])}
    df_result["APOE"] = df_result["RID"].map(d)

    df_result = df_result[
        ['RID', 'years', 'diagnosis', 'AGE_baseline', 'AGE', 'PTGENDER', 'PTEDUCAT', 'APOE', 'group',
         'CDRTOT', 'MMSE', 'diagnosis_bl', 'FLDSTRENG'] + biomarker_ver_map[hemisphere_status]] 

    # if isextra:
    #     df_abeta = pd.read_csv('./abeta/abeta_all.csv')
    #     df_tau = pd.read_csv('./tau/tau_all.csv')
    #     if hemisphere_status in ['HM', 'HS']:
    #         df_result_RID_list = []
    #         for RID in df_result['RID'].unique().tolist():
    #             df_result_RID = df_result[df_result['RID'] == RID]
    #             df_abeta_RID = df_abeta[df_abeta['RID'] == RID][
    #                 ['RID', 'years'] + biomarker_ver_map[hemisphere_status][1:]]
    #             df_abeta_RID.rename(
    #                 columns={item: 'Abeta_' + item for item in biomarker_ver_map[hemisphere_status][1:]}, inplace=True)
    #             df_tau_RID = df_tau[df_tau['RID'] == RID][
    #                 ['RID', 'years'] + biomarker_ver_map[hemisphere_status][1:]]
    #             df_tau_RID.rename(
    #                 columns={item: 'Tau_' + item for item in biomarker_ver_map[hemisphere_status][1:]}, inplace=True)
    #             df_result_RID = pd.merge_asof(df_result_RID, df_abeta_RID, on="years", by='RID', tolerance=0.5, direction='nearest')
    #             df_result_RID = pd.merge_asof(df_result_RID, df_tau_RID, on="years", by='RID', tolerance=0.5, direction='nearest')
    #             df_result_RID_list.append(df_result_RID)
    #         df_result = pd.concat(df_result_RID_list)

    #     elif hemisphere_status in ['LM', 'LS']:
    #         high_hemisphere_status = hemisphere_status.replace('L', 'H')
    #         if val:
    #             df_high = pd.read_csv(f'./output/OASIS3_input_{group}_{high_hemisphere_status}.csv')
    #         else:
    #             df_high = pd.read_csv(f'./output/OASIS3_input_{group}_{high_hemisphere_status}.csv')
    #         df_abeta_low = df_high[['RID', 'years']]
    #         df_tau_low = df_high[['RID', 'years']]
    #         for region in biomarker_ver_map[hemisphere_status][1:]:
    #             source_col = eval(f'map_{hemisphere_status}_{high_hemisphere_status}[region]')
    #             weight_lst = []
    #             item_abeta_lst = []
    #             item_tau_lst = []
    #             for col in source_col:
    #                 weight_lst.append(df_high[col])
    #                 item_abeta_lst.append(df_high['Abeta_' + col])
    #                 item_tau_lst.append(df_high['Tau_' + col])
    #             df_abeta_low['Abeta_' + region] = weighted_sum(weight_lst, item_abeta_lst)
    #             df_tau_low['Tau_' + region] = weighted_sum(weight_lst, item_tau_lst)

    #         df_result = pd.merge(df_result, df_abeta_low, how='left', on=['RID', 'years'])
    #         df_result = pd.merge(df_result, df_tau_low, how='left', on=['RID', 'years'])
    #         df_result = df_result[['RID', 'years', 'diagnosis', 'AGE_baseline', 'AGE', 'PTGENDER', 'PTEDUCAT', 'APOE', 'group',
    #      'CDRTOT', 'MMSE', 'diagnosis_bl', 'FLDSTRENG'] +
    #                               biomarker_ver_map[hemisphere_status] +
    #                               ['Abeta_' + item for item in biomarker_ver_map[hemisphere_status][1:]] +
    #                               ['Tau_' + item for item in biomarker_ver_map[hemisphere_status][1:]]]
    if match:
        df_result = match_demo(df_result, n_matches=1, with_replacement=False)
    else:
        # df_result = df_result[(df_result['AGE_baseline'] >= 69) | (df_result['group'] == 1)]
        df_result = df_result[(df_result['AGE_baseline'] >= 70) | (df_result['group'] != 0)]
    df_result.sort_values(by=["RID", "years"], inplace=True, ignore_index=True)
    df_result.dropna(subset=df_result.columns[4:8], inplace=True, ignore_index=True)
    # if isextra:
    #     df_result = delete_duplicate_row(df_result, ['Abeta_' + item for item in biomarker_ver_map[hemisphere_status][1:]])
    #     df_result = delete_duplicate_row(df_result, ['Tau_' + item for item in biomarker_ver_map[hemisphere_status][1:]])

    print('CN Group')
    print("# of subjects: ", len(df_result[df_result['group'] == 0]['RID'].unique().tolist()))
    print("# of visits: ", df_result[df_result['group'] == 0].shape[0])
    print("AGE: mean: ", np.mean(df_result[df_result['group'] == 0]['AGE']), ", std: ", np.std(df_result[df_result['group'] == 0]['AGE']))

    print('AD Group')
    print("# of subjects: ", len(df_result[df_result['group'] == 1]['RID'].unique().tolist()))
    print("# of visits: ", df_result[df_result['group'] == 1].shape[0])
    print("AGE: mean: ", np.mean(df_result[df_result['group'] == 1]['AGE']), ", std: ", np.std(df_result[df_result['group'] == 1]['AGE']))

    from scipy import stats
    t_statistic, p_value = stats.ttest_ind(df_result[df_result['group'] == 0]['AGE'], df_result[df_result['group'] == 1]['AGE'])
    print(f"T-statistic: {t_statistic}")
    print(f"P-value: {p_value}")
    df_result.to_csv(output_path, index=False)

def choose_nc_subject() -> List:
    
    df = pd.read_csv("./data/OASIS3_UDSb4_cdr.csv")
    df_nan_removed = df.dropna(subset=['dx1_code'])
    df_nc = delete_subject_discrete(df_nan_removed, "dx1_code", [1])
    subjects_nc = df_nc["OASISID"].unique().tolist()

    df_centiloid = pd.read_csv("./data/OASIS3_amyloid_centiloid.csv")
    subject_no_abeta = set(subjects_nc) - set(df_centiloid["subject_id"].unique())

    df_centiloid = df_centiloid[df_centiloid["subject_id"].isin(subjects_nc)]
    df_centiloid_av45 = df_centiloid[df_centiloid["tracer"] == "AV45"]
    df_centiloid_av45 = df_centiloid_av45[df_centiloid_av45["Centiloid_fSUVR_TOT_CORTMEAN"] >= 34.9]
    df_centiloid_pib = df_centiloid[df_centiloid["tracer"] == "PIB"]
    df_centiloid_pib = df_centiloid_pib[df_centiloid_pib["Centiloid_fSUVR_TOT_CORTMEAN"] >= 15.3]

    subject_abeta_positive = set(df_centiloid_av45["subject_id"].unique()).union(
        set(df_centiloid_pib["subject_id"].unique()))

    subject_nc_result = set(subjects_nc) - subject_abeta_positive - subject_no_abeta
    return subject_nc_result

def choose_rod1_subject() -> List:
    df = pd.read_csv("./data/OASIS3_UDSb4_cdr.csv")
    df_ad = choose_subject_discrete(df, "dx1_code", [3])
    subjects_ad = df_ad["OASISID"].unique().tolist()

    df_centiloid = pd.read_csv("./data/OASIS3_amyloid_centiloid.csv")
    df_centiloid = df_centiloid[df_centiloid["subject_id"].isin(subjects_ad)]
    df_centiloid_av45 = df_centiloid[df_centiloid["tracer"] == "AV45"]
    df_centiloid_av45 = df_centiloid_av45[df_centiloid_av45["Centiloid_fSUVR_TOT_CORTMEAN"] >= 34.9]
    df_centiloid_pib = df_centiloid[df_centiloid["tracer"] == "PIB"]
    df_centiloid_pib = df_centiloid_pib[df_centiloid_pib["Centiloid_fSUVR_TOT_CORTMEAN"] >= 15.3]

    subects_abeta_positive = list(set(df_centiloid_av45["subject_id"].unique()).union(
        set(df_centiloid_pib["subject_id"].unique())))

    subects_ad_result = list(set(subjects_ad).intersection(subects_abeta_positive))

    # print("ROD1 group # of subjects:", len(subects_ad_result))
    return subects_ad_result

def choose_rod1_val_subject() -> List:
    """
    选出诊断为AD (dx1_code=3)，但在OASIS3_amyloid_centiloid.csv中
    没有AV45或PIB示踪剂记录的受试者，且至少有一次CDR评分 >= 1。
    """
    df_dx = pd.read_csv("./data/OASIS3_UDSb4_cdr.csv")
    df_ad = choose_subject_discrete(df_dx, "dx1_code", [3])
    subjects_ad = df_ad["OASISID"].unique().tolist()
    subjects_ad_set = set(subjects_ad)

    df_ad_with_cdr = df_ad[df_ad["CDRTOT"] >= 1]
    subjects_ad_cdr = df_ad_with_cdr["OASISID"].unique().tolist()
    subjects_ad_cdr_set = set(subjects_ad_cdr)
    subjects_ad_set = subjects_ad_set.intersection(subjects_ad_cdr_set)

    df_centiloid = pd.read_csv("./data/OASIS3_amyloid_centiloid.csv")

    df_centiloid_relevant = df_centiloid[df_centiloid["subject_id"].isin(subjects_ad_set)]

    subjects_with_av45_or_pib_data = df_centiloid_relevant[
        df_centiloid_relevant["tracer"].isin(["AV45", "PIB"])
    ]["subject_id"].unique()  # 有一个指标即可选入
    subjects_with_data_set = set(subjects_with_av45_or_pib_data)

    subjects_without_data = list(subjects_ad_set - subjects_with_data_set)

    print(f"AD subjects: {len(subjects_ad)}")
    print(f"AD subjects with CDR >= 1: {len(subjects_ad_set)}")
    print(f"AD subjects with CDR >= 1 and AV45 or PIB data: {len(subjects_with_data_set)}")
    print(f"AD subjects with CDR >= 1 and without AV45 or PIB data: {len(subjects_without_data)}")

    return subjects_without_data

def main():
    os.makedirs('Result', exist_ok=True)

    subjects_ad = choose_rod1_subject()

    subjects_nc = choose_nc_subject()
    
    subjects_ad_val = choose_rod1_val_subject()

    subjects_data = {"NC": subjects_nc, "ROD1": subjects_ad, "ROD1_val": subjects_ad_val}

    for group in ["NC", "ROD1", "ROD1_val"]:
        get_volume_data(subjects_data[group], group)

    # get_abeta_form()

    for hemisphere_status in version:
        for group in ["ROD1"]:
            get_input_form(f"./Result/OASIS3_input_{group}_{hemisphere_status}.csv", group, hemisphere_status, isextra=False, val=True)


if __name__ == '__main__':
    import warnings
    warnings.filterwarnings("ignore")
    main()
