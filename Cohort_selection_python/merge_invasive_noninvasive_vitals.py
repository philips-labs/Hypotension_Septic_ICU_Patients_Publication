# module merge_invasive_noninvasive_vitals.py

import os
import sys

import pandas as pd
import numpy as np
from collections import Counter

from helper_scripts import to_parquet
from pyspark.sql.functions import max as _max, min as _min, count

######################################################################################################
def filter_bp(bp_df, bp_cols, invasive_non_invasive_flag, plausability_filters):

	for idx, bp_vital in enumerate(bp_cols):
		# Get the plausability min/max for each vital
		plaus_min_max = plausability_filters[bp_vital]

		# Plausability filters
		# Exclude readings that have,
		# 1. SBP greater than or equal to 300 or less than or equal to 20 mmHg
		# 2. SBP less than or equal to DBP + 5 mmHg
		# 3. DBP less than or equal to 5 mmHg or greater than or equal to 225 mmHg
		# 4. MBP greater than or equal to SBP or less than or equal to DBP
		# 5. MBP range derived from SBP and DBP range  

		# Replace with NaN and drop non-invasive vital outside of range
		target_idx = np.logical_or(bp_df[bp_vital] < plaus_min_max[0], bp_df[bp_vital] > plaus_min_max[1])
		bp_df.loc[target_idx, bp_vital] = np.NaN

		# Replace with NaN and drop vitals not meeting criterion
		if 'Systolic' in bp_vital:
			if 'non_invasive' in invasive_non_invasive_flag:
				target_idx = bp_df['nonInvasiveSystolic'] < bp_df['nonInvasiveDiastolic'] + 5
			else:
				target_idx = bp_df['systemicSystolic'] < bp_df['systemicDiastolic'] + 5
		elif 'Mean' in bp_vital:
			if 'non_invasive' in invasive_non_invasive_flag:
				target_idx = np.logical_or(bp_df['nonInvasiveMean'] >= bp_df['nonInvasiveSystolic'],\
							bp_df['nonInvasiveMean'] <= bp_df['nonInvasiveDiastolic'])
			else:
				target_idx = np.logical_or(bp_df['systemicMean'] >= bp_df['systemicSystolic'],\
							bp_df['systemicMean'] <= bp_df['systemicDiastolic'])

		bp_df.loc[target_idx, bp_vital] = np.NaN

	if 'non_invasive' in invasive_non_invasive_flag:
		mbp_nan_idx = bp_df['nonInvasiveMean'].isnull()
		bp_df.loc[mbp_nan_idx, 'nonInvasiveMean'] = bp_df.loc[mbp_nan_idx, 'nonInvasiveDiastolic'] +\
								(1/3 * (bp_df.loc[mbp_nan_idx, 'nonInvasiveSystolic'] - bp_df.loc[mbp_nan_idx, 'nonInvasiveDiastolic']))
	else:
		mbp_nan_idx = bp_df['systemicMean'].isnull()
		bp_df.loc[mbp_nan_idx, 'systemicMean'] = bp_df.loc[mbp_nan_idx, 'systemicDiastolic'] +\
								(1/3 * (bp_df.loc[mbp_nan_idx, 'systemicSystolic'] - bp_df.loc[mbp_nan_idx, 'systemicDiastolic']))

	# Dropping rows with one or more NaN's in the triplets
	bp_df = bp_df.loc[bp_df[bp_cols].notnull().all(axis=1), :]

	return (bp_df)

######################################################################################################
def selective_merge_invasive_noninvasive_vitals(eicu_path, eicu_path2):

	patient_id_bins = np.arange(0, 4000001, 100000)	

	invasive_cols = ['systemicSystolic', 'systemicDiastolic', 'systemicMean']
	noninvasive_cols = ['nonInvasiveSystolic', 'nonInvasiveDiastolic', 'nonInvasiveMean']
	merged_cols = ['mergedSystolic', 'mergedDiastolic', 'mergedMean']
	events_of_interest = ['MI_offset', 'AKI_offset', 'Death_offset', 'Discharge_offset']

	patient_observation_offset_cols = ['patientUnitStayID', 'observationOffset']
	plausability_filters = {'systemicSystolic': [20, 300], 'systemicDiastolic': [5, 225], 'systemicMean': [10, 250],\
				'nonInvasiveSystolic': [20, 300], 'nonInvasiveDiastolic': [5, 225], 'nonInvasiveMean': [10, 250]}
	how_much_gap = 5 # minutes

	patients_df = pd.read_csv(os.path.join(eicu_path2, 'patients_filters_first_pass.csv'), index_col=0)

	for p in np.arange(1, len(patient_id_bins)):
		print(p)
		# Check if invasive and non-invasive vitals exist for partition p
		invasive_df_filename = os.path.join(eicu_path2, 'invasive_BP_readings_batch_%s_%s.parquet' % (patient_id_bins[p-1], patient_id_bins[p]))
		noninvasive_df_filename = os.path.join(eicu_path2, 'noninvasive_BP_readings_batch_%s_%s.parquet' % (patient_id_bins[p-1], patient_id_bins[p]))

		if os.path.exists(invasive_df_filename) and os.path.exists(noninvasive_df_filename):
			print(invasive_df_filename, noninvasive_df_filename)

			# Load invasive and non-invasive vitals for partition p
			full_invasive_df = pd.read_parquet(invasive_df_filename)
			full_noninvasive_df = pd.read_parquet(noninvasive_df_filename)

			if len(full_invasive_df) > 0 and len(full_noninvasive_df) > 0:

				full_invasive_df = filter_bp(full_invasive_df, invasive_cols, 'invasive', plausability_filters)
				full_noninvasive_df = filter_bp(full_noninvasive_df, noninvasive_cols, 'non_invasive', plausability_filters)

				# Houskeeping I
				merged_df = pd.DataFrame()

				# Loop over each vital
				for idx, (inv_vitals, non_vitals) in enumerate(zip(invasive_cols, noninvasive_cols)):
					invasive_df = full_invasive_df.copy()
					noninvasive_df = full_noninvasive_df.copy()

					# Columns of interest
					cols_to_retain = patient_observation_offset_cols + [inv_vitals]
	
					selected_patients_idx = np.logical_and(patients_df['patientUnitStayID'] >= patient_id_bins[p-1],\
										patients_df['patientUnitStayID'] < patient_id_bins[p])

					# Create a df to collect time 0 to first available vital time. Note the patients are identified from patients_df NOT invasive/non-invasive df
					first_df = pd.DataFrame()
					first_df['patientUnitStayID'] = patients_df.loc[selected_patients_idx, 'patientUnitStayID'].copy()
					first_df['observationOffset'] = 0
					first_df[inv_vitals] = -999.888

					# Create a df to collect time from last available vital to time of discharge/death. Note the patients are identified from patients_df NOT invasive/non-invasive df
					last_df = pd.DataFrame()
					last_df['patientUnitStayID'] = patients_df.loc[selected_patients_idx, 'patientUnitStayID'].copy()
					last_df['observationOffset'] = patients_df.loc[selected_patients_idx, events_of_interest].max(axis=1)
					last_df[inv_vitals] = -999.888

					# Merging first, last and invasive df
					temp_df = pd.concat([invasive_df.loc[:, cols_to_retain], first_df, last_df], ignore_index=True)
					# Sorting by patients and by observation offset
					temp_df = temp_df.sort_values(by=patient_observation_offset_cols)
					temp_df.reset_index(drop=True, inplace=True)
					# Getting time difference between adjacent timestamps within each patient
					temp_df['time_diff'] = temp_df.groupby('patientUnitStayID')['observationOffset'].diff().fillna(-999.888)

					# Looking for rows in which time difference > 5
					target_idx2, = np.where(temp_df['time_diff'] > how_much_gap)
					# Getting the start and end of the gap of invasive vitals in minutes
					gap_df = temp_df.loc[target_idx2, patient_observation_offset_cols]
					gap_df.columns = ['patientUnitStayID', 'observationOffset_end']
					gap_df['observationOffset_start'] = list(temp_df.loc[target_idx2-1, 'observationOffset'])

					# Merging this gap in invasive vitals with non-invasive vitals
					selected_noninvasive_vitals_df = gap_df.merge(noninvasive_df, on='patientUnitStayID', how='inner')
					# Keeping only rows in noninvasive vitals that are inside this gap
					target_idx = np.logical_and(selected_noninvasive_vitals_df['observationOffset'] > selected_noninvasive_vitals_df['observationOffset_start'],\
								selected_noninvasive_vitals_df['observationOffset'] < selected_noninvasive_vitals_df['observationOffset_end'])
					selected_noninvasive_vitals_df = selected_noninvasive_vitals_df.loc[target_idx, :]

					# Concatenating those chosen non-invasive vitals with invasive vitals
					temp_df = pd.concat([temp_df[cols_to_retain], selected_noninvasive_vitals_df[patient_observation_offset_cols + [non_vitals]]], ignore_index=True)

					# Outer join with other vitals
					if idx == 0:
						merged_df = temp_df.copy()
					else:
						merged_df = merged_df.merge(temp_df, on=patient_observation_offset_cols, how='outer')

				merged_df.replace(-999.888, np.NaN, inplace=True)

				print('Before=', merged_df[invasive_cols].isnull().sum())
				for inv_vitals, non_vitals, merged_vitals in zip(invasive_cols, noninvasive_cols, merged_cols):
					merged_df[merged_vitals] = np.NaN
					nan_idx = merged_df[inv_vitals].isnull()
					merged_df.loc[~nan_idx, merged_vitals] = merged_df.loc[~nan_idx, inv_vitals]
					merged_df.loc[nan_idx, merged_vitals] = merged_df.loc[nan_idx, non_vitals]
				print('After=', merged_df[merged_cols].isnull().sum())

				to_parquet(merged_df, os.path.join(eicu_path2, 'invasive_BP_selectively_merged_batch_%s_%s.parquet' % (patient_id_bins[p-1], patient_id_bins[p])))
		else:
			print('Files do not exist')
			print(os.path.exists(invasive_df_filename))
			print(os.path.exists(noninvasive_df_filename))

######################################################################################################
if __name__ == '__main__':

	eicu_path = 'some/path/dataset'
	eicu_path2 = 'some/path/dataset'
	
	selective_merge_invasive_noninvasive_vitals(eicu_path, eicu_path2)

