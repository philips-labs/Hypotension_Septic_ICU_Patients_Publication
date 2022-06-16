# module fetch_noninvasive_bp_data.py

import os
import sys

import pandas as pd
import numpy as np
from collections import Counter

from load_pyspark import load_spark, assert_pyspark
spark = load_spark()
F, window = assert_pyspark()
print(spark.sparkContext.uiWebUrl)

from helper_scripts import to_parquet
from pyspark.sql.functions import max as _max, min as _min, count

######################################################################################################
def fetch_noninvasive_all_readings(eicu_path, eicu_path2):

	patient_id_bins = np.arange(0, 4000001, 100000)

	patients_df = pd.read_csv(os.path.join(eicu_path2, 'patients_filters_first_pass.csv'), index_col=0)
	vitals_aperiodic_df = spark.read.load(os.path.join(eicu_path, 'vitalAperiodic.parquet'))
	vitals_aperiodic_df = vitals_aperiodic_df.select([F.col('observationOffset').cast('int'), F.col('patientUnitStayID').cast('int'),
							F.col('nonInvasiveSystolic').cast('int'), F.col('nonInvasiveDiastolic').cast('int'),
							F.col('nonInvasiveMean').cast('int'), 'patientHealthSystemStayID', 'hospitalAdmitOffset',
							'hospitalID'])

	for p in np.arange(1, len(patient_id_bins)):
		print(patient_id_bins[p-1], patient_id_bins[p])
		target_pid = list(patients_df.loc[np.logical_and(patients_df['patientUnitStayID'] >= patient_id_bins[p-1],\
								patients_df['patientUnitStayID'] < patient_id_bins[p]), 'patientUnitStayID'])
		print('patient counts=', len(target_pid))
		if len(target_pid) > 0:
			temp_df = vitals_aperiodic_df.filter(vitals_aperiodic_df.patientUnitStayID.isin(target_pid))
			temp_df.write.mode('overwrite').parquet(os.path.join(eicu_path2, 'noninvasive_BP_readings_batch_%s_%s.parquet' % (patient_id_bins[p-1], patient_id_bins[p])))

######################################################################################################
if __name__ == '__main__':

	eicu_path = '/some/path/dataset'
	eicu_path2 = 'some/path/dataset'
	
	fetch_noninvasive_all_readings(eicu_path, eicu_path2)

