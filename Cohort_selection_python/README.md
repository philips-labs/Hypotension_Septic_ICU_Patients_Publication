# Cohort selection:

This section shows how we extracted data from the eRI dataset.
As the eRI dataset is not publicly available, the external users cannot run the following code.<br>
We provide the subset of the dataset, which consists of patients who are included in the publicly available eICU dataset, for the analysis.
See Main_analysis_R folder.

1. First set of exclusion criteria: `python cohort_selection.ipynb` -> `sepsis_patients_filters_first_pass.csv` -> # 137,197 patients
    1. Total ICU stays with sepsis: 301,447
    2. Apache IVa score missing: 73,032
    3. Invalid BMI: 12,762
    4. DNR: 3,317
    5. Mechanical cardiac support (IABP, LVAD, RVAD, BVAD, Impella, and ECMO): 4,149
    6. Length of stay < 1 day or > 14 days: 10,905
    7. Myocardial injury within 24 hours: 29,066
    8. AKI within 24 hours: 31,199

2. Fetch invasive BP readings: `python fetch_invasive_a_line_data.py` -> `invasive_BP_readings_batch_*.parquet`

3. Fetch non-invasive BP readings: `python fetch_noninvasive_bp_data.py` -> `noninvasive_BP_readings_batch_*.parquet`

4. Merge invasive and non-invasive: `python merge_invasive_noninvasive_vitals.py` -> `invasive_BP_selectively_merged_batch_*.parquet`
    1. Apply plausability filters per BP source per modality i.e. invasive or non-invasive respectively
    2. Computes MBP if missing from SBP and DBP per modality. Note we do not mix BP data between modalities
    3. Retaining only BP triplets (i.e. all three readings will need to present within each modality)
    5. Fill up gaps in invasive BP readings with non-invasive readings when gap > 5 minutes. Gaps include at the beginning (time 0 to first invasive BP reading) and at the end (last invasive BP reading to time of outcome)

5. Second set of exclusion criteria based on merged BP readings: `python cohort_selection.ipynb` -> `sepsis_patients_filters_second_pass.csv` -> # 52,261 patients
    1. Computes PP from SBP and DBP
    2. Computes max time gap between adjacent pairs of BP readings per patient
    3. Drops patients with a single gap of > 2 hours anytime during ICU stay (beginning, in-between, end)
    4. Computes minimum SBP, DBP, MBP, PP reading for each outcome that is cumulatively but non-continuously <= 2 hours (e.g., for patient N the SBP threshold of 96 meets the >= 2 hour criteria) 

6. Third set of exclusion criteria: `python cohort_selection.ipynb` -> `sepsis_patients_filters_third_pass.csv` -> # 22,494 patients
    1. Past medical history of organ dysfunction: 18,249
    2. eGFR < 60 ml/min/1.73 m2: 11,518
