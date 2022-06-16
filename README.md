# ICU Hypotension Project

This repository contains scripts for the intensive care unit (ICU) hypotension project. This is a collaborative work with Dr. Ashish K. Khanna at Wake Forest School of Medicine and Dr. Kamal Maheshwari at Cleveland Clinic. The study design meeting was started in October 2021 and a *priori* defined analysis plan was approved by Wake Forest University School of Medicine Institutional Review Board in Feburary 2022. We used the eICU Research Institute (eRI) dataset. Takahiro Kinoshita, Annamalai Natarajan, Junzi Dong, and Emma Schwager are collaborators.

# Association of Systolic, Diastolic, Mean, and Pulse Pressure with Morbidity and Mortality in Septic ICU Patients

The first study focused on the septic ICU patients. We included adult ICU patients (aged ≥ 18 years) with non-surgical sepsis diagnosis from the eRI database. We first identified adult patients with infection from either admission diagnoses or records of non-prophylactic antibiotics, defined as the antibiotic use for lasting more than 48 hours. Patients who underwent surgery before ICU admission were not included as intraoperative blood pressure records were not available. Sequential Organ Failure Assessment (SOFA) scores were collected from bedside monitor data, laboratory measurements, and medication data. Then, we identified sepsis patients as those with infection and overall SOFA score ≥ 2 on the day of admission. 

We excluded patients without APACHE IVa score, who had invalid body mass index (BMI) (≤ 10 kg/m2, ≥ 60 kg/m2, or missing height or weight), had a do not resuscitate order on admission, received mechanical circulatory support, were discharged from ICU within 24 hours or stayed in ICU longer than 14 days, or had any gaps between blood pressure readings of greater than 2 hours. Patients who had acute kidney injury (AKI) or myocardial injury (secondary outcomes) within 24 hours after ICU admission were also excluded to confirm the temporal relationship between ICU hypotension and organ system dysfunction and to remove the possibility of reverse causation. We excluded patients with a history of organ dysfunction such as chronic kidney disease, myocardial infarction, stroke, and coronary artery disease. Patients with estimated glomerular filtration rate (eGFR) less than 60 ml/min/1.73 m2 on admission based on the Cockcroft–Gault equation were also excluded.

1. For cohort selection please see Cohort_selection_python
2. For the main analysis plaese see Main_analysis_R 
