# Main analyses in R

This folder contains R scripts for the analyses.<br>
Although the entire cohort is not published, credentialed Physionet users who have finished required training can use the subset of the dataset, which consists of patients included in the publicly available eICU dataset, for the analysis.<br>
The original cohort has 22,494 patients, and 2,284 of them are available in the subcohort.<br>
Please contact Takahiro Kiritoshi <takahiro.kiritoshi@philips.com> to request the data.

1. Thin plate logistic regression spline and visualization of predicted probabilities with histogram of exposures: `Main_Analysis.R`

2. Categorizing exposure levels to different levels to explore relative risks and subgroup analyses: `Sepsis_Categorical_Analysis.R`
	<p><ol type="a">
		<li>Mean BP ≥75mmHg (reference)</li>
		<li>Mean BP <75mmHg</li>
		<li>Mean BP <65mmHg</li>
		<li>Mean BP <55mmHg</li>
		<li>Mean BP <45mmHg</li>
	</ol>
	a vs b, a vs c, a vs d, a vs e</p>
	
	<p><ol type="a">
		<li>sBP ≥110mmHg (reference)</li>
		<li>sBP <110mmHg</li>
		<li>sBP <100mmHg</li>
		<li>sBP <90mmHg</li>
		<li>sBP <80mmHg</li>
	</ol>
	a vs b, a vs c, a vs d, a vs e</p>
	
	<p><ol type="a">
		<li>dBP ≥60mmHg (reference)</li>
		<li>dBP <60mmHg</li>
		<li>dBP <50mmHg</li>
		<li>dBP <40mmHg</li>
		<li>dBP <30mmHg</li>
	</ol>
	a vs b, a vs c, a vs d, a vs e</p>
	
	<p><ol type="a">
		<li>30 ≤ PP < 40mmHg (reference)</li>
		<li>PP ≥50mmHg</li>
		<li>PP ≥40mmHg</li>
		<li>PP <30mmHg</li>
		<li>PP <20mmHg</li>
	</ol>
	a vs b, a vs c, a vs d, a vs e</p>
	
	
	subgroup analysis
	- Age (year)<br>
		<45 years<br>
		45 to 64 years<br>
		65 years or more
	- Vasopressor usage (in NEE mcg/kg/min)<br>
		0<br>
		0 to 0.05<br>
		0.05 to 0.1<br>
		0.1 or more
	- Septic shock<br>
		No<br>
		Yes
	- On ventilator<br>
		No<br>
		Yes
	- Cancer/Tumor<br>
		No<br>
		Yes
	- Admission type<br>
		Emergency department<br>
		Other ward<br>
		Elective<br>
		Other hospital

3. Threshold logistic regression analysis: `Threshold_Analysis.R`
