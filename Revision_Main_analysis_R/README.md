# Revision main analyses in R
This folder contains R scripts for the analyses.<br>
We have updated the analysis based on the revision request from Annals of Intensive Care.<br>

1. Thin plate logistic regression spline and visualization of predicted probabilities with histogram of exposures for primary outcome: `Revision_Main_Analysis.R`

2. Thin plate logistic regression spline and visualization of predicted probabilities with histogram of exposures for secondary outcomes: `Revision_Secondary_Analysis.R`

3. Categorizing exposure levels to different levels to explore relative risks and subgroup analyses: `Revision_Sepsis_Categorical_Analysis.R`<br>
        <br>
        **Categorizing exposure levels**
        <p>
        - MAP
        <ol type="a">
        <li>Mean BP ≥75mmHg (reference)</li>
        <li>Mean BP <75mmHg</li>
        <li>Mean BP <65mmHg</li>
        <li>Mean BP <55mmHg</li>
        <li>Mean BP <45mmHg</li>
        </ol>
                A vs B, A vs C, A vs D, A vs E
        </p>
        <p>
        - SBP
        <ol type="a">
        <li>sBP ≥110mmHg (reference)</li>
        <li>sBP <110mmHg</li>
        <li>sBP <100mmHg</li>
        <li>sBP <90mmHg</li>
        <li>sBP <80mmHg</li>
        </ol>
                A vs B, A vs C, A vs D, A vs E
        </p>
        <p>
        - DBP
        <ol type="a">
        <li>dBP ≥60mmHg (reference)</li>
        <li>dBP <60mmHg</li>
        <li>dBP <50mmHg</li>
        <li>dBP <40mmHg</li>
        <li>dBP <30mmHg</li>
        </ol>
                A vs B, A vs C, A vs D, A vs E
        </p>
        <p>
        - PP
        <ol type="a">
        <li>40 ≤ PP < 50mmHg (reference)</li>
        <li>PP ≥50mmHg</li>
        <li>PP <40mmHg</li>
        <li>PP <30mmHg</li>
        <li>PP <20mmHg</li>
        </ol>
                A vs B, A vs C, A vs D, A vs E
        </p>
        **Subgroup analysis**
        <p>
        - Age (year)
        <ol type="a">
        <li><45 years</li>
        <li>45 to 64 years</li>
        <li>65 years or more</li>
        </ol>
        </p>
        <p>
        - Vasopressor usage (in NEE mcg/kg/min)
        <ol type="a">
        <li>0</li>
        <li>0 to 0.05</li>
        <li>0.05 to 0.1</li>
        <li>0.1 or more</li>
        </ol>
        </p>
        <p>
        - On ventilator
        <ol type="a">
        <li>No</li>
        <li>Yes</li>
        </ol>
        </p>
        <p>
        - Cancer/Tumor
        <ol type="a">
        <li>No</li>
        <li>Yes</li>
        </ol>
        </p>
        <p>
        - Admission type
        <ol type="a">
        <li>Emergency department</li>
        <li>Other ward</li>
        <li>Elective</li>
        <li>Other hospital</li>
        </ol>
        </p>

4. Threshold logistic regression analysis: `Revision_Threshold_Analysis.R`

5. Thin plate logistic regression spline and visualization of predicted probabilities using time-weighted average under the threshold defined by the threshold logistic regression: `Revision_Time_Weighted_Average.R`
