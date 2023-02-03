#!/usr/bin/env bash

# Ensure the joins work
sqlite3 /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'delete from KEY_METRICS where Sample = "NA";'


# Dump all of the tables as CSV's
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'select * from File_Locations;' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/File_Locations.csv
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'select * from KEY_METRICS;' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/KEY_METRICS.csv
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'select * from status;' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/status.csv
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'select * from Samples;' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/Samples.csv
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'select * from Timestamps;' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/Timestamps.csv

# Early report for gauging topups and sample quality
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'SELECT Samples.Sample_True as Patient,
norm.WGS_Dedup_Coverage as Normal_Deduplex_coverage, 
tum.WGS_Dedup_Coverage as Tumour_Deduplex_coverage,
rna.Raw_Reads_Count as Raw_Reads_Count
FROM Samples 
INNER JOIN KEY_METRICS norm 
ON Samples.DNA_N = norm.Sample
INNER JOIN KEY_METRICS tum 
ON Samples.DNA_T = tum.Sample
INNER JOIN KEY_METRICS rna 
ON Samples.RNA = rna.Sample;' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/Initial_metrics.csv

# Report to give all of the metrics related to each Sample.
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'SELECT Samples.Sample_True AS Patient,
Samples.DNA_N_True AS WGS_Normal,
Samples.DNA_T_True AS WGS_Tumour,
status.dna_pipeline_execution as WGS_Complete,
Samples.RNA_True AS WTS,
status.rna_pipeline_light_execution as WTSLight_Complete,
CASE WHEN (Samples.DNA_N_True == "NA") THEN "NA" WHEN (tum.WGS_Dedup_Coverage < 80) THEN "FAIL" WHEN (norm.WGS_Dedup_Coverage < 40) THEN "FAIL" ELSE "PASS" END WGS_Coverage_Pass,
CASE WHEN (Samples.RNA_True == "NA") THEN "NA" WHEN (rna.Raw_Reads_Count < 100000000) THEN "FAIL"  ELSE "PASS" END Raw_Reads_Count_Pass,
CASE WHEN (Samples.DNA_N_True == "NA") THEN "NA" WHEN (Samples.RNA_True == "NA") THEN "NA" WHEN (tum.WGS_Dedup_Coverage < 80) THEN "FAIL" WHEN (norm.WGS_Dedup_Coverage < 40) THEN "FAIL"  WHEN (rna.Raw_Reads_Count < 100000000) THEN "FAIL"  ELSE "PASS" END WGS_WTS_PASS
FROM Samples 
LEFT JOIN status on Samples.Sample=status.patient
LEFT JOIN KEY_METRICS norm ON Samples.DNA_N = norm.Sample
LEFT JOIN KEY_METRICS tum ON Samples.DNA_T = tum.Sample
LEFT JOIN KEY_METRICS rna ON Samples.RNA = rna.Sample
;
' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/Sample_Summary.csv

# Report to give all of the metrics related to each Cohort.
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'SELECT Samples.Cohort AS Cohort,
SUM(CASE WHEN (Samples.DNA_N_True == "NA") THEN 0 ELSE 1 END) AS WGS,
SUM(CASE WHEN (Samples.RNA_True == "NA") THEN 0 ELSE 1 END) AS WTA,
SUM(CASE WHEN (status.dna_pipeline_execution == "Complete") THEN 1 ELSE 0 END) AS WGS_Complete,
SUM(CASE WHEN (status.rna_pipeline_light_execution == "Complete") THEN 1 ELSE 0 END) AS WTSLight_Complete,
SUM(CASE WHEN (status.dna_pipeline_execution == "Complete") AND (tum.WGS_Dedup_Coverage > 80) AND (norm.WGS_Dedup_Coverage > 30) THEN 1 ELSE 0 END) AS WGS_Complete_Pass,
SUM(CASE WHEN (status.rna_pipeline_light_execution == "Complete") AND (rna.Raw_Reads_Count > 100000000) THEN 1 ELSE 0 END) AS WTSLight_Complete_Pass,
SUM(CASE WHEN (status.dna_pipeline_execution == "Complete") AND (tum.WGS_Dedup_Coverage > 80) AND (norm.WGS_Dedup_Coverage > 30) AND (status.rna_pipeline_light_execution == "Complete") AND (rna.Raw_Reads_Count > 100000000) THEN 1 ELSE 0 END) AS WTS_WGS_Complete_Pass
FROM Samples
LEFT JOIN status on Samples.Sample=status.patient
LEFT JOIN KEY_METRICS norm ON Samples.DNA_N = norm.Sample
LEFT JOIN KEY_METRICS tum ON Samples.DNA_T = tum.Sample
LEFT JOIN KEY_METRICS rna ON Samples.RNA = rna.Sample
GROUP BY Samples.Cohort
;
' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/Cohort_Summary.csv

# Report to give all of the metrics related to each Instituion.
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'SELECT Samples.Instituion AS Instituion,
SUM(CASE WHEN (Samples.DNA_N_True == "NA") THEN 0 ELSE 1 END) AS WGS,
SUM(CASE WHEN (Samples.RNA_True == "NA") THEN 0 ELSE 1 END) AS WTA,
SUM(CASE WHEN (status.dna_pipeline_execution == "Complete") THEN 1 ELSE 0 END) AS WGS_Complete,
SUM(CASE WHEN (status.rna_pipeline_light_execution == "Complete") THEN 1 ELSE 0 END) AS WTSLight_Complete,
SUM(CASE WHEN (status.dna_pipeline_execution == "Complete") AND (tum.WGS_Dedup_Coverage > 80) AND (norm.WGS_Dedup_Coverage > 30) THEN 1 ELSE 0 END) AS WGS_Complete_Pass,
SUM(CASE WHEN (status.rna_pipeline_light_execution == "Complete") AND (rna.Raw_Reads_Count > 100000000) THEN 1 ELSE 0 END) AS WTSLight_Complete_Pass,
SUM(CASE WHEN (status.dna_pipeline_execution == "Complete") AND (tum.WGS_Dedup_Coverage > 80) AND (norm.WGS_Dedup_Coverage > 30) AND (status.rna_pipeline_light_execution == "Complete") AND (rna.Raw_Reads_Count > 100000000) THEN 1 ELSE 0 END) AS WTS_WGS_Complete_Pass
FROM Samples
LEFT JOIN status on Samples.Sample=status.patient
LEFT JOIN KEY_METRICS norm ON Samples.DNA_N = norm.Sample
LEFT JOIN KEY_METRICS tum ON Samples.DNA_T = tum.Sample
LEFT JOIN KEY_METRICS rna ON Samples.RNA = rna.Sample
GROUP BY Samples.Instituion
;
' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/Institution_Summary.csv


# Report to give all of the metrics related to each Instituion.
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'SELECT Samples.Sample_True as Patient,
norm.WGS_Dedup_Coverage as Normal_Deduplex_coverage,
tum.WGS_Dedup_Coverage as Tumour_Deduplex_coverage
FROM Samples
INNER JOIN KEY_METRICS norm
ON Samples.DNA_N = norm.Sample
INNER JOIN KEY_METRICS tum
ON Samples.DNA_T = tum.Sample
WHERE RNA == "NA"
;
' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/Coverage_WO_WTS.csv
