#!/usr/bin/env bash

sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'select * from File_Locations;' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/File_Locations.csv
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'select * from KEY_METRICS;' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/KEY_METRICS.csv
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'select * from STATUS;' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/STATUS.csv
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'select * from Samples;' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/Samples.csv
sqlite3 -header -csv /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/MOH_analysis.db 'select * from Timestamps;' > /lustre03/project/6007512/C3G/projects/MOH_PROCESSING/DATABASE/CSV/Timestamps.csv





