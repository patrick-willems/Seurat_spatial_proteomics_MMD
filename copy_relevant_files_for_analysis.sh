mkdir /mnt/d/Data/2023-07-CarolineAsselman-FIm/Analysis_v2/images
mkdir /mnt/d/Data/2023-07-CarolineAsselman-FIm/Analysis_v2/images/A1_ROI1
mkdir /mnt/d/Data/2023-07-CarolineAsselman-FIm/Analysis_v2/images/A1_ROI2
find /mnt/f/MACSima/PPP\ MICS_2023-08-01_Impens_moyabrain_human_custom\&REA/Impens_moyabrain_human/R6/A1/ROI1 -type f -name "*_S_*" ! -name "*C-000*" ! -name "*C-001*" ! -name "*C-002*" ! -name "*C-003*" \( -name "*C-004*" -o ! -name "*DAPI*" \) -exec cp {} /mnt/d/Data/2023-07-CarolineAsselman-FIm/Analysis_v2/images/A1_ROI1 \;
find /mnt/f/MACSima/PPP\ MICS_2023-08-01_Impens_moyabrain_human_custom\&REA/Impens_moyabrain_human/R6/A1/ROI2 -type f -name "*_S_*" ! -name "*C-000*" ! -name "*C-001*" ! -name "*C-002*" ! -name "*C-003*" \( -name "*C-004*" -o ! -name "*DAPI*" \) -exec cp {} /mnt/d/Data/2023-07-CarolineAsselman-FIm/Analysis_v2/images/A1_ROI2 \;
