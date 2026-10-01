#!/usr/bin/env bash
# ==============================================================================
# OASpepDB: Demo Repertoire Downloader for Workflow Reproducibility
# Downloads a curated set of 8 lightweight OAS repertoires (~15 KB total)
# covering Healthy controls and 3 disease cohorts (COVID-19, CLL, HIV).
# For full multi-gigabyte production download (14,433 studies), use:
# bulk_download_human_unpaired.sh
# ==============================================================================

# Healthy Controls (Soto et al., 2019 - PBMC)
wget https://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Soto_2019/csv/SRR8365444_1_Heavy_IGHM.csv.gz
wget https://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Soto_2019/csv/SRR8365442_1_Heavy_IGHE.csv.gz
wget https://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Soto_2019/csv/SRR8365439_1_Heavy_IGHA.csv.gz

# Chronic Lymphocytic Leukemia (Bashford et al., 2013 - PBMC)
wget https://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Bashford_2013/csv/ERR220437_Heavy_Bulk.csv.gz
wget https://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Bashford_2013/csv/ERR220400_Heavy_IGHM.csv.gz

# COVID-19 / SARS-CoV-2 (Woodruff et al., 2020 & Montague et al., 2021)
wget https://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Woodruff_2020/csv/SRR12113363_1_Heavy_IGHE.csv.gz
wget https://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Montague_2021/csv/SRR12190270_Heavy_IGHA.csv.gz

# HIV Infection (Zhu et al., 2013 - PBMC)
wget https://opig.stats.ox.ac.uk/webapps/ngsdb/unpaired/Zhu_2013/csv/SRR924020_1_Heavy_IGHM.csv.gz
