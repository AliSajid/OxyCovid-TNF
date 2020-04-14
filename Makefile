RS=Rscript
DATA_DIR=data/
RES_DIR=results/
FIG_DIR=figures/

.PHONY:
	clean all

clean:
	find $(DATA_DIR) -name "*.tsv" -delete
	find $(RES_DIR) -name "*.tsv" -delete
	find $(RES_DIR) -name "*.csv" -delete
	find $(FIG_DIR) -name "*.png" -delete
	find $(FIG_DIR) -name "*.pdf" -delete

curate: curate_signatures.R
	Rscript curate_signatures.R

process: process_signatures.R
	Rscript process_signatures.R

analyse: analyse_data.R
	Rscript analyse_data.R

visualise: figure*.R
	find . -name "figure*.R" | xargs Rscript

all:
	curate process analyse visualise
