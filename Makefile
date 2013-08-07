clean:
	rm -rf adapters_out barcodes_out aligned aligned_split plots collapsed_out
	mkdir adapters_out barcodes_out aligned aligned/samfiles aligned/bamfiles aligned_split aligned_split/forward aligned_split/reverse plots collapsed_out 
stash:
	tar -cvjf stashed_pipeline.tar.bz2 .
remake:
	tar -cjvf test.tar.bz2 *.json barcodes bin docs Makefile aligndepth alignplot demultiplex start_pipeline
