# for each gene in list of final genes (after all filtering steps)
cat list_final.txt | while read gene
do
	# create output folder
	mkdir codeml_out/${gene}
	# copy control file
	cp dnds.ctl codeml_out/${gene}/
	# replace initial gene ID with ID of current gene
	sed -i '' -e "s/66at32523/${gene}/g" codeml_out/${gene}/dnds.ctl
	cd codeml_out/${gene}/
	# run codeML as defined in control file, keep output in logfile
	/Users/jule/Desktop/dNdS_analysis/paml4.8/bin/codeml dnds.ctl | tee logfile.txt
	cd ../../
done