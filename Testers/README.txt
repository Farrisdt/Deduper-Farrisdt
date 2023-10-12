All files are based on the first 164 lines of /projects/bgmp/ftedder/bioinfo/Bi623/PS/QAA/Aligned.out.sam
	deduper.test.sorted.sam
		The sorted sam file for input, 100 reads long with all headers. 
		Mix of reasons for not being duplicates as well as placement of duplicates (start of chromosome, entire position, etc.)
		Includes various cigars that cause duplication. 
	deduper.test.output.sam
		Expected output file.
		It is possible more lines are removed than seen here as the original file had some natural duplicates and I may have missed some while creating the output.
	deduper.test.list.txt
		List of 100 letters. n for not a duplicate, N for the first of a duplicate, y for a duplicate that needs to be removed, and e for unknown UMIs.
		It is possible this is missing some duplicates, see above, and there may be some off by 1 errors in the lines.
		Possible future improvement if needed: adding what makes the line a duplicate or not a duplicate in case patterns emerge while debugging.

		I believe all errors in testing files will be quickly ironed out once they start being used.
