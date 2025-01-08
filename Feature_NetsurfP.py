# Import Biolib and load NetSurfP-3.0
import biolib
nsp3 = biolib.load('DTU/NetSurfP-3')
# Run NetSurfP-3.0 for our query 
nsp3_results = nsp3.cli(args='-i Phosphorylation.txt')
# Optionally save the results
nsp3_results.save_files("PhosphorylationNewBiolib_results/")
