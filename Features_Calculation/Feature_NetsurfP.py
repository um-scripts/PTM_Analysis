import biolib
nsp3 = biolib.load('DTU/NetSurfP-3')
nsp3_results = nsp3.cli(args='-i Phosphorylation.txt')
nsp3_results.save_files("PhosphorylationNewBiolib_results/")
