import sys
import pandas as pd

def merge(csv_file_input,csv_file_output):
	#Read CSV file and transform it into a data frame
	pca_input = pd.read_csv(csv_file_input,sep='\t')

	#Merged CSV output file
	pca_output = pd.DataFrame()

	#Collect column headers
	column_header = list(pca_input.columns)

	#Remove the .version of gene ID
	new_geneid = pca_input['Geneid'].apply(lambda x: x.split(".")[0])
	pca_output['Geneid'] = new_geneid

	#Collect maximal length for each gene ID
	new_length = pca_input.groupby(new_geneid,sort=False).max()['Length'].values
	pca_output['Length'] = new_length

	#Merge values of column headers
	for header in column_header[6:]:
		header_values = pca_input.groupby(new_geneid,sort=False).sum()[header].values
		pca_output[header] = header_values

	#Write new CSV file
	pca_output.to_csv(csv_file_output,sep='\t',encoding='utf-8')

	print("Done")

merge(*sys.argv[1:])


		
