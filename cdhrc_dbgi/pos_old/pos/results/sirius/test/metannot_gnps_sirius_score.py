# %%
import pandas as pd
import io 

#we defined the variable for the paths
input_sirius_path = 'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\sirius\compound_identifications.tsv'

output_sirius_path = 'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\sirius\compound_identifications_sirus_cleaned.tsv'

intput_metannot_path = 


# open compound_identifications.tsv (SIRIUS file)
df = pd.read_csv(input_sirus_path, sep='\t') #change the path to your respective file
suffix = "_sirius"

#adding _sirius suffix
for col in df.columns:
    if col != "id":
        new_col = col + suffix
        df.rename(columns={col: new_col}, inplace=True)

#for each element in column id, keep only the characters after the last '_'
df['id'] = df['id'].str.split('_').str[-1]

#Clean sirius file and save
# List of columns to keep
columns_to_keep = ['id', 'ConfidenceScore_sirius', 'CSI:FingerIDScore_sirius', 'ZodiacScore_sirius', 'SiriusScore_sirius', 'molecularFormula_sirius', 'adduct_sirius', 'InChIkey2D_sirius', 'InChI_sirius', 'name_sirius', 'smiles_sirius', 'xlogp_sirius', 'ionMass_sirius']

# Drop unwanted columns
columns_to_drop = [col for col in df.columns if col not in columns_to_keep]
df = df.drop(columns_to_drop, axis=1)

#save the new dataframe as a tsv file under the name compound_identifications_cleaned.tsv
df.to_csv(r, sep='\t', index=False)


#open the met_annont_enhancer file 
df2 = pd.read_csv(r'C:\\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\met_annot_enhancer\cdhrc_dbgi_alc_pos\cdhrc_dbgi_alc_pos_spectral_match_results_repond.tsv', sep='\t')
suffix = "_metAnnot"

for col in df2.columns:
    if col != "feature_id":
        new_col = col + suffix
        df2.rename(columns={col: new_col}, inplace=True)

#save the new dataframe as a tsv file under the name compound_identifications_cleaned.tsv
df.to_csv(r'C:\\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\met_annot_enhancer\cdhrc_dbgi_alc_pos\cdhrc_dbgi_alc_pos_spectral_match_results_repond_cleaned.tsv', sep='\t', index=False)
df2.head()

#open canopus_compound_summary.tsv 
df3 = pd.read_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\sirius\canopus_formula_summary.tsv', sep='\t')
suffix = "_canopus"
for col in df3.columns:
    if col != "id":
        new_col = col + suffix
        df3.rename(columns={col: new_col}, inplace=True)

# List of columns to keep
columns_to_keep = ['id', 'molecularFormula_canopus', 'adduct_canopus', 'NPC#pathway_canopus', 'NPC#superclass_canopus', 'NPC#class_canopus']

# Drop unwanted columns
columns_to_drop = [col for col in df3.columns if col not in columns_to_keep]
df3 = df3.drop(columns_to_drop, axis=1)
df3.head()

#for each element in column id, keep only the characters after the last '_'
df3['id'] = df3['id'].str.split('_').str[-1]

#save the new dataframe as a tsv file under the name canopus_compound_summary_cleaned.tsv
df3.to_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\sirius\canopus_formula_summary.tsv', sep='\t', index=False)


## MERGING
# Read the SIRIUS input files
df1 = pd.read_csv(r'C:\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\sirius\compound_identifications_sirus_cleaned.tsv', sep='\t')

# Merge the SIRIUS df with MetAnnot df based on the "id" or "feature_id" column
merged_df = pd.merge(df1, df2, left_on='id', right_on='feature_id', how='outer')

# Conserved only one column with the ID 
merged_df['merge_id'] = merged_df['id'].fillna(merged_df['feature_id'])

merged_df = merged_df.drop(['id', 'feature_id'], axis=1)

merged_df['merge_id'] = merged_df['merge_id'].astype(int)

# Saved the merged dataframe in a new file 
merged_df.to_csv(r'C:\\Users\Lecab\largefiles\cdhrc_dbgi\pos\merge_met_annot_sirius\merged_met_annot_sirius.tsv', sep='\t', index=False)  # Replace with the name of you output file 

df3 = pd.read_csv(r'C:\\Users\Lecab\largefiles\cdhrc_dbgi\pos\sirius_output\canopus_compound_summary_cleaned.tsv', sep='\t')
df4 = pd.read_csv(r'C:\\Users\Lecab\largefiles\cdhrc_dbgi\pos\merge_met_annot_sirius\merged_met_annot_sirius.tsv', sep='\t')  # Replace with the name of you output file 
df4.head()

# Merge the dataframes based on the "id" or "feature_id" column
merged_df = pd.merge(df3, df4, left_on='id', right_on='merge_id', how='outer')

# Conserved only one column with the ID 
merged_df['shared_id'] = merged_df['id'].fillna(merged_df['merge_id'])

merged_df = merged_df.drop(['id', 'merge_id'], axis=1)

merged_df['shared_id'] = merged_df['shared_id'].astype(int)

# Saved the merged dataframe in a new file 
merged_df.to_csv(r'C:\\Users\Lecab\largefiles\cdhrc_dbgi\pos\merge_met_annot_sirius\merged_met_annot_sirius_canopus.tsv', sep='\t', index=False)  # Replace with the name of you output file 

df5 = pd.read_csv(r'C:\\Users\Lecab\largefiles\cdhrc_dbgi\pos\merge_met_annot_sirius\merged_met_annot_sirius_canopus.tsv', sep='\t')
# %%
#open quantification_table_reformatted corresponding to your gnps job (csv)
df = pd.read_csv(r'C:\\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\met_annot_enhancer\41f6d201e79f4575b890f958f89e1084\quantification_table_reformatted\e5c06da409d8433eac9e35debe3498ea.csv') #change the path to your respective file
suffix = "_GNPS"

for col in df.columns:
    if col != "row ID":
        new_col = col + suffix
        df.rename(columns={col: new_col}, inplace=True)
df.head()

# %%
#open DB_result data corresponding to your gnps job (csv)
df1 = pd.read_csv(r'C:\\Users\Lecab\gitrepo\plantes_ornementales\ornamental\cdhrc_dbgi\pos\results\met_annot_enhancer\41f6d201e79f4575b890f958f89e1084\DB_result\18887550abb14ca89b2edad740451c83.tsv', sep='\t') #change the path to your respective file
suffix = "_GNPS"

for col in df1.columns:
    if col != "#Scan#":
        new_col = col + suffix
        df1.rename(columns={col: new_col}, inplace=True)
column_headers = df1.columns
print(column_headers)

# %%
#convert the dataframe csv to tsv
tsv_data = df.to_csv(sep='\t', index=False)
df_tsv = pd.read_csv(io.StringIO(tsv_data), sep='\t')

df_tsv.head()

# %%
#merged both GNPS dataframe with same ID and #scan#

merged_GNPS_df = pd.merge(df_tsv, df1, left_on='row ID', right_on='#Scan#', how='outer')

# Conserved only one column with the ID 
merged_GNPS_df['merge_id'] = merged_GNPS_df['row ID'].fillna(merged_GNPS_df['#Scan#'])

merged_GNPS_df = merged_GNPS_df.drop(['row ID', '#Scan#'], axis=1)

merged_GNPS_df['merge_id'] = merged_GNPS_df['merge_id'].astype(int)

# Saved the merged dataframe in a new file 
merged_GNPS_df.to_csv(r'C:\\Users\Lecab\largefiles\cdhrc_dbgi\pos\merge_met_annot_sirius\merged_GNPS.tsv', sep='\t', index=False)  # Replace with the name of you output file 
merged_GNPS_df.head()

# %%
# Read the merged GNPS with MetAnnot_sirius_canopus files dataframe 
df1 = pd.read_csv(r'C:\\Users\Lecab\largefiles\cdhrc_dbgi\pos\merge_met_annot_sirius\merged_met_annot_sirius_canopus.tsv', sep='\t')
df1.head()

# %%
# Merge the dataframes based on the "id" or "feature_id" column
merged_df = pd.merge(merged_GNPS_df, df1, left_on='merge_id', right_on='shared_id', how='outer')

# Conserved only one column with the ID 
merged_df['ID'] = merged_df['merge_id'].fillna(merged_df['shared_id'])

merged_df = merged_df.drop(['merge_id', 'shared_id'], axis=1)

merged_df['ID'] = merged_df['ID'].astype(int)

# Saved the merged dataframe in a new file 
merged_df.to_csv(r'C:\\Users\Lecab\largefiles\cdhrc_dbgi\pos\merge_met_annot_sirius\merged_met_annot_sirius_canopus_GNPS.tsv', sep='\t', index=False)  # Replace with the name of you output file 


# %%
df5 = pd.read_csv(r'C:\\Users\Lecab\largefiles\cdhrc_dbgi\pos\merge_met_annot_sirius\merged_met_annot_sirius_canopus_GNPS.tsv', sep='\t')
column_types = {'InChIkey2D_sirius': str, 'InChIKey_Planar_GNPS': str, 'short_inchikey_metAnnot': str}
for column in df5.columns:
    df5[column] = df5[column].str.replace("-", "_")


# %%
input_file = 'C:\\Users\Lecab\largefiles\cdhrc_dbgi\pos\merge_met_annot_sirius\merged_met_annot_sirius_canopus_GNPS.tsv'

# Définir les noms des colonnes à extraire
columns_to_extract = ["InChIkey2D_sirius", "InChIKey-Planar_GNPS", "short_inchikey_metAnnot"]

# Charger uniquement les colonnes spécifiées dans un DataFrame
df6 = pd.read_csv(input_file, usecols=columns_to_extract)

# Afficher le DataFrame résultant
print(df6)


# %%
#Define the input file path
input_file = 'C:\\Users\Lecab\largefiles\cdhrc_dbgi\pos\merge_met_annot_sirius\merged_met_annot_sirius_canopus_GNPS.tsv'

# Define the column names
inchikey2d_column = 'InChIkey2D_sirius'
inchikey_planar_column = 'InChIKey_Planar_GNPS'
short_inchikey_column = 'short_inchikey_metAnnot'

# Define the scoring function
def calculate_score(inchikey2d_sirius, inchikey_planar_gnps, short_inchikey_metannot):
    if inchikey2d_column and inchikey_planar_column and short_inchikey_column:
        if inchikey2d_column == inchikey_planar_column == short_inchikey_column:
            return 3
        elif inchikey2d_column == inchikey_planar_column or inchikey2d_column == short_inchikey_column or inchikey_planar_column == short_inchikey_column:
            return 2
        else:
            return 1
    else:
        return 1
    

# Apply the scoring function to create a new 'score' column
df['score'] = df.apply(lambda row: calculate_score(row[inchikey2d_column], row[inchikey_planar_column], row[short_inchikey_column]), axis=1)

# Display the updated DataFrame
print(df)




