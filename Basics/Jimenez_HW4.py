# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 21:57:11 2025

@author: Maria  G. Jimenez Velasco
Course: Post-genomic

This code take csv file and extract basic transcriptomic information
"""
import pandas as pd
import os


#####################################################################
#                         PART 1                                    #
######################################################################

path = r"C:\Users\Maria\Documents\Bioinformatics_Maria\Classes_Fall_2025\Postgenomics\HWS\Lab4\files"
def csv2tsv(path, filename, dest_filename):
    # reg_type: +1 for upregulated, -1 for down regulated
    
    
    # Read files
    file = os.path.join(path, filename)
    normal_df = pd.read_csv(file)
    print(normal_df.shape)
    print(normal_df.columns)
    
    
    #extract specific columns of interest
    normal = normal_df[['chrom', 'left','ref_seq', 'alt_seq']]

    
    # add the strand with up regulated or downregulated
    normal_pos_strand = normal.copy()
    normal_pos_strand ["strand"] = 1
    
    normal_neg_strand  = normal.copy()
    normal_neg_strand["strand"] = -1
    
    # normal_pos_strand["strand"] = normal_pos_strand["strand"].astype(str)
    # normal_neg_strand["strand"] = normal_neg_strand["strand"].astype(str)
    
    
    
    # split the chromosome column into type and id
    id_num  = normal_pos_strand["chrom"].str.replace("^chr", "", regex=True)
    normal_pos_strand.insert(1, "type", id_num)  # extract number OR letter
    normal_pos_strand["chrom"] = "chr"  # keep only "chr"
    # rename the columns 
    normal_pos_strand.columns = ["type", "id", "position", "allele1", "allele2", "strand"]
    
    
    id_num  = normal_neg_strand["chrom"].str.replace("^chr", "", regex=True)
    normal_neg_strand.insert(1, "type", id_num)  # extract number OR letter
    normal_neg_strand["chrom"] = "chr"  # keep only "chr"
    # rename the columns 
    normal_neg_strand.columns = ["type", "id", "position", "allele1", "allele2", "strand"]
    normal_neg_strand['type'] = "chromosome"
    normal_pos_strand['type'] = "chromosome" 
    
    
    print("normal pos: ", normal_pos_strand.shape)
    print("normal neg:", normal_neg_strand.shape)
          
    normal = pd.concat([normal_pos_strand, normal_neg_strand], ignore_index=True)
    print("final shape: ", normal.shape)
    
    # save them as .tsv files
    dest_path = os.path.join(path, dest_filename+"_pos_strand.tsv")
    normal_pos_strand.to_csv(dest_path, sep="\t", index=False)
    dest_path = os.path.join(path, dest_filename+"_neg_strand.tsv")
    normal_neg_strand.to_csv(dest_path, sep="\t", index=False)
    
    dest_path = os.path.join(path, dest_filename+".tsv")
    normal.to_csv(dest_path, sep="\t", index=False)
    
    
    return 

csv2tsv(path, filename = "Final_Normal.csv", dest_filename="normal")
csv2tsv(path, filename = "Final_Tumor.csv", dest_filename="tumor")

#   in progress since it is nor required in this case for the lab
# def split_batches(path, filename, batch_size):
#     file = os.path.join(path, filename)
#     df = pd.read_csv(file)
    
#     # Input TSV file
#     input_file = "variants.tsv"
    
#     # Load dataframe
#     df = pd.read_csv(input_file, sep="\t")
#     # Get base filename (without extension)
#     base_name, ext = os.path.splitext(input_file)
    
#     # Calculate number of chunks
#     num_chunks = math.ceil(len(df) / chunk_size)
    
#     # Split and save
#     for i in range(num_chunks):
#         start = i * chunk_size
#         end = start + chunk_size
#         chunk = df.iloc[start:end]
        
#         # Save with suffix (_1.tsv, _2.tsv, etc.)
#         output_file = f"{base_name}_{i+1}.tsv"
#         chunk.to_csv(output_file, sep="\t", index=False)
    
#         print(f"Saved {output_file} with {len(chunk)} variants")
    


path = r"C:\Users\Maria\Documents\Bioinformatics_Maria\Classes_Fall_2025\Postgenomics\HWS\Lab4\files"


def remove_downregulated(path, filename, reg_type, ext):
    # reg_type : '-1' for downregulated or '+1' for upregulated
    # ext: xtemsion
    file = os.path.join(path, filename+ext)
    # Load your CSV file
    df = pd.read_csv(file, sep="\t")  # assuming it's tab-separated
    df  = df.reset_index(drop=True)
    
    # Keep only rows where 'id' does NOT end with '-1'
    df_filtered = df[~df['id'].str.endswith(reg_type)]
    df_filtered  = df_filtered.reset_index(drop=True)
    
    dest_file = os.path.join(path, filename+"_filtered.csv")
    # Save the filtered dataframe to a new CSV
    df_filtered.to_csv(dest_file, sep=",", index=False)
    print(filename," before filetered " , df.shape)
    print(filename, " filtered: ", df_filtered.shape)
    
    # print(f"Original rows: {len(df)}, Filtered rows: {len(df_filtered)}")
    return df_filtered

upreg_normal = remove_downregulated(path, "normal_pervariant", '-1', ".tsv") # remove downregulated, +1remove upregulated
upreg_tumor = remove_downregulated(path, "tumor_2_pervariant", '-1', ".tsv") # 

def remove_downregulated_vcf(path, filename, reg_type, ext):
    input_vcf = os.path.join(path, filename+ext)
    # Load your CSV file
    header_lines = []
    data_rows = []
    
    # read the vcf file
    with open(input_vcf) as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                parts = line.strip().split("\t") 
                data_rows.append(parts)
    
    # remove the rows where id ends with -1
    filtered_rows = [row for row in data_rows if not row[2].endswith(reg_type)]
    
    print(f"Removed {len(data_rows) - len(filtered_rows)} rows with ID ending in '-1'")
    
    # --- Step 3: Write new VCF ---
    with open(input_vcf[:-4]+"_filtered.vcf", "w") as f:
        f.writelines(header_lines)
        for row in filtered_rows:
            f.write("\t".join(row) + "\n")
    return
    
upreg_normal_vcf = remove_downregulated_vcf(path, "gen_coords_normal", '-1', ".vcf") # remove downregulated, +1remove upregulated
upreg_tumor_vcf = remove_downregulated_vcf(path, "gen_coords_tumor2", '-1', ".vcf") # 




#  merging two files in columns
def merge_cols(path, file1, file2, filename):
    # file 1 is the name of a original file with extension
    # file 2 is a dataframe (upregulated file)
    # filename: dest file name with extensions. with expect csv
    # read the original file obtained from Lab 3
    
    file = os.path.join(path, file1)
    # Load your CSV file
    df = pd.read_csv(file, sep=",")  
    df  = df.reset_index(drop=True)
    
    file2 = file2.reset_index(drop=True)
    
    print("\n",file1, " shape : ", df.shape)
    print("filetered upregulared shape : ", file2.shape)
    merged_df = pd.concat([df, file2], axis=1) # origina with filtered file
    merged_df  = merged_df.reset_index(drop=True)
    
    print("merged ", merged_df.shape)
    
        #  merge both files
    file = os.path.join(path, filename)
    merged_df.to_csv(file, sep=",", index=False)
    
    merged = pd.merge(df, file2, left_on="left", right_on="position", how="inner")
    merged.to_csv(file[:-4]+"2.csv", sep=",", index=False)
    # print(merged["position"].duplicated().sum())
    return 

merge_cols(path, "Final_Normal.csv", upreg_normal, "merged_normal_ori_snpnexus.csv")
merge_cols(path, "Final_Tumor.csv", upreg_tumor, "merged_tumor_ori_snpnexus.csv")




# ##################################################################
#                           part  2                                #
####################################################################

# read txt file

def merge_ori_snpnexus_fathmm(path, fathm_filename, ori_snpnexus_file, dest_filename):
    fathmm_file = os.path.join(path, fathm_filename)
    # Load fathmm file
    df1 = pd.read_csv(fathmm_file, sep="\t")  
    df1  = df1.reset_index(drop=True)
    print("fathmm shape ", df1.shape)
    
    ori_nexus_file = os.path.join(path, ori_snpnexus_file)
    df2 = pd.read_csv(ori_nexus_file, sep=",")
    df2  = df2.reset_index(drop=True)
    # df2.to_csv(df2, sep=",", index=False)
    print("ori with snpnexus shape", df2.shape)
    
    merged_df = pd.concat([df2, df1], axis=1) # original with snpnexus and fathmm results
    merged_df  = merged_df.reset_index(drop=True)
    
    print("final merged  shape : ", merged_df.shape)
    file = os.path.join(path,  dest_filename)
    merged_df.to_csv(file, sep=",", index=False)
    
    return
    


print("\nNormal: ")
merge_ori_snpnexus_fathmm(path, "fathmm_result_normal_filtered.txt", "merged_normal_ori_snpnexus.csv", "Jimenez_normal_HW4.csv")

print("\nTumor: ")
merge_ori_snpnexus_fathmm(path, "fathmm_result_tumor2_filtered.txt", "merged_tumor_ori_snpnexus.csv", "Jimenez_tumor_HW4.csv")



