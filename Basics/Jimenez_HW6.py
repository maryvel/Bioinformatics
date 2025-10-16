# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 14:36:20 2025

@author: Maria G. Jimenez
Postgenomic
Homework 6
"""
import pandas as pd
import os

pd.set_option('display.max_columns', None)  # Show all column headers
pd.set_option('display.width', 0)           # Auto-detect width (don’t truncate header)
pd.set_option('display.max_colwidth', None) # Don’t truncate column names

path = r"C:\Users\Maria\Documents\Bioinformatics_Maria\Classes_Fall_2025\Postgenomics\HWS\Lab6\files"

# 1.1 Using the Normal and Tumor CSV files from Homework 3. Subset the two CSV
# files with only the columns, ["chrom", "left", "ref_seq", "alt_seq", 
# “Patient_ID”, ‘VCF_ID”]
def subset_file(path, filename, cols):
    file = os.path.join(path, filename)
    df = pd.read_csv(file, sep=",")  # assuming it's tab-separated
    df  = df.reset_index(drop=True)
    subset_file = df[cols].copy()
    print(subset_file.head())
    return subset_file
    
def extract_data(df, file_type):
    # 1.1.1	How many unique normal patients do we have?
    # 1.1.2	How many unique tumor patients do we have?
    unique_patients = df["Patient_ID"].nunique()
    print("Question 1.1.1 or 1.1.2: Unique", file_type, " patients:", unique_patients)
   
    # 1.1.3	Group by variant info, chrom, left, ref_seq, and alt_seq, let the other columns turn into list. 
    group_cols = ["chrom", "left", "ref_seq", "alt_seq"]
    print("\n", file_type, ":  1.1.3: Group by variant info, chrom, left, ref_seq, and alt_seq, ")
    df_grouped = df.groupby(group_cols).agg(list).reset_index()
        

    # # Use nunique for patient count, list for VCF_IDs
    # df_grouped = df.groupby(group_cols).agg({
    #     "Patient_ID": pd.Series.nunique,  # count unique patients
    #     "VCF_ID": list                     # keep list of VCF_IDs
    # }).reset_index()
    
    print(df_grouped.head())   
   # 
    # # 1.1.4	Create a new column with the number of patients per variant on both the normal and tumor (name the column, N# and T#, respectively). 
    print("\n", file_type, ": Question 1.1.4: ")
    if file_type == "Normal":
       df_grouped["N#"] = df_grouped["Patient_ID"].apply(lambda x: len(set(x)))
    elif file_type == "Tumor":
       df_grouped["T#"] = df_grouped["Patient_ID"].apply(lambda x: len(set(x)))
    else:
       raise ValueError("file_type must be 'Normal' or 'Tumor'")

    print("\n ",df_grouped.head())
    # 1.1.5	Rename the columns, Patient_ID and VCF_ID, to have, _Normal or _Tumor, added depending which file you are working with. 
    print("\n", file_type, ":Question 1.1.5: ")
    df_grouped.rename(columns={"Patient_ID": "Patient_ID_"+file_type, "VCF_ID": "VCF_ID_"+file_type}, inplace=True)
    print(df_grouped.head())
       
    return df_grouped


# Load your CSV file
cols = ["chrom", "left", "ref_seq", "alt_seq", "Patient_ID", "VCF_ID"]

# subset columns
filename = "Final_Normal.csv"
normal_subset = subset_file(path, filename, cols)
filename = "Final_Tumor.csv"
tumor_subset = subset_file(path, filename, cols)


print("\nInformation about Normal.csv ")
normal_grouped  = extract_data(normal_subset, "Normal")

print("\n\n ")
print("Information about Tumor.csv")
tumor_grouped = extract_data(tumor_subset, "Tumor")
 

def part2(normal_grouped, tumor_grouped, how="outer", columns =["chrom", "left", "ref_seq", "alt_seq"]):
    normal_grouped = normal_grouped.reset_index(drop=True)
    tumor_grouped = tumor_grouped.reset_index(drop=True)
    
    aml = pd.merge(normal_grouped, tumor_grouped, how="outer", on=columns, suffixes=('_normal', '_tumor'), indicator=True )
    aml.to_csv("AML.csv", index=False)


    # # 1.2.1 & 1.2.2 & 1.2.3 Unique variants and shared
    # unique_normal_variants = normal_grouped.shape[0]
    # unique_tumor_variants = tumor_grouped.shape[0]
    
    # Unique and shared variants

    print("Question 1.2.1	How many unique normal variants?")
    n_unique_normal_variants = (aml['_merge'] == 'left_only').sum()
    print("Unique normal variantes: ", n_unique_normal_variants)
    
    print("\n Question 1.2.2	How many unique tumor variants?")
    n_unique_tumor_variants = (aml['_merge'] == 'right_only').sum()
    print("Unique tumor variants:", n_unique_tumor_variants)
    
    print("\n Question 1.2.3	1.2.3	How many variants are shared between normal and tumor (common)?")
    common_variants = (aml['_merge'] == 'both').sum()
    print("Shared variants:", common_variants)
    

print("\n\n ************* Part 1.2... **********\n")

part2(normal_grouped, tumor_grouped, how="outer", columns =["chrom", "left", "ref_seq", "alt_seq"]  )

# Q 1.3: Using the Normal and Tumor files from Homework 3, concatenate these files
 # along the axis = 0, with this Expand/Explode the rows based on the CSQ columns
 # and save this file as AML_Expand.csv. Remove duplicate rows.


def part3(path, file1, file2):
    file = os.path.join(path, file1)
    normal = pd.read_csv(file, sep=",")  # assuming it's tab-separated
    normal  = normal.reset_index(drop=True)
    
    file = os.path.join(path, file2)
    tumor = pd.read_csv(file, sep=",")  # assuming it's tab-separated
    tumor  = tumor.reset_index(drop=True)
                                 
  
    # concatenate the normal and tumor csv files
    aml_concatenated = pd.concat([normal, tumor], axis=0, ignore_index=True)
    aml_merged = aml_concatenated.drop_duplicates()
    aml_merged  =  aml_merged.reset_index()
    
    print(aml_merged.columns)
    # extraction of the CSQ columns
    csq_columns = ['Allele', 'Consequence', 'IMPACT', 'SYMBOL',
    'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc',
    'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position',
    'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 'DISTANCE',
    'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID',
    'CANONICAL', 'TSL', 'APPRIS', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL',
    'UNIPARC', 'RefSeq', 'GENE_PHENO', 'SIFT', 'PolyPhen', 'DOMAINS',
    'HGVS_OFFSET', 'GMAF', 'AFR_MAF', 'AMR_MAF', 'EAS_MAF', 'EUR_MAF',
    'SAS_MAF', 'AA_MAF', 'EA_MAF', 'ExAC_MAF', 'ExAC_Adj_MAF',
    'ExAC_AFR_MAF', 'ExAC_AMR_MAF', 'ExAC_EAS_MAF', 'ExAC_FIN_MAF',
    'ExAC_NFE_MAF', 'ExAC_OTH_MAF', 'ExAC_SAS_MAF', 'CLIN_SIG', 'SOMATIC',
    'PHENO', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS',
    'MOTIF_SCORE_CHANGE', 'ENTREZ', 'EVIDENCE']
    
    for col in csq_columns:
        for row in range(aml_merged.shape[0]):
            temp = str(aml_merged[col][row])
            temp = temp.strip("[]").replace('"',"").replace("'","").replace(" ", "")
            
            aml_merged[col][row] = temp.split(",")
            if type(aml_merged[col][row])!= list:
                print(aml_merged[col][row])
                print(type(aml_merged[col][row]))
    aml_exploded = aml_merged.explode(csq_columns, ignore_index = True)
    
    
    # save the expanded CSV
    file = os.path.join(path, "AML_Expand.csv")
    aml_exploded.to_csv(file, index=False)
    
    print("1.3.1 How many rows are in this file? \n")
    print("Number of rows in expanded/exploded file:", aml_exploded.shape[0])
    
    print("1.3.1.: Subset of expanded with only the columns, ['SYMBOL', 'Gene', 'Feature']")
    print("... name this AML_gene.csv.")
    
    # Gene- related
    aml_gene = aml_exploded[['SYMBOL', 'Gene', 'Feature']].drop_duplicates()
    file = os.path.join(path, "AML_gene.csv")
    aml_gene.to_csv(file, index=False)
    
    print("1.3.2. Subset of expanded CSV with only the columns, ['chrom', 'left', \
          'right', 'ref_seq', 'alt_seq', 'Feature', 'cDNA_position', 'BIOTYPE'], name this AML_tx.csv.")
    
    # transcript - related
    aml_transcript = aml_exploded[['chrom', 'left', 'right', 'ref_seq', 'alt_seq', 
                      'Feature', 'cDNA_position', 'BIOTYPE']].drop_duplicates()
    file = os.path.join(path, "AML_tx.csv")
    aml_transcript.to_csv(file, index=False)


    return


print("\n\n ************* Part 1.3... ********PENDING**\n")

part3(path, "Final_Normal.csv", "Final_Tumor.csv")

# Allele	Consequence	IMPACT	SYMBOL	Gene	Feature_type	Feature	BIOTYPE	EXON	INTRON	HGVSc	HGVSp	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	ALLELE_NUM	DISTANCE	STRAND	FLAGS	VARIANT_CLASS	SYMBOL_SOURCE	HGNC_ID	CANONICAL	TSL	APPRIS	CCDS	ENSP	SWISSPROT	TREMBL	UNIPARC	RefSeq	GENE_PHENO	SIFT	PolyPhen	DOMAINS	HGVS_OFFSET	GMAF	AFR_MAF	AMR_MAF	EAS_MAF	EUR_MAF	SAS_MAF	AA_MAF	EA_MAF	ExAC_MAF	ExAC_Adj_MAF	ExAC_AFR_MAF	ExAC_AMR_MAF	ExAC_EAS_MAF	ExAC_FIN_MAF	ExAC_NFE_MAF	ExAC_OTH_MAF	ExAC_SAS_MAF	CLIN_SIG	SOMATIC	PHENO	PUBMED	MOTIF_NAME	MOTIF_POS	HIGH_INF_POS	MOTIF_SCORE_CHANGE	ENTREZ	EVIDENCE

def Random_Forest(t_size = 0.2, num_trees= 100, maxmun_depth= 3):
    # import required libraries
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
 
    
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import accuracy_score, confusion_matrix
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.datasets import load_iris

    #Import Dataset
    iris = load_iris()
    df = pd.DataFrame(data=iris.data, columns=iris.feature_names)
    df['target'] = iris.target
    # print(df)
    
    
    # Data preparation
    # Here we will separate the features (X) and the target variable (y)
    X = df.iloc[:, :-1].values
    y = df.iloc[:, -1].values
    
    # 4. Splitting the Dataset
    # We'll split the dataset into training and testing sets so we can train the model on one part and evaluate it on another.
    
    # X_train, y_train: 80% of the data used to train the model.
    # X_test, y_test: 20% of the data used to test the model.
    # test_size=0.2: means 20% of data goes to testing.
    # random_state=42: ensures you get the same split every time

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=t_size, random_state=42)
    
    # Feature Scaling
    # Feature scaling ensures that all the features are on a 
    # similar scale which is important for some machine learning models.
    # However Random Forest is not highly sensitive to feature scaling. But it is a good practice to scale when combining models.
    
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)
    
    # Building Random Forest Classifier
    # RandomForestClassifier(n_estimators=100, random_state=42) creates 100 trees (100 trees balance accuracy and training time).
    # classifier.fit(X_train, y_train) trains on training data.
    # classifier.predict(X_test) predicts on test data.
    # random_state=42 ensures reproducible results.

    classifier = RandomForestClassifier(n_estimators= num_trees, max_depth = maxmun_depth, random_state=42)
    classifier.fit(X_train, y_train)
    y_pred = classifier.predict(X_test)
    
    
    #7. Evaluation of the Model
    # We will evaluate the model using the accuracy score and confusion matrix.
    print("\n 2.1 What was the accuracy of the model you built?\n")
    accuracy = accuracy_score(y_test, y_pred)
    print(f'Accuracy: {accuracy * 100:.2f}%')
    
    conf_matrix = confusion_matrix(y_test, y_pred)
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(conf_matrix, annot=True, fmt='g', cmap='Blues', cbar=False, 
                xticklabels=iris.target_names, yticklabels=iris.target_names)
    
    plt.title('Confusion Matrix Heatmap')
    plt.xlabel('Predicted Labels')
    plt.ylabel('True Labels')
    plt.show()
    
    
    print("\n 2.2 What order were the important features ranked?")
    # Feature Importance
    importances = classifier.feature_importances_
    feature_names = iris.feature_names
    
    # Combine and sort
    feature_importance_df = pd.DataFrame({'Feature': feature_names, 'Importance': importances})
    feature_importance_df = feature_importance_df.sort_values(by='Importance', ascending=False)
    
    print("Feature importance ranking:")
    print(feature_importance_df)
    
    # Optional: visualize
    plt.figure(figsize=(8, 5))
    sns.barplot(x='Importance', y='Feature', data=feature_importance_df)
    # sns.barplot(x='Importance', y='Feature', data=feature_importance_df, palette='viridis')
    plt.title('Feature Importance Ranking')
    plt.show()


    return 

print("\n\n ************* Part 2:..Random Forest ********\n")
Random_Forest(t_size = 0.2, num_trees= 100, maxmun_depth= None)

print("\n2.3 Change 2 of the parameters and repeat the model generation. ")
print("Changing the max_depth parameter from None to  3, and...")
print("changing the number of trees (number of estimators) from 100 to 25")
Random_Forest(t_size = 0.2, num_trees= 25, maxmun_depth= 3)


print("\n2.3 Change 2 of the parameters and repeat the model generation. ")
print("Chaning the number of estimators (or trees) from 100 to 50, and...")
print("changing the test_size from 20% (0.2) to 50% (0.5).")
Random_Forest(t_size = 0.5, num_trees= 50, maxmun_depth= None)


print("\n\n ************* Part 3:..K-Means ********\n")

def K_Means(num_clusters = 3):
    # https://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_assumptions.html
    # Authors: The scikit-learn developers
    # SPDX-License-Identifier: BSD-3-Clause
    
    # Data Generation
    import numpy as np

    from sklearn.datasets import make_blobs
    
    n_samples = 1500
    random_state = 170
    transformation = [[0.60834549, -0.63667341], [-0.40887718, 0.85253229]]
    
    X, y = make_blobs(n_samples=n_samples, random_state=random_state)
    X_aniso = np.dot(X, transformation)  # Anisotropic blobs
    X_varied, y_varied = make_blobs(
        n_samples=n_samples, cluster_std=[1.0, 2.5, 0.5], random_state=random_state
    )  # Unequal variance
    X_filtered = np.vstack(
        (X[y == 0][:500], X[y == 1][:100], X[y == 2][:10])
    )  # Unevenly sized blobs
    y_filtered = [0] * 500 + [1] * 100 + [2] * 10

    #Data visualization

    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))
    
    axs[0, 0].scatter(X[:, 0], X[:, 1], c=y)
    axs[0, 0].set_title("Mixture of Gaussian Blobs")
    
    axs[0, 1].scatter(X_aniso[:, 0], X_aniso[:, 1], c=y)
    axs[0, 1].set_title("Anisotropically Distributed Blobs")
    
    axs[1, 0].scatter(X_varied[:, 0], X_varied[:, 1], c=y_varied)
    axs[1, 0].set_title("Unequal Variance")
    
    axs[1, 1].scatter(X_filtered[:, 0], X_filtered[:, 1], c=y_filtered)
    axs[1, 1].set_title("Unevenly Sized Blobs")
    
    plt.suptitle("Ground truth clusters").set_y(0.95)
    plt.show()
    
    # Fit models and plot results
    from sklearn.cluster import KMeans

    common_params = {
        "n_init": "auto",
        "random_state": random_state,
    }
    
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))
    
    y_pred = KMeans(n_clusters=2, **common_params).fit_predict(X)
    axs[0, 0].scatter(X[:, 0], X[:, 1], c=y_pred)
    axs[0, 0].set_title("Non-optimal Number of Clusters")
    
    y_pred = KMeans(n_clusters=num_clusters, **common_params).fit_predict(X_aniso)
    axs[0, 1].scatter(X_aniso[:, 0], X_aniso[:, 1], c=y_pred)
    axs[0, 1].set_title("Anisotropically Distributed Blobs")
    
    y_pred = KMeans(n_clusters=num_clusters, **common_params).fit_predict(X_varied)
    axs[1, 0].scatter(X_varied[:, 0], X_varied[:, 1], c=y_pred)
    axs[1, 0].set_title("Unequal Variance")
    
    y_pred = KMeans(n_clusters=num_clusters, **common_params).fit_predict(X_filtered)
    axs[1, 1].scatter(X_filtered[:, 0], X_filtered[:, 1], c=y_pred)
    axs[1, 1].set_title("Unevenly Sized Blobs")
    
    plt.suptitle("Unexpected KMeans clusters").set_y(0.95)
    plt.show()
    
    # possible solutions:
    y_pred = KMeans(n_clusters=num_clusters, **common_params).fit_predict(X)
    plt.scatter(X[:, 0], X[:, 1], c=y_pred)
    plt.title("Optimal Number of Clusters")
    plt.show()

    # To deal with unevenly sized blobs one can increase the number of random 
    # initializations. In this case we set n_init=10 to avoid finding a 
    # sub-optimal local minimum. For more details see Clustering sparse data with k-means.
    
    y_pred = KMeans(n_clusters=num_clusters, n_init=10, random_state=random_state).fit_predict(
    X_filtered)
    plt.scatter(X_filtered[:, 0], X_filtered[:, 1], c=y_pred)
    plt.title("Unevenly Sized Blobs \nwith several initializations")
    plt.show()
    
    #As anisotropic and unequal variances are real limitations of the k-means
    # algorithm, here we propose instead the use of GaussianMixture, which also
    # assumes gaussian clusters but does not impose any constraints on their variances.
    # Notice that one still has to find the correct number of blobs (see Gaussian
    # Mixture Model Selection)."
    from sklearn.mixture import GaussianMixture

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    
    y_pred = GaussianMixture(n_components=num_clusters).fit_predict(X_aniso)
    ax1.scatter(X_aniso[:, 0], X_aniso[:, 1], c=y_pred)
    ax1.set_title("Anisotropically Distributed Blobs")
    
    y_pred = GaussianMixture(n_components=num_clusters).fit_predict(X_varied)
    ax2.scatter(X_varied[:, 0], X_varied[:, 1], c=y_pred)
    ax2.set_title("Unequal Variance")
    
    plt.suptitle("Gaussian mixture clusters").set_y(0.95)
    plt.show()

    return


print("\n\n ************* Part 3:..K-MEANS ********\n")
print("\n Question 3.1: ")
K_Means()

print("\n Question 3.2: from 3 to 4 ")
K_Means(num_clusters= 4)


print("\n Question 3.2: from 3 to 2: ")
K_Means(num_clusters= 2)
