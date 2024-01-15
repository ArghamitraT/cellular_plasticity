"""
This file contains all the household common functions; built by AT
"""
import matplotlib.pyplot as plt
import pandas as pd
import csv
import os
import json
import seaborn as sns

#creates an str image name with current date and time signature
def create_image_name(name, format=".png"):
    import datetime
    import time
    crnt_tm = datetime.datetime.now()
    image_name = (name+"_" + str(crnt_tm.year) + "_" + str(crnt_tm.month) + "_" + str(crnt_tm.day) + "_"
                  + time.strftime("%H_%M_%S") + format)
    return image_name

#creates bar plot
def bar_plt_distribution(df, title='Distribution Bar Plot', image_name='Distribution Bar Plot'):
    import os

    plt.bar(df.index, df.values)

    for val in range(len(df.index)):
        df.index.values[val] = df.index[val].replace("cell", "")

    main_fldr = "figures"
    image_name = os.path.join(main_fldr, create_image_name(image_name))

    # Adding labels and title
    plt.xticks(rotation=90, fontsize=5)
    plt.tight_layout()
    plt.xlabel('Categories')
    plt.ylabel('Values')
    plt.title(title)
    plt.savefig(image_name, format='png')
    plt.show()
    print()

# Given a list and anndata, it finds out the common genes in the list and anndata and writes them in a .csv file.
def find_common_items_index(adata, other_list):

    # returns an array
    adata_list = adata.var.index.values
    common_items_with_indices = [(item, idx1, idx2) for idx1, item in enumerate(adata_list) for idx2, value in
                                 enumerate(other_list) if item == value]
    print(common_items_with_indices)

    # Specify the CSV file name
    csv_filename = create_image_name(name = "common_adata_targetGene", format=".csv")
    csv_filename = os.path.join("../files", csv_filename)

    # Write the data to the CSV file
    with open(csv_filename, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)

        # Write header
        csv_writer.writerow(["Common_gene_name", "Adata_gene_index", "Scimilarity_target_gene_index"])

        # Write data
        for item, idx1, idx2 in common_items_with_indices:
            csv_writer.writerow([item, idx1, idx2])

    print(f"Common items with indices have been written to '{csv_filename}'.")

# Given a csv file and name of a column it reads that column
def read_csv_column(file_path, desired_column):

    # Read the CSV file into a DataFrame
    df = pd.read_csv(file_path)

    # Access the desired column
    column_data = df[desired_column]

    # Print or work with the column data as needed
    return(column_data)

# Given a csv file and name of a column it reads that column
def read_csv_column_perRow(file_path, clm1, clm2, relevant_cell_types ):

    df = pd.read_csv(file_path)

    # Create an empty dictionary to store results
    results = {}

    # Iterate through relevant cell types and extract data
    for cell_type in relevant_cell_types:
        relevant_data = df[df[clm1] == cell_type]
        if not relevant_data.empty:
            # Extract the relevant column value (e.g., 'ColumnName') here
            column_value = relevant_data[clm2].values[0]
            results[cell_type] = column_value

    # Calculate the total number of cells in the relevant cells
    total_cells_in_relevant_data = sum(results.values())

    # Calculate the percentage for each cell type within the relevant cells
    for cell_type, column_value in results.items():
        if total_cells_in_relevant_data > 0:
            cell_percentage = (column_value / total_cells_in_relevant_data) * 1
        else:
            cell_percentage = 0.0  # Set to 0 if no cells in relevant_data

        # results[cell_type] = {
        #     'cell_num': column_value,
        #     'prior': cell_percentage
        # }
        results[cell_type] = cell_percentage


    return results


def make_frequency_matrix(celltype_hits_series, dist=0):
    data_dict = {}  # Dictionary to store all data

    # Iterate through each row and populate the data_dict
    for index, (idx, celltype_str) in enumerate(celltype_hits_series.items()):
        celltype_dict = json.loads(celltype_str)

        for cell, freq in celltype_dict.items():
            if cell not in data_dict:
                data_dict[cell] = [0] * len(celltype_hits_series)

            if dist:
                data_dict[cell][index] = (freq)
            else:
                data_dict[cell][index] = int(freq)

    # Convert the data_dict to a DataFrame
    frequency_df = pd.DataFrame(data_dict, index=celltype_hits_series.index)

    # This might not be necessary now, but just in case
    frequency_df = frequency_df.fillna(0)

    return frequency_df


def violin_plot_figures(cluster1_data, cluster2_data,
                    x_lable1, x_lable2, y_lable):
    plt.figure(figsize=(16, 6))

    # Plot for Cluster 1
    plt.subplot(1, 2, 1)  # (number of rows, number of columns, index of the current plot)
    sns.violinplot(data=cluster1_data, color='skyblue')
    plt.xticks(rotation=90)
    plt.xlabel(x_lable1)
    plt.ylabel(y_lable)
    plt.grid()
    #plt.title()

    # Plot for Cluster 2
    plt.subplot(1, 2, 2)
    sns.violinplot(data=cluster2_data, color='lightgreen')
    plt.xticks(rotation=90)
    plt.xlabel(x_lable2)
    plt.ylabel(y_lable)
    plt.grid()
    #plt.title("Violin plot for Cluster 2")

    plt.tight_layout()  # Ensures that the plots do not overlap
    plt.savefig(('../figures/' + create_image_name(y_lable, format='.png')), bbox_inches='tight', dpi=150)
    plt.show()