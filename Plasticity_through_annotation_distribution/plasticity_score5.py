
import os
import numpy as np
import pandas as pd
np.float = float
import scanpy as sc

import utils_AT
import utils_AT as utils
from scipy.stats import multivariate_normal #is this scimilarity? 
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns

# Term1: Calculating Entropy for Each Cell
def calculate_entropy(cell_probabilities): #is this a 2-d array? 
    """Calculate entropy for a given probability distribution."""

    entropies = []

    for cell in cell_probabilities:
        # Ensure no zero probabilities and convert to numpy array
        probabilities = np.array([p + 1e-10 for p in cell]) #based off this 2d array? 

        # Calculate entropy
        entropy = -np.sum(probabilities * np.log(probabilities)) #sum the valuea and multiply by its log - assuming some entropy function? 
        entropies.append(entropy)
    #print("printing entropies")
    #print(entropies)
    #print("printing entropies: ", entropies)
    #okay so for each cell_type there is an association 
    return entropies #return entropy for each cell -> type ->list


# Term2: Log Density and Distance Multiplication
#Term2 Previous Version 
def calculate_term2(dist_matrix, log_densities, cell_types, query_cell):
    #what exactly is each input type? 
    term2_values = []
    #print(len(query_cell))
    #print("distance matrix dimensions",len(dist_matrix),len(dist_matrix[0]))
    #print("printing distance matrix: ", dist_matrix)
    for i in range(len(query_cell)): #what holds into a query cell 
        weighted_dist = 0
        #print(log_densities) #this should be an array, but is a scalare 
        positive_log_densities = {ct: log_densities[ct][i] for ct in cell_types if log_densities[ct][i] > 0}
        for ct1, ct2 in combinations(positive_log_densities.keys(), 2):
            distance = dist_matrix[cell_types.index(ct1)][cell_types.index(ct2)]
            weighted_dist += distance * positive_log_densities[ct1] * positive_log_densities[ct2]

        term2_values.append(weighted_dist)
        #print("term2 cell no: ", i)
    return term2_values

#Term2 Version 1
def calculate_term2_updated(dist_matrix,cell_types,query_cell, freq_mat):
    #for a given query cell 
    #Sprint("query cell: ", query_cell)
    #print("dist_matrix: ", dist_matrix)
    #print("freq mat: ", freq_mat)
    #print("printing cell_types: ", cell_types)
    term2_values = []
    #goal, get the freq mat with the keys and the row value where row_val > 0

    print("printing freq mat: ",freq_mat)
    #get 1st row where vals > 0 
    #print("printing first row of freq mat", freq_mat.iloc[0]) #-> this actually did get the first val w labels 
    freq_mat_row = freq_mat.iloc[0]
    print("printing first row of freq mat",freq_mat_row)
    #now remove the zero_value 
    #for x in freq_mat_row.keys():
    freq_mat_row_pruned = freq_mat_row.loc[freq_mat_row > 0]    
    print("printing pruned freq mat",freq_mat_row_pruned)
    print("new cell types to consider", freq_mat_row_pruned.keys().tolist())

    freq = freq_mat.values.tolist()
    print("number of cells checking: ", len(freq[0]))
    #print("length of query cell is: ", len(query_cell))
    #print("distance matrix dimensions",len(dist_matrix),len(dist_matrix[0]))
    #print("printing distance matrix: ", dist_matrix)
    for i in range(len(query_cell)): #what holds into a query cell 
        #positive_log_densities = {ct: log_densities[ct][i] for ct in cell_types if log_densities[ct][i] > 0}
        if i%50 == 0:
            print("we are on iteration: ", i)
        weighted_dist = 0
        #freq_mat = freq_mat

        #breaking it down to a row work:
        desired_row = freq_mat.iloc[i] #gets row
        pruned_row = desired_row.loc[desired_row > 0] #removes values of 0 
        weight_cell_types = pruned_row.keys().tolist()
        #freq_mat.columns.tolist()
        current_cell = i
        current_freq = freq[current_cell]
        print("printing current freq table: ", current_freq)
        #for ct1, ct2 in combinations(freq_mat.keys(), 2):
        for ct1, ct2 in combinations(weight_cell_types, 2):
            distance = dist_matrix[cell_types.index(ct1)][cell_types.index(ct2)]
            #print("printing current cell type ", current_cell)
            current_cell = int(current_cell)
            weight1 = freq[current_cell][weight_cell_types.index(ct1)] #indexing into weight cells cause need the specific ones 
            weight2 = freq[current_cell][weight_cell_types.index(ct2)]
            #print(ct1, ct2)
            #print(cell_types.index(ct1), " ",cell_types.index(ct2), " ")
            #print(weight1)
            #print(weight2)
            weight = (weight1 + weight2)/2
            #print('weight: ',weight)
            #print('distance: ', distance)
            weighted_dist += distance * weight
        #print("totoal weighted is: ", weighted_dist)
        term2_values.append(weighted_dist)
        #print("term2 cell no: ", i)
    return term2_values



    #return query_cell

#Term2 Version 2
def calculate_term2_speed_updated(dist_matrix,cell_types,query_cell, freq_mat):
    term2_values = []
    freq = freq_mat.values.tolist()
    print("number of cells checking: ", len(freq[0])) #just for seeing progress purposes 
    for i in range(len(query_cell)): #what holds into a query cell 
        if i%100 == 0:
            print("we are on iteration: ", i)
        weighted_dist = 0
        #breaking it down to a row work:
        desired_row = freq_mat.iloc[i] #gets row
        pruned_row = desired_row.loc[desired_row > 0] #removes values of 0 
        weight_cell_types = pruned_row.keys().tolist()
        final_freq = pruned_row.tolist()
        #for ct1, ct2 in combinations(freq_mat.keys(), 2):
        for ct1, ct2 in combinations(weight_cell_types, 2):
            distance = dist_matrix[cell_types.index(ct1)][cell_types.index(ct2)]    
            weight1 = final_freq[weight_cell_types.index(ct1)] #indexing into weight cells cause need the specific ones 
            weight2 = final_freq[weight_cell_types.index(ct2)]
            weight = (weight1 + weight2)/2
            weighted_dist += distance * weight
        if weighted_dist < 0:
            print("happening: ", weighted_dist)
        term2_values.append(weighted_dist)
        #print("term2 cell no: ", i)
    return term2_values

#Term2 Version 3 
def calculate_term2_speed_averages(dist_matrix,cell_types,query_cell, freq_mat):
    term2_values = []
    freq = freq_mat.values.tolist()
    outlier_10 = []
    norm_10 = []
    outlier_count = 0 
    norm_count = 0 

    query_cell_names = freq_mat.index.tolist()
    for i in range(len(query_cell)): #what holds into a query cell 
        if i%100 == 0:
            print("we are on iteration: ", i)
        weighted_dist = 0
        desired_row = freq_mat.iloc[i] #gets row
        pruned_row = desired_row.loc[desired_row > 0] #removes values of 0 
        weight_cell_types = pruned_row.keys().tolist()
        final_freq = pruned_row.tolist()
        counter = 0
        combos = combinations(weight_cell_types, 2)
        for ct1, ct2 in combinations(weight_cell_types, 2):
            distance = dist_matrix[cell_types.index(ct1)][cell_types.index(ct2)]
            weight1 = final_freq[weight_cell_types.index(ct1)] #indexing into weight cells cause need the specific ones 
            weight2 = final_freq[weight_cell_types.index(ct2)]
            weight = (weight1 + weight2)/2
            weighted_dist += distance * weight
            counter += 1 
        #print("length of counter is: ", counter)
        try:
            avg_weighted_distance = weighted_dist/counter
        except: 
            avg_weighted_distance = 0 
        term2_values.append(avg_weighted_distance)
        if norm_count < 10:
            if avg_weighted_distance >= 1 and avg_weighted_distance <=3:
                norm_10.append(query_cell_names[i])
                norm_count += 1
        if  outlier_count < 10:
            if avg_weighted_distance >= 9 and avg_weighted_distance <=11:
                outlier_10.append(query_cell_names[i])
                outlier_count += 1

        #print("term2 cell no: ", i)
    print("printing normal set: ", norm_10)
    print("printing outlier set: ", outlier_10)
    outlier_write(norm_10,outlier_10,freq_mat)
    return term2_values

#Term2 Version 4
def calculate_term2_speed_multiplcation(dist_matrix,cell_types,query_cell, freq_mat):
    term2_values = []
    freq = freq_mat.values.tolist()
    print("number of cells checking: ", len(freq[0])) #just for seeing progress purposes 
    for i in range(len(query_cell)): #what holds into a query cell 
        if i%100 == 0:
            print("we are on iteration: ", i)
        weighted_dist = 0
        #freq_mat = freq_mat
        #breaking it down to a row work:
        desired_row = freq_mat.iloc[i] #gets row
        pruned_row = desired_row.loc[desired_row > 0] #removes values of 0 
        weight_cell_types = pruned_row.keys().tolist()
        final_freq = pruned_row.tolist()
        #for ct1, ct2 in combinations(freq_mat.keys(), 2):
        for ct1, ct2 in combinations(weight_cell_types, 2):
            distance = dist_matrix[cell_types.index(ct1)][cell_types.index(ct2)]
            
            weight1 = final_freq[weight_cell_types.index(ct1)] #indexing into weight cells cause need the specific ones 
            weight2 = final_freq[weight_cell_types.index(ct2)]
            weight = weight1*weight2
            weighted_dist += distance * weight
        if weighted_dist < 0:
            print("happening: ", weighted_dist)
        term2_values.append(weighted_dist)
        #print("term2 cell no: ", i)
    return term2_values




#this was just used for a little bit of analysis on outlier values ---- can be ignored unless needed again 
def outlier_write(norm_10,outlier_10, freq_mat):
    #given names in each, get dataframes of each 
    #print("printing freq mat: ", freq_mat)
    outlier_frame = freq_mat.loc[outlier_10]
    norm_frame =  freq_mat.loc[norm_10]
    
    #now get the average number where the value is greater than 0 
    #first outliers 
    o_count = len(outlier_frame)
    o_non_zero_count = 0
    for row_index in range(o_count):
            row = outlier_frame.iloc[row_index] #gets all values in that row
            for cell in row:
                if cell > 0:
                    o_non_zero_count += 1 
    try:
        o_avg_length = o_non_zero_count/o_count
    except: 
        o_avg_length = 0
    print("length of outlier frame: ", o_count)
    #second normal 
    n_count = len(norm_frame)
    print("length of normal frame: ", n_count)
    n_non_zero_count = 0
    for row_index in range(n_count):
            row = norm_frame.iloc[row_index] #gets all values in that row
            for cell in row:
                if cell > 0:
                    n_non_zero_count += 1 
    try:
        n_avg_length = n_non_zero_count/n_count
    except:
        n_avg_length = 0

    
    print("outlier data frame: ", o_avg_length)
    print("normal data frame: ", n_avg_length)

#this was an attempt to use numpy to speed up calculations
#however, bottleneck is 
def calculate_term2_numpy(dist_matrix,cell_types,query_cell, freq_mat):
    term2_values = []
    freq = freq_mat.values.tolist()
    print("number of cells checking: ", len(freq[0])) #just for seeing progress purposes 
    for i in range(len(query_cell)): #what holds into a query cell 
        if i%100 == 0:
            print("we are on iteration: ", i)
        weighted_dist = 0
        #freq_mat = freq_mat

        #breaking it down to a row work:
        desired_row = freq_mat.iloc[i] #gets row
        pruned_row = desired_row.loc[desired_row > 0] #removes values of 0 
        weight_cell_types = pruned_row.keys().tolist()
        final_freq = pruned_row.tolist()
        #for ct1, ct2 in combinations(freq_mat.keys(), 2):
        combos = combinations(weight_cell_types, 2)
        weight_vec = []
        distance_vec = []
        for i1,i2 in combos:
            weight_vec.append((final_freq[weight_cell_types.index(i1)]+final_freq[weight_cell_types.index(i2)])/2)
            distance_vec.append((dist_matrix[cell_types.index(i1)][cell_types.index(i2)]))
        #print("weight vec: ", weight_vec)
        #print("distance vec: ",distance_vec)
        term2_values.append(np.dot(np.array(weight_vec), np.array(distance_vec)))
        #print("term2 vals: ", list(term2_values))
        #print("term2 cell no: ", i)
    return term2_values





def process_clusters(Y_clusters, S_cluster):
   
    #data_path = 'data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN2_75.h5ad'
    #data_path = 'data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN2.h5ad'
    data_path = os.path.join(os.getcwd(), 'data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN2_150.h5ad')
    gaussian_path = os.path.join(os.getcwd(), 'data/scimilarity_data_division/fitted_gaussian_noGMM')
    bhatt_dist_path = os.path.join(os.getcwd(), 'data/scimilarity_data_division/bhatt_dist_NO_GMM.npy')

    adams_comp = sc.read(data_path)
    epsilon = 1e-10

    # Initialize a dictionary to store results
    plasticity_score_d = {}
    entropy_d = {}
    term2_d = {}

    for Y_cluster in Y_clusters:

        adams1 = adams_comp[adams_comp.obs['cell_type_med'].isin([Y_cluster])]#[:50] #means choosing 5 cells 
        freq_mat = utils.make_frequency_matrix(adams1.obs['sc_hits'])
        print("freq mat")
        #print(freq_mat)
        #print(type(freq_mat))
        freq_mat.columns = freq_mat.columns.str.replace(' ', '_')   #use this to calculate entropy
        #normalize the data - to be probabilities instead of counts 
        #rows = freq_mat.rows
        freq_mat_normalized = freq_mat.copy()
        num_rows = len(freq_mat)
        for row_index in range(num_rows):
            row = freq_mat.iloc[row_index] #gets all values in that row
            #get sum
            total_count = 0 
            for cell in row:
                total_count += cell
            #now normalize 
            for i in range(len(row)):
                freq_mat_normalized.iloc[row_index,i] = float(row[i])/float(total_count)
        #print("printing normalized freq mat")
        #print(freq_mat_normalized)
        freq_mat = freq_mat_normalized 
        print("freq mat normalized")
        #print(freq_mat)
        #for each row just get the total count and divide each by that value 


        
        cell_types_to_consider = freq_mat.columns.tolist()
        #print("cell_types_to_consider", cell_types_to_consider)
        #print(freq_mat)
        # load the fitted gaussian mean, variance an cell numbers
        
        means = {cell_type: np.load(os.path.join(gaussian_path, f'{cell_type}_mean.npy'))
                 for cell_type in cell_types_to_consider}
        covariances = {cell_type: np.load(os.path.join(gaussian_path, f'{cell_type}_cov.npy'))
                       for cell_type in cell_types_to_consider}
        priors = utils.read_csv_column_perRow(
            'files/scimilarity_trainingset_cellnum_AuthorLable_original.csv',
            'Cell_Type', 'Number of Samples', cell_types_to_consider)

        '''
        print(len(cell_types_to_consider))
        print("cell types to consider: ", cell_types_to_consider)
        print("og priors: ",priors) #so i guess these are already normalized 
        '''
        final_results_per_cluster = []  # Store final results for each cell in this cluster

        cell = adams1.obsm['X_scimilarity']

        # 2. Calculate Likelihood
        likelihoods = {}
        for cell_type in cell_types_to_consider:
            pdf_value = multivariate_normal.pdf(cell,
                                                mean=means[cell_type].squeeze(),
                                                cov=covariances[cell_type].squeeze())
            likelihoods[cell_type] = pdf_value + epsilon

        # Evidence for this cell
        evidence = sum(likelihoods[cell_type] * priors[cell_type] for cell_type in cell_types_to_consider)
        
        

        # 4. Calculate Posterior Probability
        posterior_probabilities = {ct: (likelihoods[ct] * priors[ct]) / evidence for ct in
                                   cell_types_to_consider}

        #print("posterior probs: ", posterior_probabilities)

        # fetch information for every cell
        num_cells = len(next(iter(posterior_probabilities.values())))
        # Initialize a dictionary to store cell-wise data
        #cell_wise_data = [[] for _ in range(num_cells)]

        # Iterate over each row and distribute its data to the respective cells
        #for row in posterior_probabilities.values():
        #    for cell_index, cell_value in enumerate(row):
        #        cell_wise_data[cell_index].append(cell_value)
        
        # 5. Calculate Entropy (Term 1)
        
        #print("cell wise data")
        cell_wise_data_simple = freq_mat.values.tolist()
        #this is the problem 
        #evidence = sum(cell_wise_data_simple[cell_type] * priors[cell_type] for cell_type in cell_types_to_consider)
        #freq_probabilities = {ct: (cell_wise_data_simple[ct] * priors[ct]) / evidence for ct in
                                   #cell_types_to_consider}
        
        freq_probabilities = np.zeros((len(cell_wise_data_simple),len(cell_wise_data_simple[0])))
        for r_index,row in enumerate(cell_wise_data_simple):
            for c_index, col in enumerate(row):
                cell_type = cell_types_to_consider[c_index] 
                #print("cell type is: ", cell_type)
                freq_probabilities[r_index][c_index] = cell_wise_data_simple[r_index][c_index]*priors[cell_type]

                # TODO # normalize these values across a query cell -> take old code s
        
        #normalizing data again 
        #print("not normailzed freq input w priors: ", freq_probabilities)
        updated_freq_probs = np.zeros((len(cell_wise_data_simple),len(cell_wise_data_simple[0])))
        for row_index,row in enumerate(freq_probabilities):
            total_count = 0 
            for cell in row:
                total_count += cell
            for i in range(len(row)):
                updated_freq_probs[row_index][i] = float(row[i])/float(total_count)

        #print("final frequency input w priors: ",updated_freq_probs)


        #print("frequency matrix: ",cell_wise_data_simple)
        #print("with priors: ", freq_probabilities)

        
        #print("priors", priors)
        '''
        print("specific cell type: ", cell_type)
        print("prior for specific cell type: ", priors[cell_type])
        print("cell wise data", freq_mat)
        '''
        #incorpoorate priors -> 


        #include priors 
        

        #print(len(cell_wise_data))
        #print(len(cell_wise_data[0]))
        #print(cell_wise_data_simple)
        #print("Cell types to consider: ", cell_types_to_consider)
        #cell_wise_with_priors = for ct in cell_types_to_consider

        cell_wise_prob_simple = [np.array(cell) for cell in updated_freq_probs]
        
        #cell_wise_prob_simple = [np.array(cell) for cell in cell_wise_data_simple]
        entropy_arr = calculate_entropy(cell_wise_prob_simple)
        
        '''likelihoods = {}
        for cell_type in cell_types_to_consider:
            pdf_value = multivariate_normal.pdf(cell)
            likelihoods[cell_type] = pdf_value + epsilon'''

        # 6. Calculate Term 2
        cell = adams1.obsm['X_scimilarity']
        likelihoods = {}
        for cell_type in cell_types_to_consider:
            pdf_value = multivariate_normal.pdf(cell,
                                                mean=means[cell_type].squeeze(),
                                                cov=covariances[cell_type].squeeze())
            likelihoods[cell_type] = pdf_value + epsilon
        log_densities = {ct: np.log(likelihoods[ct] + epsilon) for ct in cell_types_to_consider}
        dist_matrix = np.load(bhatt_dist_path)
        #print("printing dist matrix: ", dist_matrix)
        #build normalized dist_matrix
        '''dist_matrix_normalized = np.zeros((len(dist_matrix),len(dist_matrix)))
        for i,row in enumerate(dist_matrix):
            total = np.sum(row)
            for j, value in enumerate(row):
                dist_matrix_normalized[i][j] = value/total
        #print("printing dist matrix normalized: ", dist_matrix_normalized)
        #print("printing distance matrix: ", dist_matrix)
        print("norm freq matrix: ", dist_matrix)'''

        #print("printing a specific spot in freq mat: ", freq_mat_normalized)
        #print("what is freq_mat?: ", freq_mat_normalized.values.tolist()[0][0])
        #term2_arr = calculate_term2(dist_matrix, log_densities, cell_types_to_consider, cell)
        print("before term 2")
        #term2_arr = calculate_term2_updated(dist_matrix,cell_types_to_consider,cell, freq_mat)
        term2_arr = calculate_term2_speed_multiplcation(dist_matrix,cell_types_to_consider,cell, freq_mat)
        #term2_arr = calculate_term2_numpy(dist_matrix,cell_types_to_consider,cell, freq_mat)
        print("after term 2")
        #term2_arr = calculate_term2(dist_matrix, log_densities, cell_types_to_consider, cell) #apparently cell is an int 
        #print(freq_mat)
        # 7. Calculate Final Result
        plasticity_score_arr = [e * t2 for e, t2 in zip(entropy_arr, term2_arr)]

        # create dictionary for every Y cluster
        entropy_d[Y_cluster] = entropy_arr
        term2_d[Y_cluster] = term2_arr
        plasticity_score_d[Y_cluster] = plasticity_score_arr

    return entropy_d, term2_d, plasticity_score_d


S_cluster = ""
Y_clusters = ["SCLC-A", "T cell"]
entropy_d, term2_d, plasticity_score_d = process_clusters(Y_clusters, S_cluster)
print("stuff generated")


def create_df(term_dict):
    SLC = term_dict['SCLC-A']
    T = term_dict['T cell']   
    #T.append(.625)
    data = {
    'SLC-A' : np.array(SLC),
    'T-cell' : np.array(T)
    }
    maxsize = max([a.size for a in data.values()])
    data_pad = {k:np.pad(v, pad_width=(0,maxsize-v.size,), mode='constant', constant_values=np.nan) for k,v in data.items()}
    df = pd.DataFrame(data_pad)
    return df

def term2_outliers(term_dict):
    T = term_dict['T cell'] 
    outlier_10 = []
    norm_10 = []
    outlier_count = 0 
    norm_count = 0 
    print("printing length of T: ",len(T))
    print("printing T: ",T)
    '''for item in T:
        if T[item] >= 1 and T[item] <=3:
            norm_10.append(item)
            norm_count += 1
        elif T[item] >= 9 and T[item] <=11:
            outlier_10.append(item)
            outlier_count += 1
        if norm_count >= 10 and outlier_count >= 10:
            break'''
    print("printing normal set: ", norm_count)
    print("printing outlier set: ", outlier_count)
        
    

def create_two_charts(term_dict,label):
    SLC = term_dict['SCLC-A']
    T = term_dict['T cell']   
    #T.append(.625)
    sdf = pd.DataFrame.from_dict(SLC)
    tdf = pd.DataFrame.from_dict(T)
    plt.figure(figsize=(16, 10))
    colors = ['skyblue', 'lightgreen', 'violet']  # Define colors for each cluster
    sns.violinplot(data=sdf)
    plt.xlabel("SLC-A")
    plt.ylabel(label)
    plt.title("MSK Lung cancer data")
    plt.grid()
    plt.savefig("figures/"+utils.create_image_name(label+"SLC"))

    plt.figure(figsize=(16, 10))
    colors = ['orange', 'lightgreen', 'violet']  # Define colors for each cluster
    sns.violinplot(data=tdf)
    plt.xlabel("T cell")
    plt.ylabel(label)
    plt.title("MSK Lung cancer data")
    plt.grid()
    plt.savefig("figures/"+utils.create_image_name(label+"T"))
    print("done")
    return "done"

'''pd.DataFrame.from_dict(entropy_d).to_csv('entropy_d.csv', index=False)
pd.DataFrame.from_dict(term2_d).to_csv('term2_d.csv', index=False)
pd.DataFrame.from_dict(plasticity_score_d).to_csv('plasticity_score_d.csv', index=False)'''
'''results = pd.Series(entropy_d)
#results = pd.DataFrame(results_temp)
print("printing entropy")
print(results)
plt.figure(figsize=(16, 10))
colors = ['skyblue', 'lightgreen', 'violet']  # Define colors for each cluster
sns.violinplot(data=results)
plt.xticks(ticks=range(len(results.index)), labels=results.index)
plt.xlabel("Cell clusters")
plt.ylabel("entropy")
plt.title("MSK Lung cancer data")
plt.grid()
plt.savefig("figures/"+utils_AT.create_image_name("entropy"))'''
#plt.savefig("entropy_test.png")
#plt.show()th



#print(term2_d)
'''SLC = term2_d['SCLC-A']
T = term2_d['T cell']
#print(SLC)
print(len(SLC))
T.append(.625)
print(len(T))
data = {
    'SLC-A' : np.array(SLC),
    'T-cell' : np.array(T)
}

maxsize = max([a.size for a in data.values()])
data_pad = {k:np.pad(v, pad_width=(0,maxsize-v.size,), mode='constant', constant_values=np.nan) for k,v in data.items()}
df = pd.DataFrame(data_pad)
print(df)'''
#results = pd.DataFrame(SLC,T).tranpose()
#results = results.set_axis(['SCLC-A','T cell'], axis=1)
#results = pd.DataFrame.from_dict(term2_d)


#junk = create_two_charts(entropy_d,"entropy")
#junk1 = create_two_charts(plasticity_score_d,"plasticity_score")
#junk2 = create_two_charts(term2_d,"term2")


results = create_df(term2_d)
#term2_outliers(term2_d)
#print("printing term2 results")
#print(results)
plt.figure(figsize=(16, 10))
colors = ['skyblue', 'lightgreen', 'violet']  # Define colors for each cluster
sns.violinplot(data=results)
#plt.xticks(ticks=range(len(results.index)), labels=results.index)
plt.xlabel("Cell clusters")
plt.ylabel("term 2")
plt.title("MSK Lung cancer data")
plt.grid()
#plt.savefig("figures/"+utils_AT.create_image_name("term2"))
plt.savefig("figures/"+utils.create_image_name("term2"))
plt.show()

#results = pd.DataFrame.from_dict(entropy_d)
results = create_df(entropy_d)
#print("printing entropy results")
#print(results)
plt.figure(figsize=(16, 10))
colors = ['skyblue', 'lightgreen', 'violet']  # Define colors for each cluster
sns.violinplot(data=results)
#plt.xticks(ticks=range(len(results.index)), labels=results.index)
plt.xlabel("Cell clusters")
plt.ylabel("entropy")
plt.title("MSK Lung cancer data")
plt.grid()
#plt.savefig("figures/"+utils_AT.create_image_name("entropy"))
plt.savefig("figures/"+utils.create_image_name("entropy"))
#plt.savefig("entropy_test.png")
plt.show()

#results = pd.DataFrame.from_dict(plasticity_score_d)
results = create_df(plasticity_score_d)
#print("printing plasticity results")
#print(results)
plt.figure(figsize=(16, 10))
colors = ['skyblue', 'lightgreen', 'violet']  # Define colors for each cluster
sns.violinplot(data=results)
#plt.xticks(ticks=range(len(results.index)), labels=results.index)
plt.xlabel("Cell clusters")
plt.ylabel("plasticity_score")
plt.title("MSK Lung cancer data")
plt.grid()
#plt.savefig("figures/"+utils_AT.create_image_name("plasticty_score"))
plt.savefig("figures/"+utils.create_image_name("plasticty_score"))
#plt.show()
print()


# for i, cluster in enumerate(Y_clusters):
#     sns.violinplot(data=results[cluster], color=colors[i])
# Now `results` contains the final results for each cluster in Y_clusters

# Plot for Cluster 1
#plt.subplot(1, 2, 1)  # (number of rows, number of columns, index of the current plot)
#sns.violinplot(data=results, color='skyblue')
#
#sns.violinplot(x = results.index, y=results, color='skyblue')
#plt.xticks(rotation=90)
