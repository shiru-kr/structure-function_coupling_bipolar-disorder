import numpy as np
import pandas as pd 
from scipy.linalg import expm
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from netneurotools.metrics import distance_wei_floyd
import os


###### ADJUST PARAM #####

wd = os.getcwd()
wd = (wd[:-10]+'/')
path = wd+'216/' # adjuest path to be where top wd is for both storing and loading data
FC = 'fc' # pfc/fc - name of .npy for functional connectivity matrices with pearson or partial correlation
SC = 'fbc' # fbc - name of .npy for structural connectivity matrices with fiber-bundle capacity
scope = 'whole' # whole/within/between/cortical

### Load data
hc_sc_mats = np.load(f'{path}hc_{SC}_mats.npy') #load a 3D numpyarray with structural matrices (2D) for all control participant
bd_sc_mats = np.load(f'{path}bd_{SC}_mats.npy') #load a 3D numpyarray with structural matrices (2D) for all bipolar participant
hc_fc_mats = np.load(f'{path}hc_{FC}_mats.npy') #load a 3D numpyarray with functional matrices (2D) for all control participant
bd_fc_mats = np.load(f'{path}bd_{FC}_mats.npy') #load a 3D numpyarray with functional matrices (2D) for all bipolar participant
euc = np.load(f'{path}euc.npy') #load a numpyarray with Euclidean distance
leb = pd.read_csv(f'{path}atlas_labels.csv',header=None) # load labels for parcallations to be used for scope
files = SC+'_'+FC # for saving

## Graph measures functions to be used as structural connectivity predictors:

def communicability(matrix):
    """
    This function takes a 2D matrix and return its communicability matrix
    :param matrix: a matrix of the connection strength between each pair of nodes
    :return: communicability matrix
    """
    mati = np.array(matrix, dtype=float) # make sure its float numpy
    norm = np.linalg.eigvals(mati).max() # finding the largest eigenvalue of the matrix
    nmati = mati / norm # normalise by them to avoid issues with very large values and instability
    return expm(nmati) # return the matrix exponential (communicability)


def compute_weighted_mfpt(W):
    """
    This function takes a weighted adjacency matrix and returns the mean first passage time (MFPT) matrix.
    :param W: 2D weighted adjacency matrix
    :return: MFPT matrix (min-max normalized)
    """

    W = np.array(W, dtype=float)  # Ensure W is a float numpy array

    # Transition matrix P
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # Prevent divide-by-zero
    P = W / row_sums  # normalize by nodal strength
    
    # Stationary distribution pi
    pi = row_sums.flatten() / W.sum()

    # Fundamental matrix Z
    n = W.shape[0]
    I = np.eye(n)
    ones = np.ones((n, 1))
    Z = np.linalg.inv(I - P + ones @ pi.reshape(1, -1))

    # MFPT calculation
    MFPT = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            MFPT[i, j] = (Z[j, j] - Z[i, j]) / pi[j]

    # Min-max normalization: For each column, normalize the values between 0 and 1
    MFPT_minmax = (MFPT - MFPT.min(axis=0)) / (MFPT.max(axis=0) - MFPT.min(axis=0))
    return MFPT_minmax

def short_path(matrix):
    """
    This function takes a matrix and returns its shortest path length matrix.
    :param matrix: 2D matrix
    :return: shortest path matrix
    """

    mati = np.array(matrix, copy=True) #create a copy of the matrix
    mati[mati == 0] = 1e-10 #replace zeros with a tiny number in order to not mess up the inverse
    rec_mat = 1 / mati #compute the inverse of the matrix for weights (path length in terms of cost)
    np.fill_diagonal(rec_mat, 0) #fill diagonal with zeros
    short_len, pre = distance_wei_floyd(rec_mat) #compute shortest paths using Floyd-Warshall algorithm
    return short_len   


## Function to get nodal predictor vectors
def get_predictor_vectors(mats, node, scope):
    """this function takes a list of matrices and a number and return for each matrix a 1D coulmn with the values for that number index.
    depands of the scope, the vector will contain whole brain connections (to all regions), within network connections (to regions in the same network),
    between network connections (to all regions that are not in the same network), or cortical connections (only to cortical regions). 
    Important! labels here match the schafer and tian stlases, adjust if needed.
    #there's probably a cleaner way to do this#
    :param mats: list of structural matrices (same shape)
    :param node: number to use as index
    :param scope: 'between', 'within', 'whole', 'cortical'
    :return: a vector for the node number (excluding self connection)
    """
    selected = leb.iloc[node,0] # get the label of the region (the node)
    nets = ["Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default", "-"] # for subcortex I am using "-" as it's not in the schaefer labels and only apears in all the Tian
    
    if scope == 'between':
        net_indices = [node] # node number as list to be used index and later add to it more indices, if between, add self connection to be removed
    else:
        net_indices = []
    
    for net in nets: # loop over the labels of the networks
        if net in selected: 
            for i in range(len(leb)): 
                check = leb.iloc[i,0] # if the network in the loop matched the label of the node, loop over the number of regions, and pull the region label one by one
                if net in check and i != node: 
                    net_indices.append(i) # if that network is also in the label of the region being checked, and it's not the same region (self connection), add that region (its index) to the list of indices to be used for analysis scope
    net_indices = np.unique(net_indices) # avoid duplicates
    
    result = [] # store the vectors
    for mat in mats: # loop over the matrices
        row_vector = mat[node, :] # extract the column of the region by node number as index
        if scope == "between":
            row_vector = np.delete(row_vector, net_indices) # if between network, remove from the vector the indices of regions in the same network
            result.append(row_vector) # add the list of vectors (for all mats)
        elif scope == "within":
            row_vector = row_vector[net_indices] # if within network, keep only indices of regions within that network
            result.append(row_vector)
        elif scope == "whole":
            row_vector = np.delete(row_vector, node) # if whole brain, use all regions, delete self connection
            result.append(row_vector)
        elif scope == "cortical":
            if node < 200:
                row_vector = row_vector[:200] # if cortical, use indices until 200 as these are the Schaefer one
                row_vector = np.delete(row_vector, node) # delete self connections if the node is in cortical range (discard later anyway if it's over 200)
            elif node >= 200:
                row_vector = row_vector[200:] # if node is over 200, use only subcortical regions
                row_vector = np.delete(row_vector, node-200) # delete self connections (node-200 because the vector is now smaller)
            result.append(row_vector)
        else:
            print("scope does not match") 
    return result

## Function to get nodal vectors to predict
def get_predicted_vector(mat, node, scope):
    """this function takes a matrix and a number and return a 1D coulmn with the values for that number index
    :param mat: a functional matrix
    :param node: number to use as index
    :param scope: 'between', 'within', 'whole', 'cortical'
    :return: a vector for the node number (excluding self connection)
    """
    selected = leb.iloc[node,0]
    nets = ["Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default", "-"] # for subcortex I am using "-" as it's not in the schaefer labels and only apears in all the Tian
    
    if scope == 'between':
        net_indices = [node] # to be used index and later add to it more indices, if between, add self connection to be removed
    else:
        net_indices = []    

    for net in nets:
        if net in selected:
            for i in range(len(leb)):
                check = leb.iloc[i,0]
                if net in check and i != node:
                    net_indices.append(i)
    net_indices = np.unique(net_indices)
    
    row_vector = mat[node, :] # Extract the row for the node
    if scope == "between":
        row_vector = np.delete(row_vector, net_indices)
    elif scope == "within":
        row_vector = row_vector[net_indices]
    elif scope == "whole":
        row_vector = np.delete(row_vector, node)
    elif scope == "cortical":
        if node < 200:
            row_vector = row_vector[:200]
            row_vector = np.delete(row_vector, node)
        elif node >= 200:
            row_vector = row_vector[200:]
            row_vector = np.delete(row_vector, node-200)           
    else:
        raise ValueError(f"Invalid scope: {scope}. Must be 'between', 'within', 'whole', or 'cortical'.")
    return row_vector

## function to precits FC from SC for node
def linear_reg(predictors, predicted):
    """this function takes a list of structural connectivity vectors and a functional connectivity vector and performs a linrear regression to prdict functional from strctural.
    it returns the squered R between the real functional vector and the predicted functional vector.
    :param predictors: list of structural connectivity vectors
    :param functiobal: functional connectivity vector
    :return: a predicted functional matrix"""

    model = LinearRegression() 
    predictors_arr = np.transpose(np.array(predictors)) # dimensions to fit
    scaler = StandardScaler()
    predictors_scaled = scaler.fit_transform(predictors_arr) # standardise the predictors
        
    model.fit(predictors_scaled, predicted) # use linear regression model to predict functional from structural        
    r_squared = model.score(predictors_scaled, predicted) # get the R2 between observed and predicted   

    N = len(predicted) # no. of data points (regions)
    p = predictors_scaled.shape[1] # no. of regressors (IV)
    adjusted_r_squared = 1 - (((1 - r_squared) * (N - 1)) / (N - p - 1))
    return adjusted_r_squared, r_squared

## function to run all other functions for computing regional structure-function coupling 
def couple(FC_mats, SC_mats, euc, SC_measures, scope):
    """
    this function calculates regional SC-FC coupling in the form of R2 between predicted regional FC (from SC) and observed FC. 
    it takes SC matrix and euclidean distance matrix and uses the selected connectivity measures (functions) to compute additional SC predictors. 
    then it predicts the FC for each regions.it can provide R2 values in differenet scopes. "whole-brain": from each region (i) to all other regions (j/=i).
    "within-network": from each region to regions of the same network (labeled according to Schaefer 7 networks). "between-network": from each region to all other 
    regions that are in other networks. "cortical": from region to all other cotrical regions, excluding subcortical regions. 
    #all input matrices should be the same shape#
    :param FC_mats: 2D FC adjacanacy matrices
    :param SC_mats: 2D SC adjacanacy matrices
    :param euc: 2D matrix of euclidean distances
    :param SC_measures: a list of function names. these functions are designed to calculate the network measueres from the SC matrix.
    :param scope: run values for "whole", "cortical", "subcortical", "within" or "between" network.
    :return: the adjusted R2 values in a list. each represents the SC-FC coupling in a region
    """

    mats = [measure(SC_mats) for measure in SC_measures] # run each function in the list to compute the structural connectivity measures
    mats.append(euc) # append euclidean matrix to the list of predictors
    n_nodes = np.shape(FC_mats)[1] # the number of brain regions
    adj_r2_values = [] # empty list to store the adjusted R2 values

    # loop over regions to get the R2 values
    for node in range(n_nodes): #loop over brain regions
        predictors = get_predictor_vectors(mats, node, scope) # use a function to create a list of coulmn vectors to be used as predictors (SC)
        func_vec = get_predicted_vector(FC_mats, node, scope) # use a function to create a coulmn vectors to be used as predicted (FC)

        adjusted_r_squared, r_squared = linear_reg(predictors, func_vec) # use linear regression function to get the adjusted R2 value
        adj_r2_values.append(adjusted_r_squared) # Append the adjusted R2 value for the current region
    
    return adj_r2_values


#loop over participants of each group to calculate regional coupling values and store in a dictionary
hc_group_adj_r2 = {} #store in dictionary
for i in range(hc_sc_mats.shape[2]): #loop over participants
    adj_r2_values = couple(hc_fc_mats[:,:,i], hc_sc_mats[:,:,i], euc,
           SC_measures=(short_path, compute_weighted_mfpt, communicability), scope=scope) #index the 3rd dim, that pulls 2D matrix for each participant one by one
    hc_group_adj_r2[i] = adj_r2_values #add in the dict for each participant the R2 values (=the number of regions)

bd_group_adj_r2 = {}
for i in range(bd_sc_mats.shape[2]): 
    adj_r2_values = couple(bd_fc_mats[:,:,i], bd_sc_mats[:,:,i], euc,
           SC_measures=(short_path, compute_weighted_mfpt, communicability), scope=scope)
    bd_group_adj_r2[i] = adj_r2_values


# create DFs for each group with columns for participants and rows for nodes
hc_node_adj_r2 = [] #empty list to store regional R2 for all participants
for node in range(np.shape(hc_sc_mats[:,:,0])[1]): #loop over number of regions
    row = [] #
    for ind, values in hc_group_adj_r2.items(): # loop R2 dictionary, keys are participants and items are rgonal R2
        row.append(values[node]) #append as a row the regional R2 
    hc_node_adj_r2.append(row)  #append the participant regional coupling values together, so each row is a particiapnt
hc_adj_r2 = pd.DataFrame(hc_node_adj_r2) #convert to DF


bd_node_adj_r2 = []
for node in range(np.shape(hc_sc_mats[:,:,0])[1]):
    row = []
    for ind, values in bd_group_adj_r2.items():
        row.append(values[node])
    bd_node_adj_r2.append(row)  
bd_adj_r2= pd.DataFrame(bd_node_adj_r2)

#save as csv
hc_adj_r2.to_csv(f'{path}results/{files}/{scope}/hc_{files}.csv', index=False)
bd_adj_r2.to_csv(f'{path}results/{files}/{scope}/bd_{files}.csv', index=False)

