import numpy as np
import pandas as pd
import os

wd = os.getcwd()
wd = (wd[:-10]+'/')
path = wd+'216/' # adjuest path to be where top wd is for both storing and loading data
FC = 'fc'
SC = 'fbc'
PFC = 'pfc'

def strength(matrix):
    strengths = []
    for node in range(matrix.shape[0]):
        srow = np.delete(matrix[node], node)
        strength = np.sum(np.abs(srow))
        strengths.append(strength)
    return strengths

def degree(adjmat):
    adjmat = (adjmat > 0).astype(int)
    degrees = []
    for node in range(adjmat.shape[0]):
        srow = np.delete(adjmat[node], node)
        degree = np.sum(srow)
        degrees.append(degree)
    return degrees

def calculate_metric_for_all(matrices, metric_function):
    metric_results = {}
    for i in range(matrices.shape[2]):
        metric = metric_function(matrices[:, :, i])
        metric_results[i] = metric
    return metric_results

def create_node_participant_matrix(metric_results, num_nodes):
    node_participant_matrix = []
    for node in range(num_nodes):
        row = []
        for ind, values in metric_results.items():
            row.append(values[node])
        node_participant_matrix.append(row)
    return node_participant_matrix

# Load matrices
hc_sc_mats = np.load(path + 'hc_' + SC + '_mats.npy') 
bd_sc_mats = np.load(path + 'bd_' + SC + '_mats.npy')  
hc_fc_mats = np.load(path + 'hc_' + FC + '_mats.npy')  
bd_fc_mats = np.load(path + 'bd_' + FC + '_mats.npy')
hc_pfc_mats = np.load(path + 'hc_' + PFC + '_mats.npy')  
bd_pfc_mats = np.load(path + 'bd_' + PFC + '_mats.npy')

# Calculate strengths and degrees for structural and functional matrices
hc_stre_struc = calculate_metric_for_all(hc_sc_mats, strength)
bd_stre_struc = calculate_metric_for_all(bd_sc_mats, strength)
hc_stre_func = calculate_metric_for_all(hc_fc_mats, strength)
bd_stre_func = calculate_metric_for_all(bd_fc_mats, strength)
hc_stre_pfunc = calculate_metric_for_all(hc_pfc_mats, strength)
bd_stre_pfunc = calculate_metric_for_all(bd_pfc_mats, strength)

hc_deg_struc = calculate_metric_for_all(hc_sc_mats, degree)
bd_deg_struc = calculate_metric_for_all(bd_sc_mats, degree)
#hc_deg_func = calculate_metric_for_all(hc_fc_mats, degree)
#bd_deg_func = calculate_metric_for_all(bd_fc_mats, degree)
#hc_deg_pfunc = calculate_metric_for_all(hc_pfc_mats, degree)
#bd_deg_pfunc = calculate_metric_for_all(bd_pfc_mats, degree)

# Create node-participant matrices for strengths
num_nodes_sc = hc_sc_mats.shape[1]
num_nodes_fc = hc_fc_mats.shape[1]

hc_struc_strength = create_node_participant_matrix(hc_stre_struc, num_nodes_sc)
bd_struc_strength = create_node_participant_matrix(bd_stre_struc, num_nodes_sc)
hc_func_strength = create_node_participant_matrix(hc_stre_func, num_nodes_fc)
bd_func_strength = create_node_participant_matrix(bd_stre_func, num_nodes_fc)
hc_pfunc_strength = create_node_participant_matrix(hc_stre_pfunc, num_nodes_fc)
bd_pfunc_strength = create_node_participant_matrix(bd_stre_pfunc, num_nodes_fc)

# Create node-participant matrices for degrees
hc_struc_degree = create_node_participant_matrix(hc_deg_struc, num_nodes_sc)
bd_struc_degree = create_node_participant_matrix(bd_deg_struc, num_nodes_sc)
#hc_func_degree = create_node_participant_matrix(hc_deg_func, num_nodes_fc)
#bd_func_degree = create_node_participant_matrix(bd_deg_func, num_nodes_fc)
#hc_pfunc_degree = create_node_participant_matrix(hc_deg_pfunc, num_nodes_fc)
#bd_pfunc_degree = create_node_participant_matrix(bd_deg_pfunc, num_nodes_fc)

#Save to CSV
def save_to_csv(data, filename):
    df = pd.DataFrame(data).transpose()
    df.to_csv(filename, index=False)

save_to_csv(hc_struc_strength, path + 'results' +'/strength_degree' +'/hc_strength_'+ SC + '.csv')
save_to_csv(bd_struc_strength, path + 'results' +'/strength_degree' +'/bd_strength_'+ SC + '.csv')
save_to_csv(hc_func_strength, path + 'results' +'/strength_degree' +'/hc_strength_'+ FC + '.csv')
save_to_csv(bd_func_strength, path + 'results' +'/strength_degree' +'/bd_strength_'+ FC + '.csv')
save_to_csv(hc_pfunc_strength, path + 'results' +'/strength_degree' +'/hc_strength_'+ PFC + '.csv')
save_to_csv(bd_pfunc_strength, path + 'results' +'/strength_degree' +'/bd_strength_'+ PFC + '.csv')

save_to_csv(hc_struc_degree, path + 'results' + '/strength_degree' +'/hc_degree_'+ SC + '.csv')
save_to_csv(bd_struc_degree, path + 'results' + '/strength_degree' +'/bd_degree_'+ SC + '.csv')
#save_to_csv(hc_func_degree, path + 'results' + '/strength_degree' +'/hc_degree_'+ FC + '.csv')
#save_to_csv(bd_func_degree, path + 'results' + '/strength_degree' +'/bd_degree_'+ FC + '.csv')
#save_to_csv(hc_pfunc_degree, path + 'results' + '/strength_degree' +'/hc_degree_'+ PFC + '.csv')
#save_to_csv(bd_pfunc_degree, path + 'results' + '/strength_degree' +'/bd_degree_'+ PFC + '.csv')
