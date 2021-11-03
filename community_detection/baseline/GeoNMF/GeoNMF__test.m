%% Clear all things
clc; clear; close all; path(pathdef);

path_to_file = 'data/coauthorship/DBLP1_adjacency.txt';
A = read_data(path_to_file);

path_to_file = './data/coauthorship/DBLP1_community.txt';
theta = read_ground_truth(path_to_file);
[theta_hat, B] = GeoNMF(A, 6);

