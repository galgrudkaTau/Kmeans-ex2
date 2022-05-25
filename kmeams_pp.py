import numpy as np
import pandas as pd
import sys 
import math
import mykmeanssp


np.random.seed(0)

def invalid_input():
    sys.exit("Invalid Input!")

def an_error_has_occurred():
    sys.exit("An Error Has Occurred")

def files_to_dataframe(file_name_1,file_name_2):
    file1 = pd.read_csv(file_name_1)
    size1=file1.shape[1]
    file2 = pd.read_csv(file_name_1)
    size = [str(i) for i in range(file1.shape[1])]
    file1 = pd.read_csv(file_name_1, names=size)
    size = [str(i+size1-1) for i in range(file2.shape[1])]
    size[0] = str(0)
    file2 = pd.read_csv(file_name_2, names=size)

    data = pd.merge(file1, file2, on = '0')
    data = data.sort_values(by=['0'])
    return data

def kmeanspp(matrix, k):
    matrix_idx = range(len(matrix))
    if (k>len(matrix)):
        invalid_input()
    first_idx = np.random.choice(matrix_idx)
    centroids = []
    init_idx = []
    init_idx.append(first_idx)
    #print("appended first centroid "+"".join(str(first_idx)))
    centroids.append(matrix[first_idx])
    #print("first centroid "+"".join(str(matrix[first_idx])))
    while (len(centroids)<k):
        D = np.full((len(matrix)),float('inf'))
        for l,datapoint in enumerate(matrix):
            dist = [calculate_distance(centroid, datapoint) for j,centroid in enumerate(centroids)]
            D[l] = min(dist)   
        # print(D[:5])
        Dm = sum(D)
        P = [D[i]/Dm for i in matrix_idx]
        idx_chosen = np.random.choice(matrix_idx, p=P)
        init_idx.append(idx_chosen)
        #print(idx_chosen)
        centroids.append(matrix[idx_chosen])
        init_centroids = np.stack(centroids)
    print("Initial centroids to be "+"".join(str(init_idx))) 

    return init_idx

def calculate_distance(centroid, data_point):
    return sum([pow(centroid[i]-data_point[i],2) for i in range(len(centroid))])

def check_is_natural(num):
    try:
        # Convert it into float
        
        val_f = float(num)
        val_int = int(float(num))
    
    except ValueError:
        invalid_input()
    
    if val_f!=val_int or val_int<=0:
        invalid_input()

def check_is_float(num):
    try:
        # Convert it into float
        
        val_f = float(num)

    except ValueError:
        invalid_input()

def main():
    #try:
        check_is_natural(sys.argv[1])
        k=int(sys.argv[1])
        max_iter=300
        if(len(sys.argv)==5):
            check_is_float(sys.argv[2])
            epsilon=sys.argv[2]
            file_name_1=sys.argv[3]
            file_name_2=sys.argv[4]
            
        elif (len(sys.argv)==6):
            check_is_natural(sys.argv[2])
            max_iter=int(sys.argv[2])
            check_is_float(sys.argv[2])
            epsilon=sys.argv[3]
            file_name_1=sys.argv[4]
            file_name_2=sys.argv[5]
        else: 
            invalid_input()
         #k = 3
        #file_name_1 = "input_1_db_1.txt"
        #file_name_2 = "input_1_db_2.txt"
        input_data = files_to_dataframe(file_name_1, file_name_2)
        data = input_data.drop(['0'],axis=1) # data frame 
        input_matrix = data.to_numpy() #nd array 
        d=data.shape[1]
        #centroids = kmeanspp(input_matrix,k)

        idxs = kmeanspp(input_matrix,k) #initialize centroids 
        print(kmeanspp.kmeanssp(max_iter, epsilon, idxs ,input_matrix))
    #except Exception as e:
        # print("An Error Has Occurred\n")
        # exit()
    
main()
