import json
import numpy as np
import matplotlib.pyplot as plt
 
# Opening JSON file
f = open('build/results.txt')
 
save_dir = "./python_plots/" 
 
# returns JSON object as
# a dictionary
data = json.load(f)
 
# Iterating through the json
# list

#times = []

for n in [1000000]:
    processors = [2,128,384,768, 2016, 3032]
    data_points = []
    for p in processors:
            #np.amin([data_point['total time'] for data_point in data['data'] if data_point['num_local_vertices'] == n and data_point['p'] == p])
        best_time = np.amin([data_point['total time'] for data_point in data['data'] if data_point['num_local_vertices'] == n and data_point['p'] == p])
        data_points.append([data_point for data_point in data['data'] if data_point['num_local_vertices'] == n and data_point['p'] == p and data_point['total time'] == best_time][0])
    
    time_matrix = []
    steps = ["graph_umdrehen","ruler_pakete_senden","pakete_verfolgen","rekursion_vorbereiten","rekursion","finalen_ranks_berechnen"]

    for step in steps:
        time_matrix.append([np.amax(data_point[step]) for data_point in data_points])
        
    times = [0 for i in range(len(processors))]
    prefix_sum_times = [times];
    for i in range(len(steps)):
        plt.bar(processors, time_matrix[i], bottom=times, label=steps[i])
        times = np.add(times,time_matrix[i])
        prefix_sum_times.append(times)
        
        
        #plt.plot(processors,times,label=steps[i])
   
    
   
    plt.xlabel("processors")
    plt.ylabel("time in ms")
    plt.title("time of different steps in tree rooting algorithm")
    plt.legend()
    
    print("times")
    print([data_point['total time'] for data_point in data_points])
    print("dist_rulers")
    print([data_point['dist_rulers'] for data_point in data_points])

    plt.savefig(save_dir + "wood.pdf")
    plt.clf()
    
    
    communication_times = [data_point["communication"] for data_point in data_points]
    computation_times = [data_point["local_work"] for data_point in data_points]
    
    plt.plot(processors, communication_times ,label="communication")
    plt.plot(processors, computation_times, label="computation")
    plt.xlabel("processors")
    plt.ylabel("time in ms")
    plt.title("time divided in categories")
    plt.legend()
    plt.savefig(save_dir + "unnötig_wood.pdf")
    
    
    
   

   
        
#    for data_point in data['data']:
#        times.append(np.amax(data_point[step]))
#    print(times)
 


 
# Closing file
f.close()