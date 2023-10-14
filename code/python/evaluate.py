#very important: this skript handles json objects, the syntax is obviously json, but more specific 
#{"data":[data_point1, data_point2]}


import json
import numpy as np
import matplotlib.pyplot as plt
 

path1 = '../../other/supermuc_auswertung/results_pointer_doubling.txt'
path2 = '../../other/supermuc_auswertung/results_regular_ruling_set2.txt'
f = open(path2)


 

 
save_dir = "" 
 

 


#this function groups the json elements in path by the value "p" and returns the mean total_time
def get_PE_total_time_tuple(path):
    f = open(path)
    data = json.load(f)
   
    processors = [data_point['p'] for data_point in data['data']]
    processors = list(dict.fromkeys(processors))
    processors.sort() #so we have unique sorted processors array
    
    times = []
    for p in processors:
        relevant_data_points = [data_point for data_point in data['data'] if data_point['p'] == p]
        relevant_total_times = [np.amax(data_point['total_time'])  for data_point in relevant_data_points]
        times.append(np.mean(relevant_total_times))
    
    return processors, times


def generate_time_step_graph(path):
    f = open(path)
    data = json.load(f)
   
    processors = [data_point['p'] for data_point in data['data']]
    processors = list(dict.fromkeys(processors))
    processors.sort()
    
    steps = data['data'][0]['time_step_names']
    
    times = [0 for i in range(len(processors))]
    prefix_sum_times = [times];
    
    #first find for every p in processors one data_point to print
    data_points = []
    for p in processors:
        data_points.append([data_point for data_point in data['data'] if data_point['p'] == p][0])
        
    time_matrix = []
    for step in steps:
        pe = 0 #0 <= pe < p
        time_matrix.append([data_point[step][pe] for data_point in data_points])
    
    times = [0 for i in range(len(processors))]
    prefix_sum_times = [times];
    for i in range(len(steps)):
        plt.bar(processors, time_matrix[i], bottom=times, label=steps[i])
        times = np.add(times,time_matrix[i])
        prefix_sum_times.append(times)
    
    
    plt.xlabel("processors")
    plt.ylabel("time in ms")
    plt.title("time of different steps in tree rooting algorithm")
    plt.legend()
    
 

    plt.savefig(save_dir + "wood.pdf")
    plt.clf()
    
    
    


generate_time_step_graph(path2)


for n in [1000000]:
    processors = [2,4,8,16,32,64,128,256,512]
    
    
    #first plot total_times
    
    times = []
    for p in processors:
        relevant_data_points = [data_point for data_point in data['data'] if data_point['num_local_vertices'] == n and data_point['p'] == p]
        relevant_total_times = [np.amax(data_point['total_time'])  for data_point in relevant_data_points]
        times.append(np.mean(relevant_total_times))
    print(quad(2))
    
    plt.bar(processors, times, label="pointer_doubling")
    plt.xlabel("processors")
    plt.ylabel("time in ms")
    plt.title("time of different steps in tree rooting algorithm")
    plt.legend()
    plt.savefig("pointer_doubling.pdf")

    
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