#very important: this skript handles json objects, the syntax is obviously json, but more specific 
#{"data":[data_point1, data_point2]}


import json
import numpy as np
import matplotlib.pyplot as plt
 


path1 = '../../other/supermuc_auswertung/results_regular_ruling_set.txt'
path2 = '../../other/supermuc_auswertung/results_regular_ruling_set_rec.txt'
path3 = '../../other/supermuc_auswertung/results_regular_ruling_set2.txt'
path4 = '../../other/supermuc_auswertung/results_regular_ruling_set2_rec.txt'
path5 = '../../other/supermuc_auswertung/results_euler_tour.txt'
path6 = '../../other/supermuc_auswertung/results_wood_regular_ruling_set.txt'

paths = [path1, path2, path3, path4, path5, path6]
names = ["regular_ruling_set", "regular_ruling_set_rec", "regular_ruling_set2", "regular_ruling_set2_rec", "euler_tour", "wood_regular_ruling_set"]

 

 

def to_json(path):
    with open(path, 'r+') as file: 
        file_data = file.read() 
        file.seek(0, 0) 
        file.write("{\n\"data\":[\n" + file_data[:-2] + "\n}]}") 


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


def generate_time_step_graph(path, algorithm_name):
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
        plt.bar(processors, time_matrix[i], bottom=times, label=steps[i], width = 50)
        times = np.add(times,time_matrix[i])
        prefix_sum_times.append(times)
    
    
    plt.xlabel("processors")
    plt.ylabel("time in ms")
    plt.title("time of different steps in " + algorithm_name + " algorithm")
    plt.legend()
    
 

    plt.savefig("time_step_" + algorithm_name + ".pdf")
    plt.clf()
    

paths = [path1, path2, path3, path4, path5, path6]
names = ["regular_ruling_set", "regular_ruling_set_rec", "regular_ruling_set2", "regular_ruling_set2_rec", "euler_tour", "wood_regular_ruling_set"]

#for i in range(len(paths)):
    #generate_time_step_graph(paths[i], "time_step_" + names[i])

for i in range(len(paths)):    
    processors, times = get_PE_total_time_tuple(paths[i])
    
    processors = processors[7:]
    times = times[7:]
    
    width = 25
    processors_shifed = [p + i * width for p in processors]
    plt.bar(processors_shifed, times, label=names[i], width = width)

plt.xlabel("processors")
plt.ylabel("time in ms")
plt.title("time of different algorithms on random list")
plt.legend()



plt.savefig("all.pdf")



#generate_time_step_graph(path2)


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