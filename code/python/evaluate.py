#very important: this skript handles json objects, the syntax is obviously json, but more specific 
#{"data":[data_point1, data_point2]}


import json
import numpy as np
import matplotlib.pyplot as plt
 


path1 = '../../other/supermuc_auswertung/regular_ruling_set.txt'
path2 = '../../other/supermuc_auswertung/regular_ruling_set_rec.txt'
path3 = '../../other/supermuc_auswertung/regular_ruling_set2.txt'
path4 = '../../other/supermuc_auswertung/regular_ruling_set2_rec.txt'
path5 = '../../other/supermuc_auswertung/euler_tour.txt'
path6 = '../../other/supermuc_auswertung/wood_regular_ruling_set2.txt'
path7 = '../../other/supermuc_auswertung/regular_pointer_doubling.txt'


paths = [path1, path2, path3, path4, path5, path6, path7]
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

to_json(path7)
to_json(path6)
#for path in paths:
    #to_json(path)

for i in range(len(paths)):
    
    generate_time_step_graph(paths[i], "time_step_" + names[i])

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



f.close()