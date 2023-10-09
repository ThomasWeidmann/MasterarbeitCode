with open('build/results.txt', 'r+') as file: 
 file_data = file.read() 
 file.seek(0, 0) 
 file.write("{\n\"data\":[\n" + file_data[:-2] + "\n]}") 