import subprocess
import os
for i in [1, 1.25, 1.5, 1.75, 2]:
    path="result\\Maximum_Trip_Coefficient="+str(i)
    isExist = os.path.exists(path)
    if not isExist:
        # Create a new directory because it does not exist 
        os.makedirs(path)
        print("The new directory is created: "+path)
    subprocess.run(["./idarp.exe", "input/one_hundred_cust.txt", "input/ListeOfLocations_10000_Cust.txt" ,"input/Discretisation_10000_Cust.txt","input/load_file.txt", str(i)])


