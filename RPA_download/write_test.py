import numpy as np
import csv
data_folder='March_16/'
def Pi_om_in_file(qqx):
        omega1 = np.array([1,2,3,4])
        Pi_r = np.array([1,2,3,4])
        Pi_i=np.array([5,6,7,8])
        file = open(data_folder+"qx="+str(round(qqx,2))+",real.csv", "w")
        writer = csv.writer(file)

        for w in range(np.size(omega1)): #from 0-3
                writer.writerow([omega1[w], Pi_r[w]])
        file.close()

        file = open(data_folder+"qx="+str(round(qqx,2))+",imag.csv", "w")
        writer = csv.writer(file)
        for w in range(np.size(omega1)): #from 0-3
                writer.writerow([omega1[w], Pi_i[w]])

        file.close()
Pi_om_in_file(23333)
