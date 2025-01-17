import control_rod_automation
import io
from contextlib import redirect_stdout
import csv
import numpy as np
from time import sleep
 
if __name__ == "__main__": 
    """
    # Shim 1
    s1 = np.array([0, 11.1, 21.1, 29.9, 40.7, 49.3, 60.1, 70.3, 70.3, 74.5, 74.5, 78.7, 78.7, 83.8, 83.8, 92.1, 92.1, 100, 100])
    s2 = np.array([100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100])
    rr = np.array([89.5, 89.5, 89.5, 89.5, 89.5, 89.5, 89.5, 89.5, 100, 89.5, 100, 82.7, 89.5, 76.7, 82.7, 71.1, 76.7, 66.9, 71.1])
    
    # Save .csv of shim1 data
    filename = "shim1_rodworths.csv"
    with open(filename, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Headers
        writer.writerow(['Shimrod1', 'Reg Rod', 'Shimrod2', 'k-eff'])

        for i in range(len(s1)):
            # Capture OpenMC Output
            f = io.StringIO()
            with redirect_stdout(f):
                control_rod_automation.cr_sliding(s1[i], s2[i], rr[i])
            
            output = f.getvalue()
            
            # Record Combined k-eff
            output_list = output.splitlines()
            k_eff = output_list[len(output_list) - 3][31:38] # Select Combined k-eff

            # Record data
            writer.writerow([s1[i], rr[i], s2[i], k_eff])

    print("Shimrod 1 Complete!")
    """
    # Reg Rod
    s1 = np.array([90.6, 90.6, 90.6, 90.6, 90.6, 90.6, 90.6, 90.6, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100])
    s2 = np.array([94.1, 94.1, 94.1, 94.1, 94.1, 94.1, 94.1, 94.1, 100, 100, 92.3, 92.3, 85.8, 85.8, 79.1, 79.1, 74.1, 74.1, 70.3])
    rr = np.array([0, 10, 20, 30, 40, 50, 60, 66.3, 66.3, 69.1, 69.1, 73.2, 73.2, 80.2, 80.2, 88.4, 88.4, 100, 100])
    
    # Save .csv of Reg Rod data
    filename = "regrod_rodworths.csv"
    with open(filename, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Headers
        writer.writerow(['Shimrod1', 'Reg Rod', 'Shimrod2', 'k-eff'])

        for i in range(len(s1)):
            # Capture OpenMC Output
            f = io.StringIO()
            with redirect_stdout(f):
                control_rod_automation.cr_sliding(s1[i], s2[i], rr[i])
            
            output = f.getvalue()
            
            # Record Combined k-eff
            output_list = output.splitlines()
            k_eff = output_list[len(output_list) - 3][31:38] # Select Combined k-eff

            # Record data
            writer.writerow([s1[i], rr[i], s2[i], k_eff])

    print("Reg Rod Complete!")
    print("ALL DONE!")