import matplotlib.pyplot as plt
import csv
import numpy as np
import sys

def f(x):
    return (1-(np.sinh(2/x))**(-4))**(1/8)

def main():
    dt = open("Data/Data.csv", 'r')
    data = list(csv.reader(dt, delimiter = ','))
    dt.close()

    XEn, XMg, XCv, XChi = [], [], [], []
    YEn, YEnEr, YMg, YMgEr, YCv, YCvEr, YChi, YChiEr = [], [], [], [], [], [], [], []

    x = np.arange(0, 5+0.01, 0.01)
    
    for i in range(len(data)):
        XEn.append(float(data[i][0]))
        YEn.append(float(data[i][1]))
        YEnEr.append(float(data[i][2]))
        XMg.append(float(data[i][0]))
        YMg.append(float(data[i][3]))
        YMgEr.append(float(data[i][4]))
        XCv.append(float(data[i][0]))
        YCv.append(float(data[i][5]))
        YCvEr.append(float(data[i][6]))
        XChi.append(float(data[i][0]))
        YChi.append(float(data[i][7]))
        YChiEr.append(float(data[i][8]))
    
    if(sys.argv[1] == 0):        
        
        fig, axs = plt.subplots(2,2, figsize = (16,7))
        
        axs[0,0].scatter(XEn, YEn, color = 'black')
        axs[0,0].errorbar(XEn, YEn, YEnEr, color = 'black', fmt = 'o')
    
        axs[0,0].set_title(r"En. p/ part. modelo de Ising rede quadrada.")
        axs[0,0].set_xlim(0, 5)
        axs[0,0].set_ylim(-2.1, 0)
        axs[0,0].grid()
        ##########################################################3
    
        axs[0,1].plot(x, f(x), color = 'black', linestyle = 'dashed', label = 'Onsagers Solution')
        axs[0,1].scatter(XMg, YMg, color = 'black', label = 'Numerical Simulation')
        axs[0,1].errorbar(XMg, YMg, YMgEr, color = 'black', fmt = 'o')
    
        axs[0,1].set_title(r"Magnet. p/ part. modelo de Ising rede quadrada.")
        axs[0,1].set_xlim(0, 5)
        axs[0,1].set_ylim(0, 1.1)
        axs[0,1].legend()
        axs[0,1].grid()
    
        ##########################################################
    
        axs[1,0].scatter(XCv, YCv, color = 'black')
        axs[1,0].errorbar(XCv, YCv, YCvEr, color = "black")
        
        axs[1,0].set_title(r"Calor Esp. p/ part. modelo de Ising rede quadrada.")
        axs[1,0].set_xlim(0, 5)
        axs[1,0].grid()
        ##########################################################
    
        axs[1,1].scatter(XChi, YChi, color = 'black')
        axs[1,1].errorbar(XChi, YChi, YChiEr, color = "black")
    
        axs[1,1].set_title(r"Susc. Mag. p/ part. modelo de Ising rede quadrada.")
        axs[1,1].set_xlim(0, 5)
        axs[1,1].grid()
        ##########################################################
	
	
    if(sys.argv[1] == 1):        
        
        fig, axs = plt.subplots(figsize = (16,7))
        
        axs.scatter(XEn, YEn, color = 'black')
        axs.errorbar(XEn, YEn, YEnEr, color = 'black', fmt = 'o')
    
        axs.set_title(r"En. p/ part. modelo de Ising rede quadrada.")
        axs.set_xlim(0, 5)
        axs.set_ylim(-2.1, 0)
        axs.grid()
	

    if(sys.argv[1] == 2):        
        
        fig, axs = plt.subplots(figsize = (7,16))
        
        axs.plot(x, f(x), color = 'black', linestyle = 'dashed', label = 'Onsagers Solution')
        axs.scatter(XMg, YMg, color = 'black', label = 'Numerical Simulation')
        axs.errorbar(XMg, YMg, YMgEr, color = 'black', fmt = 'o')
    
        axs.set_title(r"Magnet. p/ part. modelo de Ising rede quadrada.")
        axs.set_xlim(0, 5)
        axs.set_ylim(0, 1.1)
        axs.legend()
        axs.grid()
	
	    
    if(sys.argv[1] == 3):        
        
        fig, axs = plt.subplots(figsize = (16,7))
    
        axs.scatter(XCv, YCv, color = 'black')
        axs.errorbar(XCv, YCv, YCvEr, color = "black")
        
        axs.set_title(r"Calor Esp. p/ part. modelo de Ising rede quadrada.")
        axs.set_xlim(0, 5)
        axs.grid()
	
	       
    if(sys.argv[1] == 4):        
        
        fig, axs = plt.subplots(figsize = (16,7))
    
        axs.scatter(XChi, YChi, color = 'black')
        axs.errorbar(XChi, YChi, YChiEr, color = "black")
    
        axs.set_title(r"Susc. Mag. p/ part. modelo de Ising rede quadrada.")
        axs.set_xlim(0, 5)
        axs.grid()
	
        ##########################################################
    plt.show()
    

main()
