import matplotlib.pyplot as plt
import numpy as np
from UAV import UAV

def main():
    slope, R, v = 1/np.sqrt(3), 100, 10
    t=np.arange(0, 2*np.pi*R*3/v,0.01)
    N=t.size

    uav1=UAV(lambda t:R*np.cos((v/R)*t), lambda t:R*np.sin((v/R)*t), lambda t:v*t*slope)
    uav2=UAV(lambda t:-R*np.sin((v/R)*t), lambda t:R*np.cos((v/R)*t), lambda t:v*t*slope)
    print(uav1.getAngularRate_ba_b(np.array((1,2,3,4,5))).T)

main()