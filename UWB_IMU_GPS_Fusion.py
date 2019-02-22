import matplotlib.pyplot as plt
import numpy as np
from UAV import UAV
from UAV import StateEstimator

def main():
    slope, R, v = 1/np.sqrt(3), 100.0, 10.0
    t = np.arange(0, 2*np.pi*R*3/v, 0.01)
    N = t.size

    uav1 = UAV(lambda t: R*np.cos((v/R)*t), lambda t: R*np.sin((v/R)*t), lambda t: v*t*slope)
    uav2 = UAV(lambda t: -R*np.sin((v/R)*t), lambda t: R*np.cos((v/R)*t), lambda t: v*t*slope)
    uav3 = UAV(lambda t: -R*np.cos((v/R)*t), lambda t: -R*np.sin((v/R)*t), lambda t: v*t*slope)
    uav4 = UAV(lambda t: R*np.sin((v/R)*t), lambda t: -R*np.cos((v/R)*t), lambda t: v*t*slope)
    # Initial alignment
    uav1State = StateEstimator(uav1.getR(0), uav1.getV(0), uav1.getA(0), uav1.getAltitudeQab(0))
    uav2State = StateEstimator(uav2.getR(0), uav2.getV(0), uav2.getA(0), uav2.getAltitudeQab(0))
    uav3State = StateEstimator(uav3.getR(0), uav3.getV(0), uav3.getA(0), uav3.getAltitudeQab(0))
    uav4State = StateEstimator(uav4.getR(0), uav4.getV(0), uav4.getA(0), uav4.getAltitudeQab(0))

    for time in t:
        uav1State.UpdateState(uav1.IMUmeasure(time), time, uav1.Magmeasure(time), uav1.GPSmeasure(time))
        uav2State.UpdateState(uav2.IMUmeasure(time), time, uav2.Magmeasure(time), uav2.GPSmeasure(time))
        uav3State.UpdateState(uav3.IMUmeasure(time), time, uav3.Magmeasure(time), uav3.GPSmeasure(time))
        uav4State.UpdateState(uav4.IMUmeasure(time), time, uav4.Magmeasure(time), uav4.GPSmeasure(time))

    plt.plot(t, uav1.GPSmeasure)
    plt.show()

main()
