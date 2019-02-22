import numpy as np
# The Magnatic field intensity (North, East, Down)/nT
Mag = np.array((27735, -3361, 47025))
gravity = 9.8

def DCMorthonormalization(DCM):
    """ This function normalize the column vectors in the DCM matrix and make them orthogonal 
        to each other."""
    Z = np.cross(DCM[:,0], DCM[:,1])
    Y = np.cross(Z, DCM[:,0])
    DCM = np.array((DCM[:,0],Y,Z)).T 
    ColumnNorm = np.linalg.norm(DCM, None, 0)
    DCM = DCM/ColumnNorm

    return DCM

def DCMToQuaternion(DCM):
    """ This function convert the altitude DCM representation to quaternion representation, this function assume
        that the DCM is orthonormal, that is, the column vectors in the matrix are all unit vectors and orthogonal 
        to each other."""
    Q=np.array((DCM.trace(), DCM[0,0]-DCM[1,1]-DCM[2,2], -DCM[0,0]+DCM[1,1]-DCM[2,2], -DCM[0,0]-DCM[1,1]+DCM[2,2]))
    Q=np.sqrt(Q+1)/2
    if (DCM[1,2]-DCM[2,1])<0:
        Q[1]=-Q[1]
    if (DCM[2,0]-DCM[0,2])<0:
        Q[2]=-Q[2]
    if (DCM[0,1]-DCM[1,0])<0:
        Q[3]=-Q[3]

    return Q

def QuaternionToDCM(Q):
    Q = Q/np.linalg.norm(Q)
    Qmat = np.array([[0,Q[3],-Q[2]], [-Q[3],0,Q[1]], [Q[2],-Q[1],0]])

    return np.eye(3)+2*Q[0]*Qmat+2*np.dot(Qmat,Qmat)

def QuaternionMul(Q1, Q2):
    #c1=a1*b1-a2*b2-a3*b3-a4*b4
    #c2=-a4*b3+a3*b4+a1*b2+a2*b1
    #c3=-a2*b4+a1*b3+a3*b1+a4*b2
    #c4=-a3*b2+a1*b4+a2*b3+a4*b1
    return np.array((Q1[0]*Q2[0]-Q1[1]*Q2[1]-Q1[2]*Q2[2]-Q1[3]*Q2[3],\
                     -Q1[3]*Q2[2]+Q1[2]*Q2[3]+Q1[0]*Q2[1]+Q1[1]*Q2[0],\
                     -Q1[1]*Q2[3]+Q1[0]*Q2[2]+Q1[2]*Q2[0]+Q1[3]*Q2[1],\
                     -Q1[2]*Q2[1]+Q1[0]*Q2[3]+Q1[1]*Q2[2]+Q1[3]*Q2[0],\
       ))

class UAV(object):
    """An UAV object with predefined trajectory, velocity and acceleration."""
    def __init__(self, rx, ry, rz, vx=None, vy=None, vz=None, ax=None, ay=None, az=None):
        self.hrx=rx
        self.hry=ry
        self.hrz=rz
        self.hvx=vx
        self.hvy=vy
        self.hvz=vz
        self.hax=ax
        self.hay=ay
        self.haz=az

    def getR(self, t):
        return np.row_stack((self.hrx(t), self.hry(t), self.hrz(t)))

    def getV(self, t):
        if self.hvx==None:
            vx=(self.hrx(t+0.001)-self.hrx(t-0.001))/0.002
        else:
            vx=self.hvx(t)

        if self.hvy==None:
            vy=(self.hry(t+0.001)-self.hry(t-0.001))/0.002
        else:
            vy=self.hvy(t)

        if self.hvz==None:
            vz=(self.hrz(t+0.001)-self.hrz(t-0.001))/0.002
        else:
            vz=self.hvz(t)

        return np.row_stack((vx, vy, vz))

    def getA(self, t):
        if self.hax==None:
            if self.hvx==None:
                ax=(self.hrx(t+0.001)+self.hrx(t-0.001)-2*self.hrx(t))/1e-6
            else:
                ax=(self.hvx(t+0.001)-self.hvx(t-0.001))/0.002
        else:
            ax=self.hax(t)

        if self.hay==None:
            if self.hvy==None:
                ay=(self.hry(t+0.001)+self.hry(t-0.001)-2*self.hry(t))/1e-6
            else:
                ay=(self.hvy(t+0.001)-self.hvy(t-0.001))/0.002
        else:
            ay=self.hay(t)

        if self.haz==None:
            if self.hvz==None:
                az=(self.hrz(t+0.001)+self.hrz(t-0.001)-2*self.hrz(t))/1e-6
            else:
                az=(self.hvz(t+0.001)-self.hvz(t-0.001))/0.002
        else:
            az=self.haz(t)
            
        return np.row_stack((ax, ay, az))

    def getAltitudeCab(self, t):
        if np.size(t) > 1:
            t = t[0]
            print('Only a scalar t is allowed, the following calculation will set t = t[0] =', t, '.')

        eF = self.getV(t).T
        eF = eF/np.linalg.norm(eF, None, 1)

        aInertial_a = -self.getA(t).T
        aInertial_a[0, 2] = aInertial_a[0, 2]-gravity
        eL = np.cross(eF, aInertial_a)
        eL = eL/np.linalg.norm(eL, None, 1)

        eT = np.cross(eF, eL)
        return np.row_stack((eF, eL, eT)).T

    def getAltitudeQab(self, t):
        if np.size(t) > 1:
            t = t[0]
            print('Only a scalar t is allowed, the following calculation will set t = t[0] =', t, '.')

        return DCMToQuaternion(self.getAltitudeCab(t))

    def getAngularRate_ba_b(self, t):
        Theta=np.dot(self.getAltitudeCab(t).T,(self.getAltitudeCab(t+0.0005)-self.getAltitudeCab(t-0.0005))/0.001)

        return np.array((Theta[2,1], Theta[0,2], Theta[1,0]))

    def IMUmeasure(self,T):
        #Add zero-mean Gausian noise to the acceleration and angular rate measurement data according to the datasheet.
        N = np.size(T)
        IMUdata=np.empty((6, N))
        if N < 2:
            aInertial_a=-self.getA(T)
            aInertial_a[2]=aInertial_a[2]-gravity
            IMUdata[0:3,0]=np.ravel(np.dot(self.getAltitudeCab(T).T, aInertial_a))+np.random.normal(0,0.05,3)
            IMUdata[3:6,0]=np.ravel(self.getAngularRate_ba_b(T))+np.random.normal(0,0.014,3)
        else:
            for i, t in enumerate(T):
                #a_b = Cba * a_a 
                aInertial_a=-self.getA(t)
                aInertial_a[2]=aInertial_a[2]-gravity
                IMUdata[0:3,i]=np.ravel(np.dot(self.getAltitudeCab(t).T, aInertial_a))+np.random.normal(0,0.05,3)
                IMUdata[3:6,i]=np.ravel(self.getAngularRate_ba_b(t))+np.random.normal(0,0.014,3)

        return IMUdata

    def Magmeasure(self,T):
        N = np.size(T)
        Magdata=np.empty((3, N))
        if N < 2:
            Magdata[:,0]=np.ravel(np.dot(self.getAltitudeCab(T).T, Mag))+np.random.normal(0, 350, 3)
        else:
            for i, t in enumerate(T):
                Magdata[:,i]=np.ravel(np.dot(self.getAltitudeCab(t).T, Mag))+np.random.normal(0, 350, 3)
        
        return Magdata

    def GPSmeasure(self,t):
        return self.getR(t)+np.row_stack((np.random.normal(0,2,(2,np.size(t))), np.random.normal(0,5,(1,np.size(t)))))

    def UWBmeasure(self, otherUAV, t):
        return np.abs(np.linalg.norm(self.getR(t)-otherUAV.getR(t),None,0)+np.random.normal(0,0.3,np.size(t)))

class StateEstimator(object):
    """An estimator for position, velocity, acceleration and altitude based on IMU, GPS, Compass and UWB measurement"""
    def __init__(self, R = np.array((0.0, 0.0, 0.0)), V = np.array((0.0, 0.0, 0.0)), A = np.array((0.0, 0.0, 0.0)),		\
				 Q = np.array((1.0, 0.0, 0.0, 0.0)), T = 0.0):
        self.r = R
        self.v = V
        self.a = A
        self.Qab = Q
        self.t = T

    def SetAltitude(self, altitude):
        # If we get a quaternion to represent the altitude
        if altitude.shape == (4,) or altitude.shape == (4,1):
            self.Qab = altitude/np.linalg.norm(altitude)
            if altitude[0]<0:
                self.Qab = -self.Qab
        # If we get a Direction Cosine Matrix to represent the altitude
        elif altitude.shape == (3,3):
            self.Qab = DCMToQuaternion(DCMorthonormalization(altitude))

    def GetAltitudeQ(self):
        return self.Qab

    def GetAltitudeDCM(self):
        return QuaternionToDCM(self.Qab)

    def GetP(self):
        return self.r

    def GetV(self):
        return self.v

    def GetA(self):
        return self.a

    def SetP(self, position):
        R = np.array(position)
        if R.shape==(3,):
            self.r=R
        elif R.shape==(3,1) or R.shape==(1,3):
            self.r=np.ravel(R)
        else:
            print("Only array with 3 elements is acceptable.")

    def resetState(self):
        self.r = np.array((0.0, 0.0, 0.0))
        self.v = np.array((0.0, 0.0, 0.0))
        self.a = np.array((0.0, 0.0, 0.0))
        self.Qab = np.array((1.0, 0.0, 0.0, 0.0))
        self.t = 0.0

    def UpdateState(self, IMU, Time, Mag=None, GPS=None, UWB=None):
        phy=IMU[3:6, 0]
        physquare_8=np.dot(phy,phy)/8.0
        qab=np.hstack([1-physquare_8, (physquare_8/6.0-0.5)*phy])