from packages.Imports import *
from packages.qpe_protocols.qcels import *

#This script is used to simulate QCELS applied to the ground state energy estimation of H2

def Ham(g): #Hamiltion derivied in https://arxiv.org/pdf/1512.06860.pdf
    Z = np.array([[1,0],[0,-1]])
    X = np.array([[0,1],[1,0]])
    Y = np.array([[0,-1j],[1j,0]])
    I = np.array([[1,0],[0,1]])
    H = (g[1]*np.kron(I, Z)) + (g[2]*np.kron(Z, I)) +  (g[4]*np.kron(X,X)) + (g[5]*np.kron(Y, Y))
    return H

R = np.linspace(0.2,2.85,54) #interatomic distances

#parameters for Hamiltonian coefficients based off of interatomic distance
g = np.array([[2.8489, 0.5678, -1.4508, 0.6799, 0.0791, 0.0791],
              [2.1868, 0.5449, -1.2870, 0.6719, 0.0798, 0.0798],
              [1.7252, 0.5215, -1.1458, 0.6631, 0.0806, 0.0806],
              [1.3827, 0.4982, -1.0226, 0.6537, 0.0815, 0.0815],
              [1.1182, 0.4754, -0.9145, 0.6438, 0.0825, 0.0825],
              [0.9083, 0.4534, -0.8194, 0.6336, 0.0835, 0.0835],
              [0.7381, 0.4325, -0.7355, 0.6233, 0.0846, 0.0846],
              [0.5979, 0.4125, -0.6612, 0.6129, 0.0858, 0.0858],
              [0.4808, 0.3937, -0.5950, 0.6025, 0.0870, 0.0870],
              [0.3819, 0.3760, -0.5358, 0.5921, 0.0883, 0.0883],
              [0.2976, 0.3593, -0.4826, 0.5818, 0.0896, 0.0896],
              [0.2252, 0.3435, -0.4347, 0.5716, 0.0910, 0.0910],
              [0.1626, 0.3288, -0.3915, 0.5616, 0.0925, 0.0925],
              [0.1083, 0.3149, -0.3523, 0.5518, 0.0939, 0.0939],
              [0.0609, 0.3018, -0.3168, 0.5421, 0.0954, 0.0954],
              [0.0193, 0.2895, -0.2845, 0.5327, 0.0970, 0.0970],
              [-0.0172, 0.2779, -0.2550, 0.5235, 0.0986, 0.0986],
              [-0.0493, 0.2669, -0.2282, 0.5146, 0.1002, 0.1002],
              [-0.0778, 0.2565, -0.2036, 0.5059, 0.1018, 0.1018],
              [-0.1029, 0.2467, -0.1810, 0.4974, 0.1034, 0.1034],
              [-0.1253, 0.2374, -0.1603, 0.4812, 0.1050, 0.1050],
              [-0.1452, 0.2286, -0.1413, 0.4812, 0.1067, 0.1067],
              [-0.1629, 0.2203, -0.1238, 0.4735, 0.1083, 0.1083],
              [-0.1786, 0.2123, -0.1077, 0.4660, 0.1100, 0.1100],
              [-0.1927, 0.2048, -0.0929, 0.4588, 0.1116, 0.1116],
              [-0.2053, 0.1976, -0.0792, 0.4518, 0.1133, 0.1133],
              [-0.2165, 0.1908, -0.0666, 0.4451, 0.1149, 0.1149],
              [-0.2265, 0.1843, -0.0549, 0.4386, 0.1165, 0.1165],
              [-0.2355, 0.1782, -0.0442, 0.4323, 0.1181, 0.1181],
              [-0.2436, 0.1723, -0.0342, 0.4262, 0.1196, 0.1196],
              [-0.2508, 0.1667, -0.0251, 0.4204, 0.1211, 0.1211],
              [-0.2573, 0.1615, -0.0166, 0.4148, 0.1226, 0.1226],
              [-0.2632, 0.1565, -0.0088, 0.4094, 0.1241, 0.1241],
              [-0.2684, 0.1517, -0.0015, 0.4042, 0.1256, 0.1256],
              [-0.2731, 0.1472, 0.0052, 0.3992, 0.1270, 0.1270],
              [-0.2774, 0.1430, 0.0114, 0.3944, 0.1284, 0.1284],
              [-0.2812, 0.1390, 0.0171, 0.3898, 0.1297, 0.1297],
              [-0.2847, 0.1352, 0.0223, 0.3853, 0.1310, 0.1310],
              [-0.2879, 0.1316, 0.0272, 0.3811, 0.1323, 0.1323],
              [-0.2908, 0.1282, 0.0317, 0.3769, 0.1335, 0.1335],
              [-0.2934, 0.1251, 0.0359, 0.3730, 0.1347, 0.1347],
              [-0.2958, 0.1221, 0.0397, 0.3692, 0.1359, 0.1359],
              [-0.2980, 0.1193, 0.0432, 0.3655, 0.1370, 0.1370],
              [-0.3000, 0.1167, 0.0465, 0.3620, 0.1381, 0.1381],
              [-0.3018, 0.1142, 0.0495, 0.3586, 0.1392, 0.1392],
              [-0.3035, 0.1119, 0.0523, 0.3553, 0.1402, 0.1402],
              [-0.3051, 0.1098, 0.0549, 0.3521, 0.1412, 0.1412],
              [-0.3066, 0.1078, 0.0572, 0.3491, 0.1422, 0.1422],
              [-0.3079, 0.1059, 0.0594, 0.3461, 0.1432, 0.1432],
              [-0.3092, 0.1042, 0.0614, 0.3433, 0.1441, 0.1441],
              [-0.3104, 0.1026, 0.0632, 0.3406, 0.1450, 0.1450],
              [-0.3115, 0.1011, 0.0649, 0.3379, 0.1458, 0.1458],
              [-0.3125, 0.0997, 0.0665, 0.3354, 0.1467, 0.1467],
              [-0.3135, 0.0984, 0.0679, 0.3329, 0.1475, 0.1475]])

#I'm currently handling state prep by constructing a projector from |00> to H2 ground state.
def eigensystem(h):
    l,e = np.linalg.eig(h)
    ll = np.zeros(4)
    ee = np.zeros([4,4],dtype = 'complex_')
    for i in range(4):
        m = l.argmin()
        ll[i] = l[m]
        ee[:,i] = e[:,m]
        l = np.delete(l,m)
        e = np.delete(e,m,1)
    return ll,ee

def P(h): #Creates matrix projector from basis states to hamilitonian eigenstates 
    I = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    PP = np.zeros([4,4])
    l,e = eigensystem(h)
    for i in range(4):
        PP = PP + np.outer(e[:,i],I[i,:])
    return PP

def Usp(circuit): #Turns projector into UnitaryGate object and applies to input quantum circuit
    h = Ham(g[15]) #I choose the bond distance arbitrarily although it's somewhat close to equilibirium 
    p = P(h)
    proj = UnitaryGate(p)
    #circuit.rx(pi/3,2) #rotation controls target state overlap
    circuit.append(proj,[1,2]) #Apply to qubits 1&2 as they are the system qubits, qubit 0 is ancilla qubit

def W(circuit): #exact time evolution operator (used for debugging)
    tjj = pow(2,-J-1)/(N*target)
    unitary = scipy.linalg.expm(-1j*Ham(g[15])*tjj) #hbar=1
    U = UnitaryGate(unitary)
    cU = U.control(1)
    circuit.append(cU,[0,1,2])

def W_trot(circuit):  #trotterized time evolution operator, time must be chosen so trotter error is not a limiting factor
    t = pow(2,-1)
    gTrot = g[15]
    a1 = -2*gTrot[1]*t
    a2 = -2*gTrot[2]*t
    a3 = -2*gTrot[3]*t
    a4 = -2*gTrot[4]*t
    a5 = -2*gTrot[5]*t

    #XX
    circuit.h(1)
    circuit.h(2)
    circuit.cx(2,1)
    circuit.cx(0,1)
    circuit.rz(a4/2,1)
    circuit.cx(0,1)
    circuit.cx(2,1)
    circuit.h(1)
    circuit.h(2)
    #Z1
    circuit.cx(0,1)
    circuit.rz(a1/2,1)
    circuit.cx(0,1)
    #YY
    circuit.sdg(1)
    circuit.sdg(2)
    circuit.h(1)
    circuit.h(2)
    circuit.cx(2,1)
    circuit.cx(0,1)
    circuit.rz(a5/2,1)
    circuit.cx(0,1)
    circuit.cx(2,1)
    circuit.h(1)
    circuit.h(2)
    circuit.s(1)
    circuit.s(2)
    #Z2
    circuit.cx(0,2)
    circuit.rz(a2/2,2)
    circuit.cx(0,2)

target = np.logspace(-1,-10,num=10,base=2); #Target errors
samples = 10;
N = 2 #I found RPE always has better T counts at same Tmax unless N is set to 2
Ns = 55
lc,ec = eigensystem(Ham(g[15]))
errors = np.zeros([len(target),samples])
error = np.zeros(len(target))
Tmax = np.zeros(len(target))
Ttot = np.zeros(len(target))

#depolarizing noise model
p_hardware = 0.001
b = 35*pow(p_hardware,3) #assume 15:1 T state protocol
T_rate = 11
d = 17
n_L = 19
prob_err = b + T_rate*n_L*d*0.1*pow(100*p_hardware,(d+1)/2) #probability of error per T gate
DPerror = depolarizing_error(prob_err,1);
nm = NoiseModel()
nm.add_all_qubit_quantum_error(DPerror, ['t']) #We are assuming Clifford gates have been commuted away

for l in range(len(target)):
    J = math.ceil(np.log2(1/target[l]))+1
    rNs = sqrt(Ns)
    if(rNs>4): #If state overlap=1 utilize error reduction from Ns to reduce J NOTE: come back here to implement better
    	J = max(J-2,1)
    elif(rNs>2 and rNs<4):
    	J = max(J-1,1)
    wCallsTot = pow(2,J)-1
    synPrecision = 0.9*target[l]/(Ns*wCallsTot*4) #determine total number of Z-rotations in order to determine necessary rotation sysnthesis precision
    tjj = pow(2,-1)
    for k in range(samples):
        rep = 1
        for i in range(J):
            tj = pow(2,i)*tjj
            t = np.linspace(0,tj*(N-1),N)
            sim = noisy_qcels(N,Ns,W_trot,Usp,3,rep,target[l],synPrecision);
            z = sim[0]
            if (i==0):
                est = (-1*atan2(np.imag(z[1]),np.real(z[1]))) #I found best results seeding first least squares fit with atan2
                bnds = [(-1,1),(-1,1),(-pi,pi)] #Initial bounds [-2pi, pi)
            zfit = maxF(z,N,t,est,bnds)
            rep = rep*2
            est = zfit[2] 
            bnds = [(-1,1),(-1,1),(est-pi/(2*tj),est+pi/(2*tj))] #set new bounds
        errors[l][k] = (abs(est-lc[0]))
    error[l] = np.mean(errors[l])
    Ttot[l] += Ns*sim[1]
    Tmax[l] = sim[1]/2


#np.savetxt("PaperqcelErrorp1.csv",error,delimiter=",")
#np.savetxt("PaperqcelTtotp1.csv",Ttot,delimiter=",")
#np.savetxt("PaperqcelTmaxp1.csv",Tmax,delimiter=",")



#
