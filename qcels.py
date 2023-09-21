from Imports import *
from scipy.optimize import minimize, shgo, least_squares, differential_evolution
from cmath import exp

def qcels(N, Ns, W, Usp, numQubits, repetitions,target):

    z = np.zeros(N,dtype=complex)
    wCalls = 0

    #stateprep = QuantumCircuit(numQubits)
    #Usp(stateprep)
    #stateprep = compileCT(stateprep,1e-4)
    #Tsp = stateprep.count_ops().get('t')
    for l in range(N):
        print(l)
        qcI = QuantumCircuit(numQubits)
        qcS = QuantumCircuit(numQubits)
        measurement = QuantumCircuit(numQubits,1)
        qcI.h(0)
        qcS.h(0)

        Usp(qcI)
        Usp(qcS)

        for i in range(l):
            for k in range(repetitions):
                W(qcI,target)
                W(qcS,target)
                wCalls += 2
        qcS.sdg(0)
        qcI.h(0)
        qcS.h(0)
        #print(qcI.draw())
        #qcI = compileCT(qcI,0.9*target)
        #qcS = compileCT(qcS,0.9*target)
        #numD1 = qcI.count_ops().get('t')
        #numD2 = qcS.count_ops().get('t')

        measurement.barrier(range(numQubits))
        measurement.measure(0,0)

        qcI = qcI.compose(measurement)
        qcS = qcS.compose(measurement)
        #print(qcI.draw())

        aer_sim = Aer.get_backend('aer_simulator')
        qcI_compiled = transpile(qcI, aer_sim)
        runs = Ns;
        job_sim = aer_sim.run(qcI_compiled, shots=runs)
        result_sim = job_sim.result()
        counts = result_sim.get_counts(qcI_compiled)

        if '0' in counts.keys():
            x = 2*counts['0']/Ns -1
        else:
            x = -1

        aer_sim = Aer.get_backend('aer_simulator')
        qcS_compiled = transpile(qcS, aer_sim)
        runs = Ns;
        job_sim = aer_sim.run(qcS_compiled, shots=runs)
        result_sim = job_sim.result()
        counts = result_sim.get_counts(qcS_compiled)

        if '0' in counts.keys():
            y = 2*counts['0']/Ns -1
        else:
            y = -1

        z[l] = complex(x,y)
    return [z, wCalls]

def F(x, z, N, t):
    Z_fit=np.zeros(N,dtype = 'complex_')
    Z_fit=(x[0]+1j*x[1])*np.exp(-1j*x[2]*t)
    return (np.linalg.norm(Z_fit-z)**2/N)

def maxF(z, N, t, est, bnds):
    fun = lambda x: F(x, z, N, t)
    result = minimize(fun,x0=[0.5,0,est],method = 'SLSQP',bounds=bnds)
    #result = differential_evolution(fun, bounds = [(-1,1),(-1,1),(0,2*pi)])
    return result.x


def qcelsT(N, Ns, W, Usp, numQubits, repetitions):

    z = np.zeros(N,dtype=complex)
    Tmax = 0
    T = 0

    for l in range(N):

        #print(l)
        qcI = QuantumCircuit(numQubits)
        qcS = QuantumCircuit(numQubits)
        measurement = QuantumCircuit(numQubits,1)

        qcI.h(0)
        qcS.h(0)

        Usp(qcI)
        Usp(qcS)

        for i in range(l):
            for k in range(repetitions):
                W(qcI)
                W(qcS)

        qcS.sdg(0)
        qcI.h(0)
        qcS.h(0)

        qcI = compileCT(qcI,1e-4)
        qcS = compileCT(qcS,1e-4)

        measurement.barrier(range(numQubits))
        measurement.measure(0,0)

        qcI = qcI.compose(measurement)
        qcS = qcS.compose(measurement)

        #print(qcI.draw())

        numD1 = qcI.count_ops().get('t')
        numD2 = qcS.count_ops().get('t')

        T = numD1 + numD2;
        Tmax = max(Tmax,numD1,numD2);

    return [T*Ns, Tmax]


"""
N = 10
rep = 1
t = [10*2^(-5),10*2^-4,10*2^-3,10*2^-2,10*2^-1,10];
Ns = 600

for i in range(6):

    z = qcels(N,Ns,W,Usp,3,rep);
    if (i==0):
        est = (-1*atan2(np.imag(z[1]),np.real(z[1])))%(2*pi)
        bnds = [(-1,1),(-1,1),(0,2*pi)]

    zfit = maxF(z,N,t,est,bnds)
    rep = rep*2
    est = zfit[2]
    bnds = [(-1,1),(-1,1),(est-pi/rep,est+pi/rep)]
    print(est)
    print([theta ])




print(theta)
"""



#print(t)









#
