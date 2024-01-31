from Imports import *
from scipy.optimize import minimize
from cmath import exp

#This script provides methods to simulate multi-modal quantum complex exponetenial least squares (MMQCELS) as described in https://arxiv.org/pdf/2303.05714.pdf 
#Two versions of this algorithm are provided; "MMQcels" performs noiseless simulations without compiling circuits into Clifford+Ts while "noisy_MMQcels" allows for the addition
#of a noise model as well as compiles all circuits into Clifford+Ts while tracking the Tmax and Ttot.
#Since MMQcels does not compile circuits into Clifford+Ts it runs much faster.
#Because Ns is typically much lower, MMQcels takes longer to simulate using qiskit than qcels or arpe. 

#Both MMQcels and noisy_MMQcels utilize a dataGenerator function. 
#This function generates N data pairs 
#MMQcels and noisy_MMQcels calls the dataGenerator function then performs a least squares fit on the generated data
#Additionally noisy_MMQcels also returns the Ttot and Tmax.

#N = number of data pairs, Ns = number of samples to produce each data pair, sigmaT = width of normal distribution from which Hamiltonian simulation time is selected,
#T = maximum Hamiltonian simulation time, W = function to implement time evolution operator, Usp = function to implement state prep, numQubits = number of qubits including ancilla, gen = current generation 
def dataGenerator(N, Ns, sigmaT, T, W, Usp, numQubits, gen):

    z = np.zeros(N,dtype=complex)
    t = np.zeros(N)
    k = pow(2,gen-1)

    for i in range(N):
        wCN = 0
        qcI = QuantumCircuit(numQubits)
        qcS = QuantumCircuit(numQubits)
        measurement = QuantumCircuit(numQubits,1)

        qcI.h(0)
        qcS.h(0)
        #stateprep
        Usp(qcI)
        Usp(qcS)
        #Random Time Evoloution
        while(True): #instead of sampling from a truncated Gaussian, I will throw out a t value greater than a set maximum
            t[i] = np.random.normal(0,sigmaT,1)
            if(abs(t[i]) < T):
                break

        num = abs(k*t[i])//T
        rem = abs(k*t[i])%T
        if(t[i] >= 0):
            for j in range(int(num)):
                W(qcI,T)
                W(qcS,T)
            W(qcI,rem)
            W(qcS,rem)
        elif(t[i] < 0):
            for j in range(int(num)):
                W(qcI,-1*T)
                W(qcS,-1*T)

            W(qcI,-1*rem)
            W(qcS,-1*rem)
            wCN += 2
        t[i] = t[i]*k


        qcS.sdg(0)
        qcI.h(0)
        qcS.h(0)

        measurement.barrier(range(numQubits))
        measurement.measure(0,0)

        qcI = qcI.compose(measurement)
        qcS = qcS.compose(measurement)

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

        z[i] = complex(x,y)

    return [t,z] #return times and z

def F(x, z, N, t): #Define function to minimized for least squares fit
    Z_fit=np.zeros(N,dtype = 'complex_')
    Z_fit=(x[0]+1j*x[1])*np.exp(-1j*x[2]*t)
    return (np.linalg.norm(Z_fit-z)**2/N)

def maxF(z, N, t, est, bnds): #z = time series data generated form qcels, N = number of data pairs, t = corresponding times for time series, est = prior estimate, bnds = bounds for fit
    fun = lambda x: F(x, z, N, t)
    result = minimize(fun,x0=[0.5,0,est],method = 'SLSQP',bounds=bnds)
    return result.x

#N = number of data pairs per generation, number of samples per data pair, sigmaT = initial width of normal distrubtion of Hamiltonian simulation time (doubles every subsequent generation)
#T = initial maximum Hamiltonian simulation time (doubles every subsequent generation, W = function to implement time evolution operator, Usp = function to implement state prep, numQubits = number of qubits including ancilla, gens = number of generations to run
def MMQcels(N, Ns, sigmaT, T, W, Usp, numQubits, gens):
    bnds = [(-1,1),(-1,1),(-pi,pi)]
    est = 0
    for i in range(gens):
        z = dataGenerator(N,Ns,sigmaT,T,W,Usp,numQubits,i+1)
        zfit = maxF(z[1],N,z[0],est,bnds)
        est = zfit[2]
        bnds = [(-1,1),(-1,1),(est-pi/(pow(2,i)*sigmaT),est+pi/(pow(2,i)*sigmaT))
    return est

#N = number of data pairs, Ns = number of samples to produce each data pair, sigmaT = width of normal distribution from which Hamiltonian simulation time is selected,
#T = maximum Hamiltonian simulation time, W = function to implement time evolution operator, Usp = function to implement state prep, numQubits = number of qubits including ancilla, gen = current generation 
#nm = noise model, precision = rotation synthesis precision for each Z-rotation
def noisyDataGenerator(N, Ns, sigmaT, T, W, Usp, numQubits, gen, nm, precision):

    z = np.zeros(N,dtype=complex)
    t = np.zeros(N)
    k = pow(2,gen-1)

    for i in range(N):
        wCN = 0
        qcI = QuantumCircuit(numQubits)
        qcS = QuantumCircuit(numQubits)
        measurement = QuantumCircuit(numQubits,1)

        qcI.h(0)
        qcS.h(0)
        #stateprep
        Usp(qcI)
        Usp(qcS)
        #Random Time Evoloution
        while(True): #instead of sampling from a truncated Gaussian, I will throw out a t value greater than a set maximum
            t[i] = np.random.normal(0,sigmaT,1)
            if(abs(t[i]) < T):
                break

        num = abs(k*t[i])//T
        rem = abs(k*t[i])%T
        if(t[i] >= 0):
            for j in range(int(num)):
                W(qcI,T)
                W(qcS,T)
                wCN += 2
            W(qcI,rem)
            W(qcS,rem)
            wCN += 2
        elif(t[i] < 0):
            for j in range(int(num)):
                W(qcI,-1*T)
                W(qcS,-1*T)
            W(qcI,-1*rem)
            W(qcS,-1*rem)
        t[i] = t[i]*k

        qcS.sdg(0)
        qcI.h(0)
        qcS.h(0)

        qcI = compileCT2(qcI,precision)
        qcS = compileCT2(qcS,precision)

        T_count = qcI.count_ops().get('t')+qcS.count_ops.get('t')

        measurement.barrier(range(numQubits))
        measurement.measure(0,0)

        qcI = qcI.compose(measurement)
        qcS = qcS.compose(measurement)

        aer_sim = AerSimulator(noise_model = nm)
        qcI_compiled = transpile(qcI, aer_sim)
        runs = Ns;
        job_sim = aer_sim.run(qcI_compiled, shots=runs)
        result_sim = job_sim.result()
        counts = result_sim.get_counts(qcI_compiled)

        if '0' in counts.keys():
            x = 2*counts['0']/Ns -1
        else:
            x = -1

        aer_sim = AerSimulator(noise_model = nm)
        qcS_compiled = transpile(qcS, aer_sim)
        runs = Ns;
        job_sim = aer_sim.run(qcS_compiled, shots=runs)
        result_sim = job_sim.result()
        counts = result_sim.get_counts(qcS_compiled)

        if '0' in counts.keys():
            y = 2*counts['0']/Ns -1
        else:
            y = -1

        z[i] = complex(x,y)

    return [t,z,T_count]

#N = number of data pairs per generation, number of samples per data pair, sigmaT = initial width of normal distrubtion of Hamiltonian simulation time (doubles every subsequent generation)
#T = initial maximum Hamiltonian simulation time (doubles every subsequent generation, W = function to implement time evolution operator, Usp = function to implement state prep, numQubits = number of qubits including ancilla, gens = number of generations to run
#nm = noise model, precision = rotation sysnthesis precision for each Z-rotation
def NoisyMMQcels(N, Ns, sigmaT, T, W, Usp, numQubits, gens, nm, precision):
    bnds = [(-1,1),(-1,1),(-pi,pi)]
    est = 0
    Ttot = 0
    Tmax = 0
    for i in range(gens):
        z = noisyDataGenerator(N,Ns,sigmaT,T,W,Usp,numQubits,i+1,nm,precision)
        zfit = maxF(z[1],N,z[0],est,bnds)
        est = zfit[2]
        bnds = [(-1,1),(-1,1),(est-pi/(pow(2,i)*sigmaT),est+pi/(pow(2,i)*sigmaT))]
        Ttot += Ns*z[2]
        Tmax = max(Tmax,z[2])
    return [est,Ttot,Tmax]

#




#
