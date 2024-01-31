from Imports import *
from scipy.optimize import minimize
from cmath import exp

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

    return [t,z]

def F(x, z, N, t):
    Z_fit=np.zeros(N,dtype = 'complex_')
    Z_fit=(x[0]+1j*x[1])*np.exp(-1j*x[2]*t)
    return (np.linalg.norm(Z_fit-z)**2/N)

def maxF(z, N, t, est, bnds):
    fun = lambda x: F(x, z, N, t)
    result = minimize(fun,x0=[0.5,0,est],method = 'SLSQP',bounds=bnds)
    return result.x

def MMQcels(N, Ns, sigmaT, T, W, Usp, numQubits, gens):
    bnds = [(-1,1),(-1,1),(-pi,pi)]
    est = 0
    for i in range(gens):
        z = dataGenerator(N,Ns,sigmaT,T,W,Usp,numQubits,i+1)
        zfit = maxF(z[1],N,z[0],est,bnds)
        est = zfit[2]
        bnds = [(-1,1),(-1,1),(est-pi/(pow(2,i)*sigmaT),est+pi/(pow(2,i)*sigmaT))
    return est


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

        T = qcI.count_ops().get('t')+qcS.count_ops.get('t')

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

    return [t,z,T]

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
