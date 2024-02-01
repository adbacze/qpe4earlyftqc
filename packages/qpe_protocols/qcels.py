from Imports import *
from scipy.optimize import minimize
from cmath import exp
from packages.CliffordTCompiler.CTcompiler import *

#This script provides methods to simulate Quantum Complex Exponential Least Squares (QCELS) as described in https://arxiv.org/pdf/2211.11973.pdf
#Two versions of this algorithm are provided; "qcels" performs noiseless simulations without compiling circuits into Clifford+Ts while "noisy_qcels" allows for the addition
#of a noise model as well as compiles all circuits into Clifford+Ts while tracking the Tmax and Ttot.
#Since qcels does not compile circuits into Clifford+Ts it runs much faster.

#Both qcels and noisy_qcels actually perform just a single generation of the algorithm.
#Separate functions are provide to perform the least squares fit

#N = Number of data pairs, Ns = Number of Samples, W = function to implement time evolution operator, Usp = function to implement state prep, numQubits = number of qubits including ancilla, repetitions = number of time evolution operator applications per data pair
def qcels(N, Ns, W, Usp, numQubits, repetitions):

    z = np.zeros(N,dtype=complex)

    for l in range(N):
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
            x = -1 #If measurement outcomes are all 1 the set x = -1

        aer_sim = Aer.get_backend('aer_simulator')
        qcS_compiled = transpile(qcS, aer_sim)
        runs = Ns;
        job_sim = aer_sim.run(qcS_compiled, shots=runs)
        result_sim = job_sim.result()
        counts = result_sim.get_counts(qcS_compiled)

        if '0' in counts.keys():
            y = 2*counts['0']/Ns -1
        else:
            y = -1 #If measurement outcomes are all 1 the set y = -1

        z[l] = complex(x,y) #return time series  
    return z

def F(x, z, N, t): #Define function to minimized for least squares fir 
    Z_fit=np.zeros(N,dtype = 'complex_')
    Z_fit=(x[0]+1j*x[1])*np.exp(-1j*x[2]*t)
    return (np.linalg.norm(Z_fit-z)**2/N)

#z = time series data generated form qcels, N = number of data pairs, t = corresponding times for time series, est = prior estimate, bnds = bounds for fit
def maxF(z, N, t, est, bnds): #performs least squares fit and returns eigenvalue estimates
    fun = lambda x: F(x, z, N, t)
    result = minimize(fun,x0=[0.5,0,est],method = 'SLSQP',bounds=bnds)
    return result.x


#N = Number of data pairs, Ns = Number of Samples, W = function to implement time evolution operator, Usp = function to implement state prep, numQubits = number of qubits including ancilla, repetitions = number of time evolution operator applications per data pair, precision = rotation synthesis 
def noisy_qcels(N, Ns, W, Usp, numQubits, repetitions, nm, precision):

    z = np.zeros(N,dtype=complex)

    for l in range(N):
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

        qcI = compileCT2(qcI,precision)
        qcS = compileCT2(qcS,precision)
        numD1 = qcI.count_ops().get('t')
        numD2 = qcS.count_ops().get('t')
        T = numD1 + numD2

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

        z[l] = complex(x,y)
    return [z, T] #returns time series along with the T count for a single shot of a generation of qcels


#
