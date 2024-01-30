from Imports import *
from scipy.optimize import minimize, shgo, least_squares, differential_evolution
from cmath import exp

def qcels(N, Ns, W, Usp, numQubits, repetitions,target):

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
                W(qcI,target)
                W(qcS,target)

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

        z[l] = complex(x,y)
    return z

def F(x, z, N, t):
    Z_fit=np.zeros(N,dtype = 'complex_')
    Z_fit=(x[0]+1j*x[1])*np.exp(-1j*x[2]*t)
    return (np.linalg.norm(Z_fit-z)**2/N)

def maxF(z, N, t, est, bnds):
    fun = lambda x: F(x, z, N, t)
    result = minimize(fun,x0=[0.5,0,est],method = 'SLSQP',bounds=bnds)
    return result.x


def noisy_qcels(N, Ns, W, Usp, numQubits, repetitions, nm, precision, target):

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
                W(qcI,target)
                W(qcS,target)

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
    return [z, T]


#
