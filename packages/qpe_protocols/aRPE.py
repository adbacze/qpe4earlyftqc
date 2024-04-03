from Imports import *
from packages.CliffordTCompiler.CTcompiler import *
from cmath import exp

#This script provides methods to simulate adapted robust phase estimation (aRPE) as described in https://arxiv.org/pdf/2302.02454.pdf 
#Two versions of this algorithm are provided; "arpe" performs noiseless simulations without compiling circuits into Clifford+Ts while "noisy_arpe" allows for the addition
#of a noise model as well as compiles all circuits into Clifford+Ts while tracking the Tmax and Ttot.
#Since arpe does not compile circuits into Clifford+Ts it runs much faster.

#Both arpe and noisy_arpe utilize a dataGenerator function. This function performs a single generation of aRPE and returns the estimate from that generation.
#arpe and noisy_arpe track the estimates from each generation and returns the final estimate. 
#Additionally noisy_arpe also returns the Ttot and Tmax.

#Ns = Number of Samples, W = function to implement time evolution operator, Usp = function to implement state prep, numQubits = number of qubits including ancilla, gen = current generation
def dataGenerator(Ns, W, Usp, numQubits, gen): 

    qcI = QuantumCircuit(numQubits)
    qcS = QuantumCircuit(numQubits)
    measurement = QuantumCircuit(numQubits,1)

    qcI.h(0)
    qcS.h(0)

    Usp(qcI)
    Usp(qcS)

    k = pow(2,gen-1)
    for i in range(k):
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
    else: #if measurement outcome is all 1's set x = -1
        x = -1

    aer_sim = Aer.get_backend('aer_simulator')
    qcS_compiled = transpile(qcS, aer_sim)
    runs = Ns;
    job_sim = aer_sim.run(qcS_compiled, shots=runs)
    result_sim = job_sim.result()
    counts = result_sim.get_counts(qcS_compiled)

    if '0' in counts.keys():
        y = 2*counts['0']/Ns -1
    else: #if measurement outcome is all 1's set y = -1
        y = -1

    return (atan2(y,x)%(2*pi))/k #estimate is in range [0, (2*pi)/k)

#Ns = Number of Samples, Ns_final = Number of samples on final generation, W = function to implement time evolution operator, Usp = function to implement state prep, numQubits = number of qubits including ancilla, J = total generations to run
def arpe(Ns, Ns_final, W, Usp, numQubits, J):

    for i in range(J):

        if (i<(J-1)):
            argz = dataGenerator(int(Ns), W, Usp, numQubits, i+1)

        elif (i==(J-1)):
            argz = dataGenerator(int(Ns_final), W, Usp, numQubits, i+1)

        if(i==0):
            est = argz
        else:
            #shift estimate to be in correct range based off previous generation
            while ( (abs(argz - est)%(2*pi)) >  (abs((argz + 2*pi/(pow(2,i+1))) - est))%(2*pi) ):
                argz += 2*pi/(pow(2,i+1))
            est = argz #update estimate

    #Shift to (-pi,pi]
    if(est>pi):
        est = est - (2*pi)
    return est

#Ns = Number of Samples, W = function to implement time evolution operator, Usp = function to implement state prep, numQubits = number of qubits including ancilla, gen = current generation, noiseModel = qiskit noise model, precision = rotation synthesis precision for every Z-rotation
def noisyDataGenerator(Ns, W, Usp, numQubits, gen, noiseModel, precision):

    qcI = QuantumCircuit(numQubits)
    qcS = QuantumCircuit(numQubits)
    measurement = QuantumCircuit(numQubits,1)

    qcI.h(0)
    qcS.h(0)

    Usp(qcI)
    Usp(qcS)

    k = pow(2,gen-1)
    for i in range(k):
        W(qcI)
        W(qcS)

    qcS.sdg(0)
    qcI.h(0)
    qcS.h(0)

    qcI = compileCT2(qcI,precision)
    qcS = compileCT2(qcS,precision)

    T = qcI.count_ops().get('t')+qcS.count_ops().get('t')

    measurement.barrier(range(numQubits))
    measurement.measure(0,0)

    qcI = qcI.compose(measurement)
    qcS = qcS.compose(measurement)

    aer_sim = AerSimulator(noise_model = noiseModel)
    qcI_compiled = transpile(qcI, aer_sim)
    runs = Ns;
    job_sim = aer_sim.run(qcI_compiled, shots=runs)
    result_sim = job_sim.result()
    counts = result_sim.get_counts(qcI_compiled)

    if '0' in counts.keys():
        x = 2*counts['0']/Ns -1
    else: #if measurement outcome is all 1's set x = -1
        x = -1

    aer_sim = AerSimulator(noise_model = noiseModel)
    qcS_compiled = transpile(qcS, aer_sim)
    runs = Ns;
    job_sim = aer_sim.run(qcS_compiled, shots=runs)
    result_sim = job_sim.result()
    counts = result_sim.get_counts(qcS_compiled)

    if '0' in counts.keys():
        y = 2*counts['0']/Ns -1
    else: #if measurement outcome is all 1's set y = -1
        y = -1

    return [(atan2(y,x)%(2*pi))/k, T] #estimate is in range [0, (2*pi)/k), T currently does not include all the shots

#Ns = Number of Samples, Ns_final = Number of samples on final generation, W = function to implement time evolution operator, Usp = function to implement state prep, numQubits = number of qubits including ancilla, J = total generations to run, noiseModel = qiskit noise model, precision = rotation synthesis precision for every Z-rotation
def noisy_arpe(Ns, Ns_final, W, Usp, numQubits, J,  noiseModel, precision):

    Ttot = 0
    TMax = 0

    for i in range(J):

        if (i<(J-1)):
            sim = noisyDataGenerator(int(Ns), W, Usp, numQubits, i+1,  noiseModel, precision)
            argz = sim[0]
            Ttot += sim[1]*Ns #include number of shots into Ttot
            TMax = max(TMax,sim[1]/2) #divide by 2 as real and imaginary componets of Hadamard tests are separate circuits with same Ts 
        elif (i==(J-1)):
            sim = noisyDataGenerator(int(Ns_final), W, Usp, numQubits, i+1,  noiseModel, precision)
            argz = sim[0]
            Ttot += sim[1]*Ns_final #include number of shots into Ttot
            TMax = max(TMax,sim[1]/2) #divide by 2 as real and imaginary componets of Hadamard tests are separate circuits with same Ts

        if(i==0):
            est = argz
        else:
            #shift estimate to be in correct range based off previous generation
            while ( (abs(argz - est)%(2*pi)) >  (abs((argz + 2*pi/(pow(2,i+1))) - est))%(2*pi) ):
                argz += 2*pi/(pow(2,i+1))
            est = argz #update estimate
    #Shift to (-pi,pi]
    if(est>pi):
        est = est - (2*pi)
    return [est, Ttot, TMax]




#
