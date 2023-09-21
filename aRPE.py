from Imports import *
from cmath import exp

def dataGenerator(Ns, W, Usp, numQubits, gen, noiseModel):

    Tmax = 0
    T = 0
    k = pow(2,gen-1)
    wCalls = 0

    qcI = QuantumCircuit(numQubits)
    qcS = QuantumCircuit(numQubits)
    measurement = QuantumCircuit(numQubits,1)

    qcI.h(0)
    qcS.h(0)

    Usp(qcI)
    Usp(qcS)

    for i in range(k):
        W(qcI)
        W(qcS)
        wCalls += 2

    qcS.sdg(0)
    qcI.h(0)
    qcS.h(0)

    #qcI = compileCT(qcI,1e-4)
    #qcS = compileCT(qcS,1e-4)

    measurement.barrier(range(numQubits))
    measurement.measure(0,0)

    qcI = qcI.compose(measurement)
    qcS = qcS.compose(measurement)

    #numD1 = qcI.count_ops().get('t')
    #numD2 = qcS.count_ops().get('t')

    #T = numD1 + numD2;
    #Tmax = max(Tmax,numD1,numD2);
    #print(qcI.draw())

    aer_sim = AerSimulator(noise_model = noiseModel)
    qcI_compiled = transpile(qcI, aer_sim)
    runs = Ns;
    #print(runs)
    job_sim = aer_sim.run(qcI_compiled, shots=runs)
    result_sim = job_sim.result()
    counts = result_sim.get_counts(qcI_compiled)

    if '0' in counts.keys():
        x = 2*counts['0']/Ns -1
    else:
        x = -1


    aer_sim = AerSimulator(noise_model = noiseModel)
    qcS_compiled = transpile(qcS, aer_sim)
    runs = Ns;
    job_sim = aer_sim.run(qcS_compiled, shots=runs)
    result_sim = job_sim.result()
    counts = result_sim.get_counts(qcS_compiled)

    #print(sum(counts.values()))
    if '0' in counts.keys():
        y = 2*counts['0']/Ns -1
    else:
        y = -1

    return [(atan2(y,x)%(2*pi))/k, wCalls]

def arpe(Ns, Ns_final, W, Usp, numQubits, J,  noiseModel):

    wCalls = 0
    wCallsMax = 0

    for i in range(J):

        if (i<(J-1)):
            sim = dataGenerator(int(Ns), W, Usp, numQubits, i+1,  noiseModel)
            argz = sim[0]
            wCalls += sim[1]*Ns
            wCallsMax = max(wCallsMax,sim[1]/2)
        elif (i==(J-1)):
            sim = dataGenerator(int(Ns_final), W, Usp, numQubits, i+1,  noiseModel)
            argz = sim[0]
            wCalls += sim[1]*Ns_final
            wCallsMax = max(wCallsMax,sim[1]/2)

        if(i==0):
            est = argz
        else:
            while ( (abs(argz - est)) >  (abs((argz + 2*pi/(pow(2,i+1))) - est)) ):
                #print("hi")
                argz += 2*pi/(pow(2,i+1))
            est = argz
            #print(est)
    #Shift to (-pi,pi]
    if(est>pi):
        est = est - (2*pi)
    return [est, wCalls, wCallsMax]


#
