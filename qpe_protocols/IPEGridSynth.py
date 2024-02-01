from Imports import *

#Compiling IPE circuits into Clifford+Ts isn't straightforward due to the measurement and classically controlled operations
#
def gsIPE(Ns, W, Usp, numQubits, digits, precision, nm):

    q = QuantumRegister(numQubits)
    c = ClassicalRegister(digits)
    ipe = QuantumCircuit(q,c)
    prepCirc = QuantumCircuit(numQubits)

    Usp(prepCirc)
    prepCirc = compileCT2(prepCirc, precision)
    ipe = ipe.compose(prepCirc)

    for i in range(digits):
        wCirc = QuantumCircuit(numQubits)
        ipe.h(0)
        for j in range(pow(2,digits-i-1)):
            W(wCirc,precision)

        #wCirc = compileCT2(wCirc,precision)
        ipe = ipe.compose(wCirc)

        for l in range(pow(2,i)-1):
            ipe.rz(-((pi)/pow(2,i))*float(l+1),0).c_if(c,l+1)
            #epsilonClCRZ(ipe,0,-((pi)/pow(2,i))*float(l+1),0.9*precision,c,l+1) #Synthesizes classically controlled rotation

        ipe.h(0)
        ipe.measure(0,i)
        ipe.reset(0)

    #T = ipe.count_ops().get('t')
    #print(ipe.draw())

    aer_sim = AerSimulator(noise_model = nm)
    ipe_compiled = transpile(ipe, aer_sim)
    runs = Ns;
    job_sim = aer_sim.run(ipe_compiled, shots=runs)
    result_sim = job_sim.result()
    counts = result_sim.get_counts(ipe_compiled)

    #Convert results from binary representation to decimal
    keyDec = [str(int(key,2)/pow(2,digits)) for key in list(counts.keys())]
    countsDec = dict(zip(keyDec, counts.values()))

    est = float(max(countsDec, key=countsDec.get)) #majority vote
    if(est<=0.5):
        est = est*2*pi
    elif(est>0.5):
        est = est*2*pi - 2*pi

    return [est, wCalls]














#
