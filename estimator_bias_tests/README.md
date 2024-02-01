This folder contains scripts used to calculate the bias on the phase estimate of RPE and QCELS.
The bias on a estimator can be defined as the error in the limit that the number of samples goes to infinity.
The bias is calculated here by increasing the number of samples until the resulting error does not decrease, that error is then the calculated bias.
Calculations were done using the 3-qubit H2 Hamilitonian, however, the estimator bias should be indepenet of the problem instance and only dependent on the target state overlap.

These scripts should be run from the parent directory.
