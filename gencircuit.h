#pragma once

#include "qcircuit.h"

#define random(x) (rand()%x)

//
// Quantum circuit generators
//

QCircuit test();
QCircuit QFT(int numQubits);
QCircuit QFT_Quirk(int numQubits);
QCircuit QAOA(int numQubits);
QCircuit Grover(int numQubits);
QCircuit VQC(int numQubits);
QCircuit VQC1(int numQubits);
QCircuit VQC2(int numQubits);

// Generate a circuit with regular state vector using H, Z, X, CX gates
QCircuit RandomRegular(int numQubits, int numDepths);

// Generate a circuit with relatively random state vector using H, Z, X, RY, CX gates
QCircuit RandomMedium(int numQubits, int numDepths);

// Generate a circuit with random state vector using RY, CX gates
QCircuit RandomRandom(int numQubits, int numDepths);

