#pragma once

#include "qcircuit.h"

/**
 * @brief State vector simulation of a gate on the state vector
 *
 * @param sv    the initial state vector
 * @param svlen the length of the state vector
 * @param qc    a quantum circuit
 */
void StateVectorSim(Matrix<DTYPE>& sv, long long svlen, QCircuit& qc);

//
// Utility functions
//

void svsim2x2GateMatrix(Matrix<DTYPE>& sv, long long svlen, QGate& gate);
void apply2x2GateMatrix(shared_ptr<Matrix<DTYPE>> gmat, DTYPE &a, DTYPE &b);

void svsim4x4GateMatrix(Matrix<DTYPE>& sv, long long svlen, QGate& gate);
void apply4x4GateMatrix(shared_ptr<Matrix<DTYPE>> gmat, DTYPE &a, DTYPE &b, DTYPE &c, DTYPE &d);
