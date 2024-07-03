#pragma once

#include "qcircuit.h"

/**
 * @brief State vector simulation of a quantum circuit on the state vector
 *
 * @param sv the state vector
 * @param qc a quantum circuit
 */
void StateVectorSim(Matrix<DTYPE>& sv, QCircuit& qc);

//
// Utility functions
//

/**
 * @brief State vector simulation for a quantum gate
 * 
 * @param sv   the state vector
 * @param gate the processing gate
 */
void svsimForGate(Matrix<DTYPE>& sv, QGate& gate);

/**
 * @brief Check if the index of a amplitude contains a legal control mask of the gate
 * 
 * @param amp  the amplitude index
 * @param gate the processing gate
 * @return int 0: illegal control mask; 1: legal control mask
 */
int checkControlMask(ll amp, QGate& gate);