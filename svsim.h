#pragma once

#include <omp.h>
#include "qcircuit.h"

/**
 * @brief State vector simulation of a quantum circuit on the state vector
 *
 * @param sv the state vector
 * @param qc a quantum circuit
 */
void SVSim(Matrix<DTYPE>& sv, QCircuit& qc);
void SVSim(Matrix<DTYPE>& sv, vector<QGate>& gateSeq);
void SVSimNoFuser(Matrix<DTYPE>& sv, QCircuit& qc);
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
void svsimForGateNoFuser(Matrix<DTYPE>& sv, QGate& gate);
void applyPhase(Matrix<DTYPE>& sv, QGate& gate);
void applyDiagonal(Matrix<DTYPE>& sv, QGate& gate);
void apply1Targ(Matrix<DTYPE>& sv, QGate& gate);
void applySwap(Matrix<DTYPE>& sv, QGate& gate);
void apply2Targs(Matrix<DTYPE>& sv, QGate& gate);
void apply3Targs(Matrix<DTYPE>& sv, QGate& gate);
void apply4Targs(Matrix<DTYPE>& sv, QGate& gate);
void apply5Targs(Matrix<DTYPE>& sv, QGate& gate);
void applyMultiTargs(Matrix<DTYPE>& sv, QGate& gate);

/**
 * @brief Check if the index of an amplitude is a legal control pattern of the gate
 * 
 * @param ampidx the amplitude index
 * @param gate the processing gate
 * @return true  ampidx is a legal control pattern
 * @return false ampidx is an illegal control pattern
 */
bool isLegalControlPattern(ll ampidx, QGate& gate);