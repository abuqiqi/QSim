#pragma once

#include <omp.h>
#include "qcircuit.h"

// 1. 根据numTargets()分别处理
// diagonal指代是否启用diagonal gate的优化
void SVSim_targ(Matrix<DTYPE>& sv, QCircuit& qc, bool diagonal = false);
void SVSim_targ(Matrix<DTYPE>& sv, vector<QGate>& gateSeq, bool diagonal = false);

// 2. 不根据numTargets()分别处理，比collapse慢，但比没有collapse的targ快
void SVSim(Matrix<DTYPE>& sv, QCircuit& qc, bool diagonal = false);
void SVSim(Matrix<DTYPE>& sv, vector<QGate>& gateSeq, bool diagonal = false);

// Utility functions

void svsimForGate_targ(Matrix<DTYPE>& sv, QGate& gate, bool diagonal = false);
void svsimForGate(Matrix<DTYPE>& sv, QGate& gate, bool diagonal = false);

void apply1Targ(Matrix<DTYPE>& sv, QGate& gate, bool diagonal = false);
void apply2Targs(Matrix<DTYPE>& sv, QGate& gate, bool diagonal = false);
void apply3Targs(Matrix<DTYPE>& sv, QGate& gate, bool diagonal = false);
void apply4Targs(Matrix<DTYPE>& sv, QGate& gate, bool diagonal = false);
void apply5Targs(Matrix<DTYPE>& sv, QGate& gate, bool diagonal = false);
void applyMultiTargs(Matrix<DTYPE>& sv, QGate& gate, bool diagonal = false);

void applyPhase(Matrix<DTYPE>& sv, QGate& gate);
void applySwap(Matrix<DTYPE>& sv, QGate& gate);

bool isLegalControlPattern(ll ampidx, QGate& gate);
