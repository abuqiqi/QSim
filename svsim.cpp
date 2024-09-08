#include "svsim.h"

#define OMP_ENABLED

void SVSimNoFuser(Matrix<DTYPE>& sv, QCircuit& qc) {
    for (int j = 0; j < qc.numDepths; ++ j) {
        for (int qid = 0; qid < qc.numQubits; ++ qid) {
            QGate& gate = qc.gates[j][qid];
            if (gate.isIDE() || gate.isMARK()) {
                continue;
            }
            // gate.print();
            svsimForGateNoFuser(sv, gate);
            // sv.print();
        }
    }
}

/**
 * @brief State vector simulation of a quantum circuit on the state vector
 * 
 * @param sv    the initial state vector
 * @param qc    a quantum circuit
 */
void SVSim(Matrix<DTYPE>& sv, QCircuit& qc) {
    for (int j = 0; j < qc.numDepths; ++ j) {
        for (int qid = 0; qid < qc.numQubits; ++ qid) {
            QGate& gate = qc.gates[j][qid];
            if (gate.isIDE() || gate.isMARK()) {
                continue;
            }
            // gate.print();
            svsimForGate(sv, gate);
            // sv.print();
        }
    }
}

void SVSim(Matrix<DTYPE>& sv, vector<QGate>& gateSeq) {
    for (auto& gate : gateSeq) {
        if (gate.isIDE() || gate.isMARK()) {
            continue;
        }
        // cout << "[DEBUG] Processing gate: " << gate.gname << endl;
        // gate.printInfo();
        svsimForGate(sv, gate);
    }
}

void svsimForGateNoFuser(Matrix<DTYPE>& sv, QGate& gate) {
    if (gate.numTargets() == 1) {
        apply1Targ(sv, gate);
    } else if (gate.gname == "SWAP") {
        applySwap(sv, gate);
    } else if (gate.numTargets() == 2) {
        apply2Targs(sv, gate);
    } else if (gate.numTargets() == 3) {
        apply3Targs(sv, gate);
    } else if (gate.numTargets() == 4) {
        apply4Targs(sv, gate);
    } else if (gate.numTargets() == 5) {
        apply5Targs(sv, gate);
    } else {
        applyMultiTargs(sv, gate);
    }
}

/**
 * @brief State vector simulation for a quantum gate
 * 
 * @param sv    the state vector
 * @param gate  the processing gate
 */
void svsimForGate(Matrix<DTYPE>& sv, QGate& gate) {
    if (gate.isPhase()) {
        applyPhase(sv, gate);
    } else if (gate.numTargets() == 1) {
        apply1Targ(sv, gate);
    } else if (gate.gname == "SWAP") {
        applySwap(sv, gate);
    } else if (gate.numTargets() == 2) {
        apply2Targs(sv, gate);
    } else if (gate.numTargets() == 3) {
        apply3Targs(sv, gate);
    } else if (gate.numTargets() == 4) {
        apply4Targs(sv, gate);
    } else if (gate.numTargets() == 5) {
        apply5Targs(sv, gate);
    } else {
        applyMultiTargs(sv, gate);
    }
    // applyMultiTargs(sv, gate);
}

void applyPhase(Matrix<DTYPE>& sv, QGate& gate) {
    DTYPE phase = gate.gmat->data[gate.gmat->row-1][gate.gmat->col-1];
    ll pattern = 0; // pattern: ...11...1...
    for (int i = 0; i < gate.numQubits(); ++ i) {
        pattern |= (1 << gate.qubits()[i]);
    }
    for (ll i = 0; i < sv.row; ++ i) {
        if ((i & pattern) == pattern) {
            sv.data[i][0] *= phase;
        }
    }
}

void apply1Targ(Matrix<DTYPE>& sv, QGate& gate) {
    int qid = gate.targetQubits[0];
#ifdef OMP_ENABLED
#pragma omp parallel for // collapse(2)
#endif
    for (ll i = 0; i < sv.row; i += (1<<(qid+1))) {
        for (ll j = 0; j < (1<<qid); ++ j) {
            auto p = i | j;
            if (! isLegalControlPattern(p, gate))
                continue;
            auto q0 = sv.data[p][0];
            auto q1 = sv.data[p|(1<<qid)][0];
            if (gate.isDiagonal) {
                sv.data[p][0] = gate.gmat->data[0][0] * q0;
                sv.data[p|(1<<qid)][0] = gate.gmat->data[1][1] * q1;
            } else {
                sv.data[p][0] = gate.gmat->data[0][0] * q0 + gate.gmat->data[0][1] * q1;
                sv.data[p|(1<<qid)][0] = gate.gmat->data[1][0] * q0 + gate.gmat->data[1][1] * q1;
            }
        }
    }
}

void applySwap(Matrix<DTYPE>& sv, QGate& gate) {
    int q0 = gate.targetQubits[0];
    int q1 = gate.targetQubits[1];
#ifdef OMP_ENABLED
#pragma omp parallel for // collapse(3)
#endif
    for (ll i = 0; i < sv.row; i += (1<<(q1+1))) {
        for (ll j = 0; j < (1<<q1); j += (1<<(q0+1))) {
            for (ll k = 0; k < (1<<q0); ++ k) {
                auto p = i | j | k;
                if (! isLegalControlPattern(p, gate))
                    continue;
                auto amp1 = sv.data[p|(1<<q0)][0];
                sv.data[p|(1<<q0)][0] = sv.data[p|(1<<q1)][0];
                sv.data[p|(1<<q1)][0] = amp1;
            }
        }
    }
}

// void apply2Targs(Matrix<DTYPE>& sv, QGate& gate) {
//     int q0 = gate.targetQubits[0];
//     int q1 = gate.targetQubits[1];
//     vector<ll> masks;
//     int numAmps = 4;
//     masks.push_back(0);
//     masks.push_back(1<<q0);
//     masks.push_back(1<<q1);
//     masks.push_back((1<<q0)|(1<<q1));
// #pragma omp parallel for collapse(3)
//     for (ll i = 0; i < sv.row; i += (1<<(q1+1))) {
//         for (ll j = 0; j < (1<<q1); j += (1<<(q0+1))) {
//             for (ll k = 0; k < (1<<q0); ++ k) {
//                 auto p = i | j | k;
//                 if (! isLegalControlPattern(p, gate))
//                     continue;
//                 vector<DTYPE> amps_vec(numAmps); // save the involved amplitudes
//                 for (int idx = 0; idx < numAmps; ++ idx) {
//                     amps_vec[idx] = sv.data[p | masks[idx]][0];
//                 }
//                 if (gate.isDiagonal) {
//                     for (int idx = 0; idx < numAmps; ++ idx) {
//                         sv.data[p | masks[idx]][0] = gate.gmat->data[idx][idx] * amps_vec[idx];
//                     }
//                 } else {
//                     for (int i = 0; i < numAmps; ++ i) {
//                         sv.data[p | masks[i]][0] = 0;
//                         for (int j = 0; j < numAmps; ++ j) {
//                             sv.data[p | masks[i]][0] += gate.gmat->data[i][j] * amps_vec[j];
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }

void apply2Targs(Matrix<DTYPE>& sv, QGate& gate) {
    int q0 = gate.targetQubits[0];
    int q1 = gate.targetQubits[1];
#ifdef OMP_ENABLED
#pragma omp parallel for // collapse(3)
#endif
    for (ll i = 0; i < sv.row; i += (1<<(q1+1))) {
        for (ll j = 0; j < (1<<q1); j += (1<<(q0+1))) {
            for (ll k = 0; k < (1<<q0); ++ k) {
                auto p = i | j | k;
                if (! isLegalControlPattern(p, gate))
                    continue;
                if (gate.isDiagonal) {
                    sv.data[p][0] *= gate.gmat->data[0][0];
                    sv.data[p|(1<<q0)][0] *= gate.gmat->data[1][1];
                    sv.data[p|(1<<q1)][0] *= gate.gmat->data[2][2];
                    sv.data[p|(1<<q0)|(1<<q1)][0] *= gate.gmat->data[3][3];
                } else {
                    auto q00 = sv.data[p][0];
                    auto q01 = sv.data[p|(1<<q0)][0];
                    auto q10 = sv.data[p|(1<<q1)][0];
                    auto q11 = sv.data[p|(1<<q0)|(1<<q1)][0];
                    sv.data[p][0] = gate.gmat->data[0][0] * q00 + gate.gmat->data[0][1] * q01 + gate.gmat->data[0][2] * q10 + gate.gmat->data[0][3] * q11;
                    sv.data[p|(1<<q0)][0] = gate.gmat->data[1][0] * q00 + gate.gmat->data[1][1] * q01 + gate.gmat->data[1][2] * q10 + gate.gmat->data[1][3] * q11;
                    sv.data[p|(1<<q1)][0] = gate.gmat->data[2][0] * q00 + gate.gmat->data[2][1] * q01 + gate.gmat->data[2][2] * q10 + gate.gmat->data[2][3] * q11;
                    sv.data[p|(1<<q0)|(1<<q1)][0] = gate.gmat->data[3][0] * q00 + gate.gmat->data[3][1] * q01 + gate.gmat->data[3][2] * q10 + gate.gmat->data[3][3] * q11;
                }
            }
        }
    }
}

// void apply3Targs(Matrix<DTYPE>& sv, QGate& gate) {
//     int q0 = gate.targetQubits[0];
//     int q1 = gate.targetQubits[1];
//     int q2 = gate.targetQubits[2];
//     vector<ll> masks;
//     int numAmps = 8;
//     masks.push_back(0);
//     masks.push_back(1<<q0);
//     masks.push_back(1<<q1);
//     masks.push_back((1<<q0)|(1<<q1));
//     masks.push_back(1<<q2);
//     masks.push_back((1<<q0)|(1<<q2));
//     masks.push_back((1<<q1)|(1<<q2));
//     masks.push_back((1<<q0)|(1<<q1)|(1<<q2));
// #pragma omp parallel for collapse(4)
//     for (ll i = 0; i < sv.row; i += (1<<(q2+1))) {
//         for (ll j = 0; j < (1<<q2); j += (1<<(q1+1))) {
//             for (ll k = 0; k < (1<<q1); k += (1<<(q0+1))) {
//                 for (ll l = 0; l < (1<<q0); ++ l) {
//                     auto p = i | j | k | l;
//                     if (! isLegalControlPattern(p, gate))
//                         continue;
//                     vector<DTYPE> amps_vec(numAmps); // save the involved amplitudes
//                     for (int idx = 0; idx < numAmps; ++ idx) {
//                         amps_vec[idx] = sv.data[p | masks[idx]][0];
//                     }
//                     if (gate.isDiagonal) {
//                         for (int idx = 0; idx < numAmps; ++ idx) {
//                             sv.data[p | masks[idx]][0] = gate.gmat->data[idx][idx] * amps_vec[idx];
//                         }
//                     } else {
//                         // for (int i = 0; i < numAmps; ++ i) {
//                         //     sv.data[p | masks[i]][0] = 0;
//                         //     for (int j = 0; j < numAmps; ++ j) {
//                         //         sv.data[p | masks[i]][0] += gate.gmat->data[i][j] * amps_vec[j];
//                         //     }
//                         // }
//                         sv.data[p | masks[0]][0] = gate.gmat->data[0][0] * amps_vec[0] + gate.gmat->data[0][1] * amps_vec[1] + gate.gmat->data[0][2] * amps_vec[2] + gate.gmat->data[0][3] * amps_vec[3];
//                     }
//                 }
//             }
//         }
//     }
// }

void apply3Targs(Matrix<DTYPE>& sv, QGate& gate) {
    int q0 = gate.targetQubits[0];
    int q1 = gate.targetQubits[1];
    int q2 = gate.targetQubits[2];
#ifdef OMP_ENABLED
#pragma omp parallel for // collapse(4)
#endif
    for (ll i = 0; i < sv.row; i += (1<<(q2+1))) {
        for (ll j = 0; j < (1<<q2); j += (1<<(q1+1))) {
            for (ll k = 0; k < (1<<q1); k += (1<<(q0+1))) {
                for (ll l = 0; l < (1<<q0); ++ l) {
                    auto p = i | j | k | l;
                    if (! isLegalControlPattern(p, gate))
                        continue;
                    if (gate.isDiagonal) {
                        sv.data[p][0] *= gate.gmat->data[0][0];
                        sv.data[p|(1<<q0)][0] *= gate.gmat->data[1][1];
                        sv.data[p|(1<<q1)][0] *= gate.gmat->data[2][2];
                        sv.data[p|(1<<q0)|(1<<q1)][0] *= gate.gmat->data[3][3];
                        sv.data[p|(1<<q2)][0] *= gate.gmat->data[4][4];
                        sv.data[p|(1<<q0)|(1<<q2)][0] *= gate.gmat->data[5][5];
                        sv.data[p|(1<<q1)|(1<<q2)][0] *= gate.gmat->data[6][6];
                        sv.data[p|(1<<q0)|(1<<q1)|(1<<q2)][0] *= gate.gmat->data[7][7];
                    } else {
                        auto q000 = sv.data[p][0];
                        auto q001 = sv.data[p|(1<<q0)][0];
                        auto q010 = sv.data[p|(1<<q1)][0];
                        auto q011 = sv.data[p|(1<<q0)|(1<<q1)][0];
                        auto q100 = sv.data[p|(1<<q2)][0];
                        auto q101 = sv.data[p|(1<<q0)|(1<<q2)][0];
                        auto q110 = sv.data[p|(1<<q1)|(1<<q2)][0];
                        auto q111 = sv.data[p|(1<<q0)|(1<<q1)|(1<<q2)][0];
                        sv.data[p][0] = gate.gmat->data[0][0] * q000 + gate.gmat->data[0][1] * q001 + gate.gmat->data[0][2] * q010 + gate.gmat->data[0][3] * q011 + gate.gmat->data[0][4] * q100 + gate.gmat->data[0][5] * q101 + gate.gmat->data[0][6] * q110 + gate.gmat->data[0][7] * q111;
                        sv.data[p|(1<<q0)][0] = gate.gmat->data[1][0] * q000 + gate.gmat->data[1][1] * q001 + gate.gmat->data[1][2] * q010 + gate.gmat->data[1][3] * q011 + gate.gmat->data[1][4] * q100 + gate.gmat->data[1][5] * q101 + gate.gmat->data[1][6] * q110 + gate.gmat->data[1][7] * q111;
                        sv.data[p|(1<<q1)][0] = gate.gmat->data[2][0] * q000 + gate.gmat->data[2][1] * q001 + gate.gmat->data[2][2] * q010 + gate.gmat->data[2][3] * q011 + gate.gmat->data[2][4] * q100 + gate.gmat->data[2][5] * q101 + gate.gmat->data[2][6] * q110 + gate.gmat->data[2][7] * q111;
                        sv.data[p|(1<<q0)|(1<<q1)][0] = gate.gmat->data[3][0] * q000 + gate.gmat->data[3][1] * q001 + gate.gmat->data[3][2] * q010 + gate.gmat->data[3][3] * q011 + gate.gmat->data[3][4] * q100 + gate.gmat->data[3][5] * q101 + gate.gmat->data[3][6] * q110 + gate.gmat->data[3][7] * q111;
                        sv.data[p|(1<<q2)][0] = gate.gmat->data[4][0] * q000 + gate.gmat->data[4][1] * q001 + gate.gmat->data[4][2] * q010 + gate.gmat->data[4][3] * q011 + gate.gmat->data[4][4] * q100 + gate.gmat->data[4][5] * q101 + gate.gmat->data[4][6] * q110 + gate.gmat->data[4][7] * q111;
                        sv.data[p|(1<<q0)|(1<<q2)][0] = gate.gmat->data[5][0] * q000 + gate.gmat->data[5][1] * q001 + gate.gmat->data[5][2] * q010 + gate.gmat->data[5][3] * q011 + gate.gmat->data[5][4] * q100 + gate.gmat->data[5][5] * q101 + gate.gmat->data[5][6] * q110 + gate.gmat->data[5][7] * q111;
                        sv.data[p|(1<<q1)|(1<<q2)][0] = gate.gmat->data[6][0] * q000 + gate.gmat->data[6][1] * q001 + gate.gmat->data[6][2] * q010 + gate.gmat->data[6][3] * q011 + gate.gmat->data[6][4] * q100 + gate.gmat->data[6][5] * q101 + gate.gmat->data[6][6] * q110 + gate.gmat->data[6][7] * q111;
                        sv.data[p|(1<<q0)|(1<<q1)|(1<<q2)][0] = gate.gmat->data[7][0] * q000 + gate.gmat->data[7][1] * q001 + gate.gmat->data[7][2] * q010 + gate.gmat->data[7][3] * q011 + gate.gmat->data[7][4] * q100 + gate.gmat->data[7][5] * q101 + gate.gmat->data[7][6] * q110 + gate.gmat->data[7][7] * q111;
                    }
                }
            }
        }
    }
}

void apply4Targs(Matrix<DTYPE>& sv, QGate& gate) {
    int q0 = gate.targetQubits[0];
    int q1 = gate.targetQubits[1];
    int q2 = gate.targetQubits[2];
    int q3 = gate.targetQubits[3];

    vector<ll> masks;
    int numAmps = 16;
    masks.push_back(0);
    masks.push_back(1<<q0);
    masks.push_back(1<<q1);
    masks.push_back((1<<q0)|(1<<q1));
    masks.push_back(1<<q2);
    masks.push_back((1<<q0)|(1<<q2));
    masks.push_back((1<<q1)|(1<<q2));
    masks.push_back((1<<q0)|(1<<q1)|(1<<q2));
    masks.push_back(1<<q3);
    masks.push_back((1<<q0)|(1<<q3));
    masks.push_back((1<<q1)|(1<<q3));
    masks.push_back((1<<q0)|(1<<q1)|(1<<q3));
    masks.push_back((1<<q2)|(1<<q3));
    masks.push_back((1<<q0)|(1<<q2)|(1<<q3));
    masks.push_back((1<<q1)|(1<<q2)|(1<<q3));
    masks.push_back((1<<q0)|(1<<q1)|(1<<q2)|(1<<q3));

#ifdef OMP_ENABLED
#pragma omp parallel for // collapse(5)
#endif
    for (ll i = 0; i < sv.row; i += (1<<(q3+1))) {
        for (ll j = 0; j < (1<<q3); j += (1<<(q2+1))) {
            for (ll k = 0; k < (1<<q2); k += (1<<(q1+1))) {
                for (ll l = 0; l < (1<<q1); l += (1<<(q0+1))) {
                    for (ll m = 0; m < (1<<q0); ++m) {
                        auto p = i | j | k | l | m;
                        if (!isLegalControlPattern(p, gate))
                            continue;
                        vector<DTYPE> amps_vec(numAmps); // save the involved amplitudes
                        for (int idx = 0; idx < numAmps; ++ idx) {
                            amps_vec[idx] = sv.data[p | masks[idx]][0];
                        }
                        if (gate.isDiagonal) {
                            for (int idx = 0; idx < numAmps; ++ idx) {
                                sv.data[p | masks[idx]][0] = gate.gmat->data[idx][idx] * amps_vec[idx];
                            }
                        } else {
                            for (int i = 0; i < numAmps; ++ i) {
                                sv.data[p | masks[i]][0] = 0;
                                for (int j = 0; j < numAmps; ++ j) {
                                    sv.data[p | masks[i]][0] += gate.gmat->data[i][j] * amps_vec[j];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void apply5Targs(Matrix<DTYPE>& sv, QGate& gate) {
    int q0 = gate.targetQubits[0];
    int q1 = gate.targetQubits[1];
    int q2 = gate.targetQubits[2];
    int q3 = gate.targetQubits[3];
    int q4 = gate.targetQubits[4];

    vector<ll> masks;
    int numAmps = 32;
    masks.push_back(0);
    masks.push_back(1<<q0);
    masks.push_back(1<<q1);
    masks.push_back((1<<q0)|(1<<q1));
    masks.push_back(1<<q2);
    masks.push_back((1<<q0)|(1<<q2));
    masks.push_back((1<<q1)|(1<<q2));
    masks.push_back((1<<q0)|(1<<q1)|(1<<q2));
    masks.push_back(1<<q3);
    masks.push_back((1<<q0)|(1<<q3));
    masks.push_back((1<<q1)|(1<<q3));
    masks.push_back((1<<q0)|(1<<q1)|(1<<q3));
    masks.push_back((1<<q2)|(1<<q3));
    masks.push_back((1<<q0)|(1<<q2)|(1<<q3));
    masks.push_back((1<<q1)|(1<<q2)|(1<<q3));
    masks.push_back((1<<q0)|(1<<q1)|(1<<q2)|(1<<q3));
    masks.push_back(1<<q4);
    masks.push_back((1<<q0)|(1<<q4));
    masks.push_back((1<<q1)|(1<<q4));
    masks.push_back((1<<q0)|(1<<q1)|(1<<q4));
    masks.push_back((1<<q2)|(1<<q4));
    masks.push_back((1<<q0)|(1<<q2)|(1<<q4));
    masks.push_back((1<<q1)|(1<<q2)|(1<<q4));
    masks.push_back((1<<q0)|(1<<q1)|(1<<q2)|(1<<q4));
    masks.push_back((1<<q3)|(1<<q4));
    masks.push_back((1<<q0)|(1<<q3)|(1<<q4));
    masks.push_back((1<<q1)|(1<<q3)|(1<<q4));
    masks.push_back((1<<q0)|(1<<q1)|(1<<q3)|(1<<q4));
    masks.push_back((1<<q2)|(1<<q3)|(1<<q4));
    masks.push_back((1<<q0)|(1<<q2)|(1<<q3)|(1<<q4));
    masks.push_back((1<<q1)|(1<<q2)|(1<<q3)|(1<<q4));
    masks.push_back((1<<q0)|(1<<q1)|(1<<q2)|(1<<q3)|(1<<q4));

#ifdef OMP_ENABLED
#pragma omp parallel for // collapse(6)
#endif
    for (ll i = 0; i < sv.row; i += (1<<(q4+1))) {
        for (ll j = 0; j < (1<<q4); j += (1<<(q3+1))) {
            for (ll k = 0; k < (1<<q3); k += (1<<(q2+1))) {
                for (ll l = 0; l < (1<<q2); l += (1<<(q1+1))) {
                    for (ll m = 0; m < (1<<q1); m += (1<<(q0+1))) {
                        for (ll n = 0; n < (1<<q0); ++n) {
                            auto p = i | j | k | l | m | n;
                            if (!isLegalControlPattern(p, gate))
                                continue;
                            vector<DTYPE> amps_vec(numAmps); // save the involved amplitudes
                            for (int idx = 0; idx < numAmps; ++ idx) {
                                amps_vec[idx] = sv.data[p | masks[idx]][0];
                            }
                            if (gate.isDiagonal) {
                                for (int idx = 0; idx < numAmps; ++ idx) {
                                    sv.data[p | masks[idx]][0] = gate.gmat->data[idx][idx] * amps_vec[idx];
                                }
                            } else {
                                for (int i = 0; i < numAmps; ++ i) {
                                    sv.data[p | masks[i]][0] = 0;
                                    for (int j = 0; j < numAmps; ++ j) {
                                        sv.data[p | masks[i]][0] += gate.gmat->data[i][j] * amps_vec[j];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// void applyMultiTargs(Matrix<DTYPE>& sv, QGate& gate) {
//     bool isAccessed[sv.row];
//     memset(isAccessed, 0, sv.row*sizeof(bool));

//     ll numAmps = (1 << gate.numTargets()); // the number of amplitudes involved in matrix-vector multiplication
//     Matrix<DTYPE> amps_vec(numAmps, 1); // save the involved amplitudes

//     // 1. Calculate the strides for the involved amplitudes
//     /*  <e.g.>
//         (1) If there is only one target qubit, two amplitudes are involved. 
//             If the target qubit is q[k], strides = {0, 2^k}. 
//         (2) If there are two target qubits, four amplitudes are involved. 
//                       c  t
//             idx[0] |..0..0..> -+ strides[0] = 0 --+------------------+
//                                |                  |                  |
//             idx[1] |..0..1..> -+ strides[1] = 2^t |                  |
//                                                   |                  |
//             idx[2] |..1..0..> --------------------+ strides[2] = 2^c |
//                                                                      |
//             idx[3] |..1..1..> ---------------------------------------+ strides[3] = 2^c + 2^t
//         (3) If there are 'gate.numTargets()' target qubits, there are 'numAmps' amplitudes. 
//     */
//     vector<ll> strides(numAmps, 0);
//     for (ll idx = 0; idx < numAmps; ++ idx) {
//         ll stride = 0;
//         for (int j = 0; j < gate.numTargets(); ++ j) {
//             if (idx & (1 << j)) { // if the j-th bit of idx is 1
//                 stride += (1 << gate.targetQubits[j]);
//             }
//         }
//         strides[idx] = stride;
//     }

//     // 2. Iterate over all amplitudes
//     for (ll ampidx = 0; ampidx < sv.row; ++ ampidx) {
//         // 2.1. Skip the amplitude if it is already accessed
//         if (isAccessed[ampidx]) continue;

//         // 2.2. Save the involved amplitudes to amps_vec and mark them as accessed
//         for (ll idx = 0; idx < numAmps; ++ idx) {
//             if (ampidx + strides[idx] >= sv.row) {
//                 cout << "[ERROR] Exceed the length of the state vector." << endl;
//                 exit(1);
//             }
//             amps_vec.data[idx][0] = sv.data[ampidx + strides[idx]][0];
//             isAccessed[ampidx + strides[idx]] = true;
//         }

//         // 3. Check the control bits of the current amplitude
//         //    If the control bits are not satisfied, skip this amplitude group
//         if (! isLegalControlPattern(ampidx, gate)) {
//             continue;
//         }

//         // 4. Apply the gate matrix and update the state vector in place
//         amps_vec = (*gate.gmat) * amps_vec;
//         for (ll idx = 0; idx < numAmps; ++ idx) {
//             sv.data[ampidx + strides[idx]][0] = amps_vec.data[idx][0];
//         }
//     }
// }

void applyMultiTargs(Matrix<DTYPE>& sv, QGate& gate) {
    int numTargets = gate.numTargets();
    int numAmps = (1 << numTargets); // the number of amplitudes involved in matrix-vector multiplication

    // 1. Calculate the strides for the involved amplitudes
    vector<ll> strides(numAmps, 0);
#ifdef OMP_ENABLED
#pragma omp parallel for
#endif
    for (int idx = 0; idx < numAmps; ++ idx) {
        ll stride = 0;
        for (int j = 0; j < gate.numTargets(); ++ j) {
            if (idx & (1 << j)) { // if the j-th bit of idx is 1
                stride += (1 << gate.targetQubits[j]);
            }
        }
        strides[idx] = stride;
    }

    // 2. Iterate over all amplitudes
#ifdef OMP_ENABLED
#pragma omp parallel for
#endif
    for (ll ampidx = 0; ampidx < sv.row; ++ ampidx) {
        // 2.1. Skip some of the amplitudes
        bool isStart = true;
        for (const auto& qid : gate.targetQubits) {
            if (ampidx & (1 << qid)) { // ...0...0...
                isStart = false;
                break;
            }
        }
        if (! isStart || ! isLegalControlPattern(ampidx, gate)) continue;

        // 2.2. Save the involved amplitudes to amps_vec
        Matrix<DTYPE> amps_vec(numAmps, 1); // save the involved amplitudes
        for (int idx = 0; idx < numAmps; ++ idx) {
            amps_vec.data[idx][0] = sv.data[ampidx | strides[idx]][0];
        }

        // 2.3. Apply the gate matrix and update the state vector in place
        // amps_vec = (*gate.gmat) * amps_vec;
        if (gate.isDiagonal) {
            for (int idx = 0; idx < numAmps; ++ idx) {
                sv.data[ampidx | strides[idx]][0] = gate.gmat->data[idx][idx] * amps_vec.data[idx][0];
            }
        } else {
            for (int i = 0; i < numAmps; ++ i) {
                sv.data[ampidx | strides[i]][0] = 0;
                for (int j = 0; j < numAmps; ++ j) {
                    sv.data[ampidx | strides[i]][0] += gate.gmat->data[i][j] * amps_vec.data[j][0];
                }
            }
        }
    }
}

/**
 * @brief Check if the index of an amplitude is a legal control pattern of the gate
 * 
 * @param ampidx the amplitude index
 * @param gate the processing gate
 * @return true  ampidx is a legal control pattern
 * @return false ampidx is an illegal control pattern
 */
bool isLegalControlPattern(ll ampidx, QGate& gate) {
    int ctrl;
    ll ctrlmask = 0;
    for (int i = 0; i < gate.numControls(); ++ i) {
        // Check the control qubits of the gate ////////////////
        // [HINT] If the i-th bit of amp is 0 and q_i is a |1> control qubit of gate, return false. 
        ctrl = gate.controlQubits[i];
        ctrlmask = (1 << ctrl);

        // suppose a negative ctrl means 0-controlled and a positive ctrl means 1-controlled
        // [TODO] 0-controlled and the control qubit of amp is 1
        // if (ctrl < 0 && (ampidx & ctrlmask) == 1) {
        //     return false;
        // }

        // 1-controlled and the control qubit of amp is 0
        if (ctrl >= 0 && (ampidx & ctrlmask) == 0) {
            return false;
        }
    }
    return true;
}
