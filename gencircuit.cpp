#include "gencircuit.h"

vector<QGate> testkqubits(int k) {
    QCircuit subqc(k, "subqc");
    vector<int> targs;
    for (int i = 0; i < k; ++ i) {
        subqc.ry((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
        targs.push_back(i);
    }
    // subqc.print();
    Matrix<DTYPE> subU = getOperationMatrix(subqc);
    QGate kqubitgate(subqc.cmatkey(), {}, targs, subU);
    // kqubitgate.print();

    vector<QGate> gateSeq;
    // for (int i = 0; i < 120 / k; ++ i) {
    for (int i = 0; i < 100; ++ i) {
        gateSeq.push_back(kqubitgate);
    }
    return gateSeq;
}

vector<QGate> testkfusion(int k) {
    QCircuit subqc(k, "subqc");
    vector<int> targs;
    for (int i = 0; i < k; ++ i) {
        subqc.ry((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
        targs.push_back(i);
    }
    // subqc.print();
    Matrix<DTYPE> subU = getOperationMatrix(subqc);
    QGate kqubitgate(subqc.cmatkey(), {}, targs, subU);
    // kqubitgate.print();

    vector<QGate> gateSeq;
    for (int i = 0; i < 120 / k; ++ i) {
        gateSeq.push_back(kqubitgate);
    }
    return gateSeq;
}

// a controlled gate with k control qubits
vector<QGate> testkctrls(int k) {
    // generate a controlled gate with k control qubits
    vector<int> ctrls(k);
    for (int i = 0; i < k; ++ i) {
        ctrls[i] = i;
    }

    vector<QGate> gateSeq;
    QGate kctrlgate("H", ctrls, {k + 1});
    // kqubitgate.print();

    for (int i = 0; i < 100; ++ i) {
        gateSeq.push_back(kctrlgate);
    }
    return gateSeq;
}

QCircuit test(int numQubits) {
    // test circuit
    QCircuit qc(numQubits, "test");
    qc.h(0);
    qc.x(0);
    qc.cx(1, 0);
    qc.h(0);
    qc.x(0);
    qc.h(1);
    qc.x(1);
    qc.print();
    return qc;
}

QCircuit QFT(int numQubits) {
    QCircuit qc(numQubits, "QFT");
    for (int i = numQubits - 1; i >= numQubits / 2; -- i) {
        qc.h(numQubits-1-i);
    }
    for (int i = numQubits - 1; i >= 0; -- i) {
        qc.h(numQubits-1-i);
        for (int j = i - 1; j >= 0; -- j) {
            int k = i - j + 1;
            qc.cu1(2 * acos(-1.0) / pow(2, k), numQubits-1-j, numQubits-1-i);
        }
    }

    // qc.printInfo();
    qc.print();
    return qc;
}

QCircuit Adder(int numQubits) {
    QCircuit qc(numQubits, "Adder");
    qc.h(0);
    qc.h(1);
    // QFT
    for (int i = numQubits/2; i < numQubits; ++ i) {
        qc.h(i);
        for (int j = i+1; j < numQubits; ++ j) {
            int k = j - i + 1;
            qc.cu1(2 * acos(-1.0) / pow(2, k), j, i);
        }
    }
    // addition
    for (int i = 0; i < numQubits/2; ++ i) {
        int targ = i + numQubits/2;
        int k = 0;
        for (int j = i; j < numQubits/2; ++ j) {
            ++ k;
            qc.cu1(2 * acos(-1.0) / pow(2, k), j, targ);
        }
    }
    // QFT^(-1)
    qc.h(numQubits-1);
    for (int i = numQubits-2; i > numQubits/2-1; -- i) {
        for (int j = numQubits-1; j > i; -- j) {
            int k = j - i + 1;
            qc.cu1(-2 * acos(-1.0) / pow(2, k), j, i);
        }
        qc.h(i);
    }
    // qc.print();
    return qc;
}

QCircuit QFT_Quirk(int numQubits) {
    QCircuit qc(numQubits, "QFT_Quirk");

    qc.h(numQubits-1);
    qc.h(numQubits-2);

    int ii = 0;
    int jj = numQubits - 1;
    while (ii < jj) {
        qc.swap(ii, jj);
        ii ++;
        jj --;
    }

    qc.h(0);
    qc.cs(1, 0);
    qc.barrier();

    qc.h(1);
    qc.cx(2, 0);
    qc.cx(2, 1);
    qc.barrier();

    int n = numQubits - 1;
    int i = 2; // H gate
    while (i < n) {
        qc.h(i);
        int t = i;
        int tt = 0;

        while (t >= 2) {
            qc.cx(i+1, tt);
            tt += 1;
            t -= 1;
        }
        qc.cx(i+1, i-1);
        qc.cx(i+1, i);
        i += 1;
        qc.barrier();
    }
    qc.h(n);

    // qc.print();
    return qc;
}

/**
 * @brief Generate a QAOA quantum circuit
 * 
 * @param numQubits #Qubits
 * @return QCircuit  the generated quantum circuit
 */
QCircuit QAOA(int numQubits) {
    QCircuit qc(numQubits, "QAOA");

    for (int p = 0; p < 1; ++ p) {
        for (int i = 0; i < numQubits; ++ i)
            qc.h(i);
        
        for (int j = 0; j < numQubits; ++ j) {
            for (int i = j + 1; i < numQubits; ++ i) {
                qc.cx(j, i);
                qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
                qc.cx(j, i);
            }
        }

        for (int i = 0; i < numQubits; ++ i)
            qc.h(i);
        for (int i = 0; i < numQubits; ++ i)
            qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
        for (int i = 0; i < numQubits; ++ i)
            qc.h(i);
    }
    qc.printInfo();
    return qc;
}

// IQP
QCircuit IQP(int numQubits) {
    QCircuit qc(numQubits, "IQP");
    for (int i = 0; i < numQubits; ++ i) {
        qc.h(i);
    }

    for (int i = 0; i < numQubits; ++ i) {
        for (int j = i + 1; j < numQubits; ++ j) {
            qc.cu1(random(100) * acos(-1.0) / 2, j, i);
        }
    }

    for (int i = 0; i < numQubits; ++ i) {
        qc.u1(random(100) * acos(-1.0) / 8, i);
    }

    for (int i = 0; i < numQubits; ++ i) {
        qc.h(i);
    }

    // qc.print();
    return qc;
}

/**
 * @brief Generate a Grover's circuit
 * 
 * @param numQubits #Qubits
 * 
 * @return a quantum circuit
 */
QCircuit Grover(int numQubits, int k) {
    int numIterations = 5;
    int totalIterations = floor(std::acos(-1.0) / 4 * sqrt(pow(2, numQubits) / pow(2, k)));

    QCircuit qc(numQubits, "Grover_"+to_string(numIterations)+"-"+to_string(totalIterations));

    for (int i = 0; i < numQubits; ++ i) {
        qc.h(i);
    }

    vector<int> oracle, controlQubits;
    for (int i = 0; i < numQubits-k-1; ++ i) {
        oracle.push_back(i);
    }
    for (int i = 0; i < numQubits-1; ++ i) {
        controlQubits.push_back(i);
    }

    for (int itr = 0; itr < numIterations; ++ itr) {
        qc.ccz(oracle, numQubits-k-1); // oracle
        qc.barrier();
        for (int i = 0; i < numQubits; ++ i) {
            qc.h(i);
        }
        for (int i = 0; i < numQubits; ++ i) {
            qc.x(i);
        }
        qc.ccz(controlQubits, numQubits-1);
        qc.barrier();
        for (int i = 0; i < numQubits; ++ i) {
            qc.x(i);
        }
        for (int i = 0; i < numQubits; ++ i) {
            qc.h(i);
        }
    }

    qc.print();
    // qc.printInfo();
    return qc;
}

QCircuit VQC_NN(int numQubits) { // nearest neighbour
    QCircuit qc = QCircuit(numQubits, "VQC_NN");

#if 1
    for (int j = 0; j < 1; ++ j) {
        for (int i = 0; i < numQubits; ++ i)
            qc.rx((double)rand() / RAND_MAX * 2 * acos(-1.0), i);

        for (int i = 0; i < numQubits; ++ i)
            qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);

        // for (int i = 0; i + 1 < numQubits; i += 2)
        //     qc.crx((double)rand() / RAND_MAX * 2 * acos(-1.0), i, i+1);

        for (int i = numQubits - 2; i > 0; i -= 2)
            qc.crx((double)rand() / RAND_MAX * 2 * acos(-1.0), i, i-1);
            // qc.cx(i, i-1);

        for (int i = numQubits - 1; i > 0; i -= 2)
            qc.crx((double)rand() / RAND_MAX * 2 * acos(-1.0), i-1, i);
            // qc.cx(i, i-1);

        for (int i = 0; i < numQubits; ++ i)
            qc.rx((double)rand() / RAND_MAX * 2 * acos(-1.0), i);

        for (int i = 0; i < numQubits; ++ i)
            qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);

        // for (int i = 1; i + 1 < numQubits; i += 2)
        //     qc.crx((double)rand() / RAND_MAX * 2 * acos(-1.0), i, i+1);


        // for (int i = 1; i + 1 < numQubits; i += 2) {
        //     qc.crz((double)rand() / RAND_MAX * 2 * acos(-1.0), i, i+1);
        // }

        // for (int i = numQubits - 1; i > 0; i -= 2) {
        //     qc.crz((double)rand() / RAND_MAX * 2 * acos(-1.0), i, i-1);
        // }
        // for (int i = numQubits - 2; i > 0; i -= 2) {
        //     qc.crz((double)rand() / RAND_MAX * 2 * acos(-1.0), i, i-1);
        // }

        // // levels of CU
        // for (int i = 0; i < numQubits-1; ++ i) {
        //     // qc.cx(i, i+1);
        //     qc.crz((double)rand() / RAND_MAX * 2 * acos(-1.0), i, i+1);
        //     qc.barrier();
        // }

        // for (int i = 1; i < numQubits-1; ++ i) {
        //     // qc.cx(i, i+1);
        //     qc.crz((double)rand() / RAND_MAX * 2 * acos(-1.0), i, i+1);
        //     qc.barrier();
        // }

        // for (int i = 0; i < numQubits; ++ i)
        //     qc.ry((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
        // for (int i = 0; i < numQubits; ++ i)
        //     qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
    }
#endif
    qc.print();
    return qc;
}

// VQC circuit-block
QCircuit VQC_CB(int numQubits) {
    QCircuit qc = QCircuit(numQubits, "VQC_CB");

    for (int j = 0; j < 3; ++ j) {
        for (int i = 0; i < numQubits; ++ i)
            qc.rx((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
        for (int i = 0; i < numQubits; ++ i)
            qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);

        for (int i = 0; i < numQubits-1; ++ i) {
            // qc.cx(i, i+1);
            qc.crz((double)rand() / RAND_MAX * 2 * acos(-1.0), i, i+1);
            qc.barrier();
        }
        qc.crz((double)rand() / RAND_MAX * 2 * acos(-1.0), numQubits-1, 0);
    }

    // qc.print();
    return qc;
}

QCircuit VQC_AA(int numQubits) {
    QCircuit qc = QCircuit(numQubits, "VQC_AA");
    for (int _ = 0; _ < 1; ++ _) {
        for (int i = 0; i < numQubits; ++ i)
            qc.rx((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
        for (int i = 0; i < numQubits; ++ i)
            qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);

        for (int i = 0; i < numQubits; ++ i) {
            for (int j = 0; j < numQubits; ++ j) {
                if (i == j)
                    continue;
                qc.barrier();
                qc.crz((double)rand() / RAND_MAX * 2 * acos(-1.0), i, j);
            }
        }
    }
    return qc;
}

/**
 * @brief Generate VQC
 * 
 * @param numQubits #Qubits
 * 
 * @return a quantum circuit
 */
QCircuit VQC(int numQubits) {
    QCircuit qc = QCircuit(numQubits, "VQC");

    for (int j = 0; j < 1; ++ j) {
        // 2 levels of RY
        for (int k = 0; k < 2; ++ k)
            for (int i = 0; i < numQubits; ++ i)
                qc.ry((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
                // qc.ry(i, i);
        // levels of CX
        for (int i = 0; i < numQubits-1; ++ i) {
            // qc.cx(i, i+1);
            qc.crz((double)rand() / RAND_MAX * 2 * acos(-1.0), i, i+1);
            qc.barrier();
        }
        // 2 levels of RY
        for (int k = 0; k < 2; ++ k)
            for (int i = 0; i < numQubits; ++ i)
                qc.ry((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
                // qc.ry(i, i);
    }

    qc.printInfo();
    return qc;
}

/**
 * @brief Generate VQC1 (a small number of inseparable levels)
 * 
 * @param numQubits #Qubits
 * 
 * @return a quantum circuit
 */
QCircuit VQC1(int numQubits) {
    QCircuit qc = QCircuit(numQubits, "VQC1");

    for (int j = 0; j < 1; ++ j) {
        for (int k = 0; k < 2; ++ k)
            for (int i = 0; i < numQubits; ++ i)
                qc.ry((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
        for (int i = 0; i < numQubits-1; ++ i)
            qc.cx(i, i+1);
        for (int k = 0; k < 2; ++ k)
            for (int i = 0; i < numQubits; ++ i)
                qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
    }

    for (int j = 0; j < 1; ++ j) {
        for (int k = 0; k < 2; ++ k)
            for (int i = 0; i < numQubits; ++ i)
                qc.ry((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
        for (int i = 1; i < numQubits; ++ i)
            qc.cx(i, i-1);
        for (int k = 0; k < 2; ++ k)
            for (int i = 0; i < numQubits; ++ i)
                qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
    }

    for (int j = 0; j < 1; ++ j) {
        for (int k = 0; k < 2; ++ k)
            for (int i = 0; i < numQubits; ++ i)
                qc.ry((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
        for (int i = 0; i < numQubits-1; ++ i)
            qc.cx(i, i+1);
        for (int k = 0; k < 2; ++ k)
            for (int i = 0; i < numQubits; ++ i)
                qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
    }

    // qc.print();
    return qc;
}

/**
 * @brief Generate VQC2 (circuit 5, a larger number of inseparable levels)
 * 
 * @param numQubits #Qubits
 * 
 * @return a quantum circuit
 */
QCircuit VQC2(int numQubits) {
    QCircuit qc = QCircuit(numQubits, "VQC2");

    // for (int _ = 0; _ < 2; ++ _) {
    for (int i = 0; i < numQubits; ++ i)
        qc.rx((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
    for (int i = 0; i < numQubits; ++ i)
        qc.rx((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
    for (int i = 0; i < numQubits; ++ i)
        qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
    for (int i = 0; i < numQubits; ++ i)
        qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
    for (int i = 0; i < numQubits; ++ i) {
        for (int j = 0; j < numQubits; ++ j) {
            if (i == j)
                continue;
            qc.barrier();
            qc.crz((double)rand() / RAND_MAX * 2 * acos(-1.0), i, j);
        }
    }
    for (int i = 0; i < numQubits; ++ i)
        qc.rx((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
    for (int i = 0; i < numQubits; ++ i)
        qc.rx((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
    for (int i = 0; i < numQubits; ++ i)
        qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
    for (int i = 0; i < numQubits; ++ i)
        qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
    // }
    // for (int k = 0; k < 2; ++ k)
    //     for (int i = 0; i < numQubits; ++ i)
    //         qc.rx((double)rand() / RAND_MAX * 2 * acos(-1.0), i);

    // for (int i = 0; i < numQubits; ++ i) {
    //     for (int k = 0; k < 1; ++ k)
    //         for (int ii = 0; ii < numQubits; ++ ii)
    //             qc.ry((double)rand() / RAND_MAX * 2 * acos(-1.0), ii);

    //     for (int j = 0; j < numQubits; ++ j) {
    //         if (i == j)
    //             continue;
    //         qc.barrier();
    //         qc.cx(i, j);
    //     }
    // }

    // for (int k = 0; k < 2; ++ k)
    //     for (int i = 0; i < numQubits; ++ i)
    //         qc.ry((double)rand() / RAND_MAX * 2 * acos(-1.0), i);

    // qc.print();
    return qc;
}

QCircuit RandomRegular(int numQubits, int numDepths) {
    QCircuit qc = QCircuit(numQubits, "RandomRegular");
    int gTyp;

    while (true) {
        // add two levels of CX gates
        if (qc.numDepths % 10 == 2) {
            for (int j = numQubits - 1; j > 0; j -= 2) {
                qc.cx(j, j-1);
            }
            for (int j = numQubits - 2; j > 0; j -= 2) {
                qc.cx(j, j-1);
            }
        }
        // add one level of single-qubit gates
        else {
            // random single-qubit gates
            for (int j = 0; j < numQubits; ++ j) {
                gTyp = random(3);
                if (gTyp == 0) {
                    qc.h(j);
                }
                else if (gTyp == 1) {
                    qc.z(j);
                }
                else {
                    qc.x(j);
                }
            }
        }

        if (qc.numDepths >= numDepths) {
            break;
        }
    }

    // qc.print();
    return qc;
}

QCircuit RandomMedium(int numQubits, int numDepths) {
    QCircuit qc = QCircuit(numQubits, "RandomMedium");
    int gTyp;

    while (true) {
        // add two levels of CX gates
        if (qc.numDepths % 10 == 2) {
            for (int j = numQubits - 1; j > 0; j -= 2) {
                qc.cx(j, j-1);
            }
            for (int j = numQubits - 2; j > 0; j -= 2) {
                qc.cx(j, j-1);
            }
        }
        // add one level of single-qubit gates
        else {
            // random single-qubit gates
            for (int j = 0; j < numQubits; ++ j) {
                gTyp = random(4);
                if (gTyp == 0) {
                    qc.h(j);
                }
                else if (gTyp == 1) {
                    qc.z(j);
                }
                else if (gTyp == 2) {
                    qc.x(j);
                } else {
                    qc.ry((double)rand() / RAND_MAX * 2 * acos(-1.0), j);
                }
            }
        }

        if (qc.numDepths >= numDepths) {
            break;
        }
    }

    qc.printInfo();
    return qc;
}

QCircuit RandomRX(int numQubits, int numDepths) {
    QCircuit qc = QCircuit(numQubits, "RandomRX");
    int gTyp;

    while (true) {
        // add two levels of CX gates
        if (qc.numDepths % 10 == 2) {
            for (int j = numQubits - 1; j > 0; j -= 2) {
                qc.cx(j, j-1);
            }
            for (int j = numQubits - 2; j > 0; j -= 2) {
                qc.cx(j, j-1);
            }
        }
        // add one level of single-qubit gates
        else {
            // random single-qubit gates
            for (int j = 0; j < numQubits; ++ j) {
                gTyp = random(4);
                if (gTyp == 0) {
                    qc.h(j);
                }
                else if (gTyp == 1) {
                    qc.z(j);
                }
                else if (gTyp == 2) {
                    qc.x(j);
                } else {
                    qc.rx((double)rand() / RAND_MAX * 2 * acos(-1.0), j);
                }
            }
        }

        if (qc.numDepths >= numDepths) {
            break;
        }
    }

    qc.printInfo();
    return qc;
}

QCircuit RandomRZ(int numQubits, int numDepths) {
    QCircuit qc = QCircuit(numQubits, "RandomRZ");
    int gTyp;

    while (true) {
        // add two levels of CX gates
        if (qc.numDepths % 10 == 2) {
            for (int j = numQubits - 1; j > 0; j -= 2) {
                qc.cx(j, j-1);
            }
            for (int j = numQubits - 2; j > 0; j -= 2) {
                qc.cx(j, j-1);
            }
        }
        // add one level of single-qubit gates
        else {
            // random single-qubit gates
            for (int j = 0; j < numQubits; ++ j) {
                gTyp = random(4);
                if (gTyp == 0) {
                    qc.h(j);
                }
                else if (gTyp == 1) {
                    qc.z(j);
                }
                else if (gTyp == 2) {
                    qc.x(j);
                } else {
                    qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), j);
                }
            }
        }

        if (qc.numDepths >= numDepths) {
            break;
        }
    }

    qc.printInfo();
    return qc;
}

QCircuit RandomRandom(int numQubits, int numDepths) {
    QCircuit qc = QCircuit(numQubits, "RandomRandom");

    int numGates = 0;

    while (true) {
        // add one level of CX gates
        // if (qc.numDepths % 10 == 2) {
        //     for (int j = numQubits - 1; j > 0; j -= 2) {
        //         qc.cx(j, j-1);
        //         ++ numGates;
        //     }
            // for (int j = numQubits - 2; j > 0; j -= 2) {
            //     qc.cx(j, j-1);
            // }
        // }
        // add one level of single-qubit gates
        // else {
            // random single-qubit gates
            for (int j = 0; j < numQubits; ++ j) {
                if (j % 2 == 0)
                    qc.rx((double)rand() / RAND_MAX * 2 * acos(-1.0), j);
                else
                    // qc.rx((double)j / numQubits * 2 * acos(-0.1), j);
                    qc.ry((double)rand() / RAND_MAX * 2 * acos(-1.0), j);
                ++ numGates;
            }
        // }

        if (qc.numDepths >= numDepths) {
            break;
        }
    }

    // qc.print();
    return qc;
}

QCircuit Hyperbolic(int numQubits) {
    QCircuit qc = QCircuit(numQubits, "Hyperbolic");

    for (int _ = 0; _ < 3; ++ _) {
        for (int i = 0; i < numQubits; ++ i) {
            qc.ry((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
        }

        for (int i = 0; i < numQubits; ++ i) {
            qc.rz((double)rand() / RAND_MAX * 2 * acos(-1.0), i);
        }

        for (int k = 0; k < numQubits; ++ k) {
            for (int i = 1; i + (1 << k) < numQubits; i += 2*(1 << k)) {
                // qc.cx(i + (1 << k), i);
                qc.crz((double)rand() / RAND_MAX * 2 * acos(-1.0), i + (1 << k), i);
            }
        }
    }

    qc.print();
    return qc;
}

QCircuit QPE(int numQubits) {
    QCircuit qc = QCircuit(numQubits, "QPE");

    qc.x(0);
    for (int i = 1; i < numQubits; ++ i) {
        qc.h(i);
    }
    qc.swap(0, numQubits-1);
    qc.swap(0, numQubits-1);

    for (int i = numQubits-1; i > 0; -- i) {
        qc.cu1(acos(-1.0) / pow(2, numQubits-i), i, 0);
    }

    // IQFT
    qc.h(numQubits-1);
    qc.swap(0, numQubits-1);
    qc.swap(0, numQubits-1);
    for (int i = numQubits-2; i > 0; -- i) {
        for (int j = numQubits-1; j > i; -- j) {
            int k = j - i + 1;
            qc.cu1(-2 * acos(-1.0) / pow(2, k), j, i);
        }
        qc.h(i);
        qc.swap(0, numQubits-1);
        qc.swap(0, numQubits-1);
    }

    qc.printInfo();
    qc.print();
    return qc;
}
