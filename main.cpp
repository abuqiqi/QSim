#include "qcircuit.h"

/**
 * @brief Calculate the number of multiplication when simulating a quantum circuit
 * 
 * @param qc 
 * @return int 
 */
int cost(QCircuit qc) {
    int numMult = 0;
    ll N = 1 << qc.numQubits;

    for (int j = 0; j < qc.numDepths; ++ j) {
        for (int i = 0; i < qc.gates[j].size(); ++ i) {
            if (qc.gates[j][i].isSingle()) {
                numMult += 2 * N;
            } else if (qc.gates[j][i].is2QubitControlled()) {
                numMult += N;
            }
        }
    }
    cout << "[INFO] Number of multiplication: " << numMult << "\n";
    return numMult;
}

/**
 * @brief Merge some gates
 * 
 * @param qc 
 * @return QCircuit 
 */
QCircuit merge(QCircuit qc) {
    
}

int main() {
    int numQubits = 4;
    QCircuit qc = QAOA(numQubits);

    cost(qc);

    QCircuit qc2 = qc;

    return 0;
}