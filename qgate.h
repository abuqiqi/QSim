#pragma once

#include "matrix.h"

class QGate {
public:
    string gname;
    vector<int> controlQubits;
    vector<int> targetQubits;
    shared_ptr<Matrix<DTYPE>> gmat;
    
    QGate(string gname_, vector<int> controls_, vector<int> targets_);
    QGate(string gname_, vector<int> controls_, vector<int> targets_, double theta);
    QGate(const QGate& other);

    QGate& operator=(const QGate& other);

    int numQubits();
    int numControlQubits();
    int numTargetQubits();

    bool isIDE();
    bool isSingle();
    bool is2QubitControl(int qid);
    bool is2QubitControlled(int qid);
    bool isTargetQubit(int qid);
    bool isControlQubit(int qid);

    void print();

    ~QGate();
};