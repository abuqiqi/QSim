#include "qgate.h"

// string to_string(double d) {
//     ostringstream os;
//     os << fixed << setprecision(6) << d;
//     return os.str();
// }

QGate::QGate(string gname_, vector<int> controls_, vector<int> targets_) {
    gname = gname_;
    controlQubits = controls_;
    targetQubits = targets_;

    gmat = Matrix<DTYPE>::MatrixDict[gname];
    if (gmat == nullptr) {
        cout << "[ERROR] Gate " << gname << " not found in MatrixDict" << endl;
        exit(1);
    }
}

QGate::QGate(string gname_, vector<int> controls_, vector<int> targets_, double theta) {
    gname = gname_;
    controlQubits = controls_;
    targetQubits = targets_;

    string matkey = gname + to_string(theta);
    gmat = Matrix<DTYPE>::MatrixDict[matkey];

    if (gmat != nullptr) { // the gate matrix already exists
        cout << "[DEBUG] Matrix already exists: " << matkey << ", " << gmat << endl;
        return;
    }

    Matrix<DTYPE> mat;

    if (gname == "RX") {
        mat.rotationX(theta);
    } else if (gname == "RY") {
        mat.rotationY(theta);
    } else if (gname == "RZ") {
        mat.rotationZ(theta);
    } else {
        cout << "[ERROR] Gate " << gname << " not implemented" << endl;
        exit(1);
    }
    Matrix<DTYPE>::MatrixDict[matkey] = make_shared<Matrix<DTYPE>>(move(mat));
    gmat = Matrix<DTYPE>::MatrixDict[matkey];
}

QGate::QGate(const QGate& other) {
    gname = other.gname;
    controlQubits = other.controlQubits;
    targetQubits = other.targetQubits;
    gmat = other.gmat;
}

QGate& QGate::operator=(const QGate& other) {
    gname = other.gname;
    controlQubits = other.controlQubits;
    targetQubits = other.targetQubits;
    gmat = other.gmat;
    return *this;
}

int QGate::numQubits() {
    return controlQubits.size() + targetQubits.size();
}

int QGate::numControlQubits() {
    return controlQubits.size();
}

int QGate::numTargetQubits() {
    return targetQubits.size();
}

bool QGate::isIDE() {
    return gname == "IDE";
}

bool QGate::isTargetQubit(int qid) {
    return find(targetQubits.begin(), targetQubits.end(), qid) != targetQubits.end();
}

bool QGate::isControlQubit(int qid) {
    return find(controlQubits.begin(), controlQubits.end(), qid) != controlQubits.end();
}

void QGate::print() {
    cout << "===== Gate: " << gname << " =====" << endl;
    cout << "Control qubits: ";
    for (vector<int>::size_type i = 0; i < controlQubits.size(); i++) {
        cout << controlQubits[i] << " ";
    }
    cout << endl;
    cout << "Target qubits: ";
    for (vector<int>::size_type i = 0; i < targetQubits.size(); i++) {
        cout << targetQubits[i] << " ";
    }
    cout << endl;
    gmat->print();
}

QGate::~QGate() {
    return;
}