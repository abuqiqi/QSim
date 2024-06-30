# QSim

This repo provides some basic data structures and functions necessary for quantum circuit simulation. 

## 1. Basic Data Structures

### 1.1. The Matrix Structure

> matrix.[h/cpp]

In `matrix.[h/cpp]`, we implement the matrix structure. The gate matrix of common quantum gates are recorded. Related operations, including matrix initialization, tensor-product, addition and multiplication are implemented. 

The data type of the matrix elements can be set in `matrix.h` as follows. 

```cpp
#define DTYPE complex<double>
```

### 1.2. Quantum Gates

> qgate.[h/cpp]
//
In `qgate.[h/cpp]`, we implement the structure of quantum gates `QGate`, which has four data members, i.e., `gname`, `controlQubits`, `targetQubits`, `gmat`. To save memory footprints, `gmat` is a shared pointer to a global gate matrix map `MatrixDict` defined in `matrix.[h/cpp]`. 
Some common gate matrices are already initialized. 
When creating a new gate, we can first add a new gate matrix entry to `MatrixDict` by modifying function `void Matrix<T>::initMatrixDict()`. Then, make the `gmat` member of this new gate point to this entry. 

### 1.3. Quantum Circuits

> qcircuit.[h/cpp]

In `qcircuit.[h/cpp]`, we implement the structure of quantum circuits and provide an interface for creating a quantum circuit and adding gates. 

## 2. Quantum Circuit Simulation

### 2.1. State Vector Simulation

> svsim.[h/cpp]

State vector simulation is a useful algorithm for simulation of quantum circuits, especially in distributed situations. Its implementation is in `svsim.[h/cpp]`. 

### 2.2. Generate 2-qubit Controlled Gate Matrices

> gen_gate.[h/cpp]

In `gen_gate.cpp`, we implement the generation of the operation-matrix for arbitrary controlled gates. 

In *QuanPath* and *QuanTrans*, gates on low-order qubits are simulated using the indexing algorithm locally on each node (which will be discussed below), while operations on the high-order `h` qubits require constructing the operation-matrix. Therefore, we can relabel the high-order `h` qubits from lower qubits to higher qubits as `0`, `1`, `2`, ..., `h-1`.

**Case 1: a higher qubit controls a lower qubit**


**Case 2: a lower qubit controls a higher qubit**
