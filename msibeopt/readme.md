## Compile

``
mkdir build
``

``
cd ./build
``

``
cmake ..
``

``
make ..
``

## Usage

``
msibeopt fp alpha beta theta rd=[0,1,2] sf=[0,1] pf=[0,1] sp=[0,1] sm=[0,1,2,3,4]
``

This executable accepts the following parameters:

- "fp": the path of the selected dataset

- "alpha, beta, theta": the three constraints

- "rd": "0" for executing algorithm directly, "1" for performing reduction then running algorithm, and "2" for performing reduction only

- "sf, pf, sp": flags for the pruning-aware expanding strategy, the powerful pruning rule, and the special case, respectively

- "sm": the order of vertices, "0" for original order, "1" for degree order, "2" for positive degree order, "3" for core order, and "4" for positive core order
