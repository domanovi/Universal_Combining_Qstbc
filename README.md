# Universal_Combining_Qstbc

Steps:
===================================================================================================
  1. open main.py
  2. change the values of 'n_t' and 'n_r' to be the number of the transmit and receive antennas.
  3. run the program, the results are printed in the screen.
===================================================================================================

Results representation:
===================================================================================================
Symbolic representation includes:
  Target: EA-QSTBC and the "desired stream" (stream of symbols used for detection) in the MISO channel.
  Result: the proposed scheme for the MIMO channel. includes: transmission matrix, combining scheme, and channel coefficients transformation.
  Check Results: Comparison between the signal part of the stream used for detection in both the MISO and the MIMO channels.
  
  * Note 1: all the presented matrices are symbolic (i.e. with indexed letters, such that the indices begin from 0)
  * Note 2: the time index here is denoted with paranthesis.
  
 Constant matrices representation includes:
  Transmission: N_t matrices representing the transmission scheme.
  Combining: N_r matrices representing the combining scheme.
  
  * Note 1:
    for the transmission matrices, each matrix is multiplied from right by the information symbols vector: [[ x₀   x₁   x₀⃰  x₁⃰]]ᵀ.
  * Note 2:
    The stream received in rx i (and its conjugates terms):
    [[ rⱼ(0)   rⱼ(1)   rⱼ(0)⃰   rⱼ(1)⃰ ]]ᵀ
    The stream received in rx j (and its conjugates terms) is multiplied
    from left by the corresponding combining matrix and then the results are summed up.

===================================================================================================

Example:
===================================================================================================

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 1x0:: Combining Scheme~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 ================================= Target ==================================
Nt = 1
Nr = 2
Target QSTBC:
[[ x₀   x₁ ]
 [-x₁⃰  x₀⃰]]
Desired Stream:
ŝ = [[ x₀   x₁ ] [[ ĥ₀ ]   =  [[+   x₀  •  ĥ₀    +   x₁  •  ĥ₁    ] 
     [-x₁⃰  x₀⃰]] [ ĥ₁ ]]      [-   x₁⃰ •  ĥ₀    +   x₀⃰ •  ĥ₁    ]]

 ================================= Result ==================================
Transmission Matrix:
[[ x₀ ]
 [ x₁⃰]]
Combining Scheme:
[[+  r₀(0)   +  r₁⃰(1)   ]
 [-  r₀(1)   +  r₁⃰(0)   ]]
Channel coefficients transformation:
[[ ĥ₀ ]  = [[ h₀₋₀ ] 
 [ ĥ₁ ]]    [ h₀₋₁⃰]]

 =============================== Check Result ==============================
Received Signals Matrix:
                               0                    1
R    =    0  +   x₀  •  h₀₋₀      +   x₀  •  h₀₋₁    
          1  +   x₁⃰ •  h₀₋₀      +   x₁⃰ •  h₀₋₁    
Received Stream (after combining):
s = [[+   x₀  •  h₀₋₀    +   x₁  •  h₀₋₁⃰   ] 
     [-   x₁⃰ •  h₀₋₀    +   x₀⃰ •  h₀₋₁⃰   ]]
Compare 'Desired Stream' and 'Received Stream':
[[ True]
 [ True]]

 ========================== Comparison Succeeded !!! ==========================
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 1x0:: Constant Matrices~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
Nt = 1
Nr = 2

 =============================== Transmission Matrices ==============================
Each matrix is multiplied from right by the information symbols vector:
[[ x₀   x₁   x₀⃰  x₁⃰]]ᵀ
Tx 0:              
[[1. 0. 0. 0.]     
 [0. 0. 0. 1.]]    

 =============================== Combining Matrices ==============================
The stream received in rx i (and its conjugates terms):
[[ rⱼ(0)   rⱼ(1)   rⱼ(0)⃰   rⱼ(1)⃰ ]]ᵀ
The stream received in rx i (and its conjugates terms) is multiplied
 from left by the corresponding combining matrix and then the results are summed up.
Rx 0:                  Rx 1:              
[[ 1.  0.  0.  0.]     [[0. 0. 0. 1.]     
 [ 0. -1.  0.  0.]]     [0. 0. 1. 0.]]    

