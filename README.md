%------------------------**Convex programming for GME-dimension**\------------------------%

Multipartite GME states with high entangled degrees of freedom are interesting from a fundamental quantum information science perspective and offer advantages in a range of applications. The GME-dimention quantifies the entangled degrees of freedom in a high-dimensional GME state as the lowest Schmidt number needed to generate the state by classical mixture over all bipartitions.

Based on the methods outlined in \[1\] we build on the code of \[4\] to detect the GME-dimention. This work introduces:

* An SDP to compute the critical visibility of any target state.  
* An SDP to simulate a state under projective measurements using mutually unbiased bases (MUBs).  
* Generating d-partite d-dimensional supersinglet states.

%\--------------------------------------\-**List of the functions**\---------------------------------------%

* SDPVisibility  
  * VisibilitySDP  
    Function evaluates the SemidefiniteProgram (SDP) for a target n-partite d-dimensional state using YALMIPn

  * reductionMap  
    Applies the general reduction map

  * SetPartition\[2\]

  * Stirling2nd\[2\]

* MeasurementStatisticsSDP  
  * MeasurementStatisticsSDP  
    Function evaluates the measurement statistic SDP for a target n-partite d-dimensional state using YALMIP

  * GenerateMUB  
    Generates all possible MUBs

  * GenAllK  
    This function generates all iterations of n numbers ranging from d

  * reductionMap  
    Applies the general reduction map

  * SetPartition\[2\]

  * Stirling2nd\[2\]

  * mub \[3\]

* SuperSinglet  
  * SuperSinglet  
    Creates a supersinglet state

  * LeviCivita  
    Evaluates the sign of a permutation

\[1\]: Gabriele Cobucci, Armin Tavakoli (2024). Characterising and detecting genuinely high-dimensional genuine multipartite entanglement ([https://arxiv.org/abs/2402.06234](https://arxiv.org/abs/2402.06234))  
\[2\]: Bruno Luong (2024). Set partition ([https://www.mathworks.com/matlabcentral/fileexchange/24133-set-partition](https://www.mathworks.com/matlabcentral/fileexchange/24133-set-partition)), MATLAB Central File Exchange. Retrieved May 12, 2025\.  
\[3\]: Sebastien Designolle (2024). MUB  
([https://github.com/sebastiendesignolle/MUB](https://github.com/sebastiendesignolle/MUB)),  
GitHub. Retrieved May 12, 2025\.
\[4\]: Gabriele Cobucci (2024). GMEN-mixer ([https://github.com/GabrieleCobucci/GMEN-mixer/tree/GME-dimension](https://github.com/GabrieleCobucci/GMEN-mixer/tree/GME-dimension))  
