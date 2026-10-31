## Dynamic mesh refinement
A mesh refinement algorithm where time is a free variable in multiple shooting with piecewise constant control.

### Constraints on time
$$
\begin{bmatrix}
t_1 & t_2 & \dots & t_{N+1}
\end{bmatrix}
$$

$$
t_0=t_i
$$

$$
t_{N+1}=t_f
$$

$$
\alpha_1 \frac{(t_f-t_i)}{N} \leq t_{k+1}-t_k \leq \alpha_2 \frac{(t_f-t_i)}{N}
$$

1. $\alpha_i$ chosen heuristically ($\alpha_1= 0.1$ and $\alpha_2= 2$)
2. $\alpha_2 =1$, it becomes a uniform grid

### Problems
1. Ally Chen problem
2. batch fed bio reactor
3. Bressnan problem
4. Jennings problem
5. Catalyst mixing problem
6. chemical reactor

### References
All the OCPs in the repo and its solution can be found in the paper listed below.

Chen, Weifeng, and Lorenz T. Biegler. “Nested Direct Transcription Optimization for Singular Optimal Control Problems.” AIChE Journal, vol. 62, no. 10, Oct. 2016, pp. 3611–27. DOI.org (Crossref), https://doi.org/10.1002/aic.15272.
