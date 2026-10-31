# Dynamic mesh refinement
A mesh refinement algorithm where  time is a free variable.

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
