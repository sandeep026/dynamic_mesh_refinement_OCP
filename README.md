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
