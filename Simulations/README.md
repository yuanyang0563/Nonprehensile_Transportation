The continuous-time equations of motion of the manipulator's end-effector frame are known by

$$
\begin{align*}
\dot{\mathbf{x}}=&\mathbf{\upsilon}\textrm{,}\\
\dot{\mathbf{R}}=&\mathbf{R}\mathbf{\omega}^\times\textrm{,}
\end{align*} 
$$

where $\fbm{x}$ is the position vector, $\mathbf{R}$ is the rotation matrix representation of the orientation, $\mathbf{\upsilon}$ is the linear spatial velocity, and $\bmathbf{\omega}$ is the angular body velocity.

Let $\Delta T$ be the sampling time period. Then, in discrete time, the equations of motion of the manipulator's end-effector frame can be expressed as

$$
\begin{align*}
\mathbf{x}_{k+1}=&\mathbf{x}_{k}+\Delta T\mathbf{\upsilon}_{k}\textrm{,}\\
\mathbf{R}_{k+1}=&\mathbf{R}_{k}\exp\left(\Delta T\mathbf{\omega}^\times_k\right)\textrm{,}
\end{align*}
$$

at any time step~$k=0,1,2,\cdots$, where $\fbm{x}_{\ast}$ and $\fbm{R}_{\ast}$ denote the position vectors and rotation matrices at the time steps~$\ast=k$ and $\ast=k+1$, and $\bm{\upsilon}_{k}$ and $\bm{\omega}_{k}$ denote the linear spatial and the angular body velocities at the time step~$k$, respectively.
