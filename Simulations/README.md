The continuous-time equations of motion of the manipulator's end-effector frame are known by

$$
E = mc^2
$$

where $\fbm{x}$ is the position vector, $\fbm{R}$ is the rotation matrix representation of the orientation, $\bm{\upsilon}$ is the linear spatial velocity, and $\bm{\omega}$ is the angular body velocity.

Let $\Delta T$ be the sampling time period. Then, in discrete time, the equations of motion of the manipulator's end-effector frame can be expressed as
\begin{align*}
\fbm{x}_{k+1}=&\fbm{x}_{k}+\Delta T\bm{\upsilon}_{k}\textrm{,}\\
\fbm{R}_{k+1}=&\fbm{R}_{k}\exp\left(\Delta T\bm{\omega}^\times_k\right)\textrm{,}
\end{align*}
at any time step~$k=0,1,2,\cdots$, where $\fbm{x}_{\ast}$ and $\fbm{R}_{\ast}$ denote the position vectors and rotation matrices at the time steps~$\ast=k$ and $\ast=k+1$, and $\bm{\upsilon}_{k}$ and $\bm{\omega}_{k}$ denote the linear spatial and the angular body velocities at the time step~$k$, respectively.
