## Constants in the paper [^1]
### Constants


| Code variable | Variable | Value   | Units                          | Description               |
| ------------- | -------- | ------- | ------------------------------ | ------------------------- |
| F             | $F$      | $96485$ | $\mathrm{C \  mol^{-1}}$       | Faraday constant          |
| R             | $R$      | $8.314$ | $\mathrm{J\ mol^{-1}\ K^{-1}}$ | Universal gas constant    |
| T_ref         | $T_0$    | $298$   | $\mathrm{K}$                   | Reference temperature     |
| r0            | $r_0$    | 1e-6    | $\mathrm{m}$                   | Reference particle radius |

### Reference scales
$i=n,e,p,sep$, where
- $e$ - electrolyte phase
- $n$ - negative electrode
- $p$ - positive electrode
- $sep$ - separator

| Code variable | Variable  | Formula                                  | Units                  | Description                                                                 |
| ------------- | --------- | ---------------------------------------- | ---------------------- | --------------------------------------------------------------------------- |
| a0            | $a_0$     | $a_0=\dfrac{1}{r_0}$                     | $\mathrm{m^{-1}}$      | specific surface area scale                                                 |
| L             | $L_0$     | $L_0=L_n+L_{sep}+L_p$                    | $\mathrm{m}$           | total cell thickness scale                                                  |
| phi_scale     | $\phi_0$  | $\phi_0=\dfrac{R \  T_0}{F}$             | $\mathrm{V}$           | thermal voltage scale                                                       |
| j_scale       | $j_0$     | $j_0=\dfrac{I_0}{a_0 A L_0}$             | $\mathrm{A\ m^{-2}}$   | volumetric current density scale                                            |
| tn<br>tp      | $t_{i0}$  | $t_{i0}=\dfrac{F A L_0 c_{i,\max}}{I_0}$ | $\mathrm{s}$           | characteristic solid $(i=n,p)$ diffusion time scale for electrode particles |
| te            | $t_{e0}$  | $t_{e0}=\dfrac{F A L_0 c_{e0}}{I_0}$     | $\mathrm{s}$           | characteristic electrolyte diffusion time scale                             |
| sig_scale     | -         | $\dfrac{I_0 L_0}{\phi_0 A}$              | $\mathrm{S \  m^{-1}}$ | electrode conductivity scale                                                |
| kappa_scale   | -         | $\dfrac{I_0 L_0}{\phi_0 A}$              | $\mathrm{S \  m^{-1}}$ | electrolyte conductivity scale                                              |


### Parameters

| Variable             | Units                             | Description                                  |
| -------------------- | --------------------------------- | -------------------------------------------- |
| $\kappa$             | $\mathrm{S \  m^{-1}}$            | conductivity of electrolyte                  |
| $\phi_i$             | $\mathrm{V}$                      | potential of phase $i$                       |
| $\sigma_i$           | $\mathrm{S \  m^{-1}}$            | conductivity of phase $i$                    |
| $A$                  | $\mathrm{m^2}$                    | cell surface area                            |
| $a_i$                | $\mathrm{m^{-1}}$                 | specific surface area of phase $i$           |
| $c_i$                | $\mathrm{mol\  m^{-3}}$           | concentration of lithium in phase $i$        |
| $D_i$                | $\mathrm{m^2 \  s^{-1}}$          | lithium diffusion coefficient of phase $i$   |
| $D_i^{\mathrm{eff}}$ | $\mathrm{m^2 \  s^{-1}}$          | effective diffusion coefficient of phase $i$ |
| $I$                  | $\mathrm{A}$                      | electronic current                           |
| $j$                  | $\mathrm{A\ m^{-3}}$              | reaction current density per volume          |
| $L$                  | $\mathrm{m}$                      | cell thickness                               |
| $L_i$                | $\mathrm{m}$                      | thickness of domain $i$                      |
| $m_i$                | $\mathrm{A\ m^{2.5}\ mol^{-1.5}}$ | reaction rate in phase $i$                   |
| $R_i$                | $\mathrm{\mu m}$                  | particle radius of phase $i$                 |
| $T$                  | $\mathrm{K}$                      | absolute temperature                         |
| $t_+^{0}$ or $t^+$   | —                                 | transference number of lithium-ion           |
| $U_i$                | $\mathrm{V}$                      | open circuit potential of phase $i$          |



### Scaled (i.e. adimensionalised) parameters


| Variable                          | Formula                                                                                             | Description                                         |
| --------------------------------- | --------------------------------------------------------------------------------------------------- | --------------------------------------------------- |
| $\bar{c}_e$                       | $`\bar{c}_e = \dfrac{c_e}{c_{e0}}`$                                                                 | scaled electrolyte concentration                    |
| $\bar{x}$                         | $\bar{x} = \dfrac{x}{L_0}$                                                                          | scaled spatial coordinate                           |
| $\bar{a}_i$                       | $\bar{a}_i = \dfrac{a_i}{a_0}$                                                                      | scaled specific surface area of phase $i$           |
| $\bar{D}_e^{\mathrm{eff}}$        | $`\bar{D}_e^{\mathrm{eff}} = D_e^{\mathrm{eff}}\dfrac{t_{e0}}{L_0^{2}}`$                            | scaled effective electrolyte diffusivity            |
| $\bar{t}_e$                       | $`\bar{t}_e = \dfrac{t}{t_{e0}}`$                                                                   | scaled time                                         |
| $\bar{r}$                         | $\bar{r} = \dfrac{r}{r_0}$                                                                          | scaled radial coordinate                            |
| $\bar{R}_i$                       | $\bar{R}_i = \dfrac{R_i}{r_0}$                                                                      | scaled particle radius of phase $i$                 |
| $\bar{c}_i$                       | $`\bar{c}_i = \dfrac{c_i}{c_{i,\max}}`$                                                             | scaled lithium concentration in phase $i$           |
| $\bar{D}_i^{\mathrm{eff}}$        | $`\bar{D}_i^{\mathrm{eff}} = D_i^{\mathrm{eff}}\dfrac{t_{i0}}{r_0^{2}}`$                            | scaled effective diffusion coefficient in phase $i$ |
| $\bar{t}_i$                       | $`\bar{t}_i = \dfrac{t}{t_{i0}}`$                                                                   | scaled time for phase $i$ diffusion                 |
| $\bar{T}$                         | $\bar{T}=\dfrac{T}{T_0}$                                                                            | scaled temperature                                  |
| $\bar{\phi}_e$                    | $`\bar{\phi}_e=\dfrac{\phi_e}{\phi_0}`$                                                             | scaled electrolyte potential                        |
| $\bar{\phi}_i$                    | $`\bar{\phi}_i=\dfrac{\phi_i}{\phi_0}`$                                                             | scaled potential of phase $i$                       |
| $\bar{\kappa}^{\mathrm{eff}}$     | $`\bar{\kappa}^{\mathrm{eff}}=\kappa^{\mathrm{eff}}\left[\dfrac{\phi_0 A}{I_0 L_0}\right]`$         | scaled effective electrolyte conductivity           |
| $\bar{\kappa}^{\mathrm{eff}}_{D}$ | $`\bar{\kappa}^{\mathrm{eff}}_{D}=2\ \bar{T}\ \bar{\kappa}^{\mathrm{eff}}(1-t_+)`$                  | scaled “diffusional” conductivity factor            |
| $\bar{\sigma}^{\mathrm{eff}}_{i}$ | $`\bar{\sigma}^{\mathrm{eff}}_{i}=\sigma^{\mathrm{eff}}_{i}\left[\dfrac{\phi_0 A}{I_0 L_0}\right]`$ | scaled effective conductivity of phase $i$          |
| $\bar{j}$                         | $\bar{j}=\dfrac{j}{j_0}$                                                                            | scaled reaction current density per volume          |
| $\bar{I}$                         | $\bar{I}=\dfrac{I}{I_0}$                                                                            | scaled applied/electronic current                   |
| $\bar{U}_i$                       | $\bar{U}_i=\dfrac{U_i}{\phi_0}$                                                                     | scaled open-circuit potential of phase $i$          |
| $\bar{m}_i$                       | $`\bar{m}_i=\dfrac{m_i \ (c_{i,\max}\ c_{e0}^{0.5})}{j_0}`$                                         | scaled reaction-rate parameter in phase $i$, $\left[\dfrac{\mathrm{A\ m^{2.5}\ mol^{-1.5}}\left(\mathrm{mol\ m^{-3}}\right)^{3/2}}{\mathrm{A\ m^{-2}}}\right]$ |

## Parameters of the LG M50 battery cells [^2]

#### Negative electrode
| Code variable                | Parameter                     | Value               | Units                                    | Description                                  |
| ---------------------------- | ----------------------------- | ------------------  | ---------------------------------------- | -------------------------------------------- |
| cn0                          | $c_0$                         |               29866 | $\mathrm{mol}/\mathrm{m}^3$              | initial electrode lithium concentration      |
| cnmax                        | $c_{\text{max}}$              |               33133 | $\mathrm{mol}/\mathrm{m}^3$              | maximum lithium concentration                |
|                              |                               |                     |                                          |                                              |
| rn                           | $R_k$                         | $5.86\times10^{-6}$ | $\mathrm{m}$                             | particle radius, negative electrode          |
| ln                           | $L_k$                         | $85.2\times10^{-6}$ | $\mathrm{m}$                             | negative electrode thickness                 |
|                              |                               |                     |                                          |                                              |
| eps_n                        | $\epsilon_{\mathrm{e}}$       |                  25 | \%                                       | electrolyte volume fraction                  |
| eps_n_s                      | $\epsilon_k$                  |                  75 | \%                                       | active material volume fraction              |
|                              |                               |                     |                                          |                                              |
| kn_dim                       | $m_k$                         | $6.48\times10^{-7}$ | $\mathrm{A\ m}^{2.5}/\mathrm{mol}^{1.5}$ | reaction rate                                |
| Dn                           | $D_k$                         | $3.3\times10^{-14}$ | $\mathrm{m}^2/\mathrm{s}$                | electrode diffusivity                        |
| An                           | $a_k=3 \epsilon_{k}/R_k$      |                     | $\mathrm{m}^{-1}$                        | specific surface area of negative electrode  |
|                              |                               |                     |                                          |                                              |
|  sig_n                       | $\sigma_k$                    |                 215 | $\mathrm{S}/\mathrm{m}$                  | negative electrode conductivity              |
|                              |                               |                     |                                          |                                              |
| Un()                         | $U_{-}()$                     |                     | $\mathrm{V}$                             | OCV curve for the negative electrode         |

#### Positive electrode

| Code variable                | Parameter                     | Value               | Units                                    | Description                                  |
| ---------------------------- | ----------------------------- | ------------------  | ---------------------------------------- | -------------------------------------------- |
| cp0                          | $c_0$                         |               17038 | $\mathrm{mol}/\mathrm{m}^3$              | initial electrode lithium concentration      |
| cpmax                        | $c_{\text{max}}$              |               63104 | $\mathrm{mol}/\mathrm{m}^3$              | maximum lithium concentration                |
|                              |                               |                     |                                          |                                              |
| rp                           | $R_k$                         | $5.22\times10^{-6}$ | $\mathrm{m}$                             | particle radius, negative electrode          |
| lp                           | $L_k$                         | $75.6\times10^{-6}$ | $\mathrm{m}$                             | positive electrode thickness                 |
|                              |                               |                     |                                          |                                              |
| eps_p                        | $\epsilon_{\mathrm{e}}$       |                33.5 | \%                                       | electrolyte volume fraction                  |
| eps_p_s                      | $\epsilon_k$                  |                66.5 | \%                                       | active material volume fraction              |
|                              |                               |                     |                                          |                                              |
| kp_dim                       | $m_k$                         | $3.42\times10^{-6}$ | $\mathrm{A\ m}^{2.5}/\mathrm{mol}^{1.5}$ | reaction rate                                |
| Dp                           | $D_k$                         | $4.0\times10^{-15}$ | $\mathrm{m}^2/\mathrm{s}$                | electrode diffusivity                        |
| Ap                           | $a_k=3 \epsilon_{k}/R_k$      |                     | $\mathrm{m}^{-1}$                        | specific surface area of positive electrode  |
|                              |                               |                     |                                          |                                              |
|  sig_p                       | $\sigma_k$                    |                0.18 | $\mathrm{S}/\mathrm{m}$                  | negative electrode conductivity              |
|                              |                               |                     |                                          |                                              |
| Up()                         | $U_{+}()$                     |                     | $\mathrm{V}$                             | OCV curve for the positive electrode         |

#### Separator

| Code variable                | Parameter                     | Value               | Units                                    | Description                                  |
| ---------------------------- | ----------------------------- | ------------------  | ---------------------------------------- | -------------------------------------------- |
| ls                           | $L_k$                         | $12\times10^{-6}$   | $\mathrm{m}$                             | separator thickness                          |
| tplus                        | $t^{+}$                       | 0.2594              | -                                        | transference number of lithium ions          |
| eps_s                        | $\epsilon_{\mathrm{s}}$       |                  47 | \%                                       | electrolyte volume fraction, separator       |

#### Electrolyte

| Code variable                | Parameter                     | Value               | Units                                    | Description                                    |
| ---------------------------- | ----------------------------- | ------------------  | ---------------------------------------- | ---------------------------------------------- |
|  ce0                         | $c_{\mathrm{e}0}$             |   1000              | $\mathrm{mol}/\mathrm{m}^3$              | initial electrolyte concentration              |
|  De()                        | $`D_{\mathrm{e}}()`$          |                     | $\mathrm{m}^2/\mathrm{s}$                | diffusivity of lithium ions in the electrolyte |
|  kappa()                     | $`\sigma_{\mathrm{e}}()`$     |                     | $\mathrm{S}/\mathrm{m}$                  | electronic conductivity                        |

#### Cell parameters

| Code variable                | Parameter                     | Value               | Units                                    | Description                                  |
| ---------------------------- | ----------------------------- | ------------------  | ---------------------------------------- | -------------------------------------------- |
|  llayer                      |                               |  1.58               |  $\mathrm{m}$                            | wound layer length                           |
|  wlayer                      |                               |  $6.5\times10^{-2}$ |  $\mathrm{m}$                            | wound layer width                            |

#### Other parameters

| Code variable                | Parameter                     | Value               | Units                                    | Description                                  |
| ---------------------------- | ----------------------------- | ------------------  | ---------------------------------------- | -------------------------------------------- |
|  I_typ                       |                               |  5                  |  $\mathrm{A}$                            | typical current                              |
|  brugg                       |                               |  1.5                |  -                                       | Bruggeman coefﬁcient (theoretical value)     |


## SPMe model

According to [^1],
the averaged exchange current is

$$
\bar{j}_{i, \mathrm{ex}}=\frac{1}{\bar{L}_i} \int_0^{\bar{L}_i} \bar{m}_i \bar{c}_{\mathrm{e}}^{0.5} \bar{c}_{i, \text { surf }}^{0.5}\left(\bar{c}_{i, \text { surf }}-\bar{c}_{i, \max }\right)^{0.5} \mathrm{~d} \bar{L}
$$

the cell voltage is

$$
\bar{V}=\bar{U}_{\mathrm{p}}\left(\bar{c}_{\mathrm{p}, \text { surf }}\right)-\bar{U}_{\mathrm{n}}\left(\bar{c}_{\mathrm{n}, \text { surf }}\right)+2 \bar{T} \sinh ^{-1}\left(\frac{\bar{I}}{2 \bar{a}_{\mathrm{p}} \bar{L}_{\mathrm{p}} \bar{j}_{\mathrm{p}, \mathrm{ex}}}\right)-2 \bar{T} \sinh ^{-1}\left(\frac{-\bar{I}}{2 \bar{a}_{\mathrm{n}} \bar{L}_{\mathrm{n}} \bar{j}_{\mathrm{n}, \mathrm{ex}}}\right)-\bar{V}_{\mathrm{e}} .
$$

where the voltage drop in the electrode is approximated by

$$
\bar{V}_{\mathrm{e}}=-2 \bar{T}\left(1-t_{+}\right) \ln \frac{\bar{c}_{\mathrm{e}}(\bar{x}=0)}{\bar{c}_{\mathrm{e}}(\bar{x}=1)}+\frac{\bar{I}}{2}\left(\frac{\bar{L}_{\mathrm{n}}}{3 \bar{\kappa}_{\mathrm{n}}^{\text {eff }}}+\frac{\bar{L}_{\mathrm{sep}}}{\bar{\kappa}_{\mathrm{sep}}^{\text {eff }}}+\frac{\bar{L}_{\mathrm{p}}}{3 \bar{\kappa}_{\mathrm{p}}^{\text {eff }}}\right) \qquad (1)
$$
The logarithm can be approximated as follows
$$
\ln \bar{c}_{\mathrm{e}}(\bar{x}=0) = \ln (\pm c_{e0} + \bar{c}_{\mathrm{e}}(\bar{x}=0) )
= \ln \left(1 + \frac{- c_{e0} +\bar{c}_{\mathrm{e}}(\bar{x}=0)}{c_{e0}} \right) + \ln c_{e0}\approx \frac{- c_{e0} +\bar{c}_{\mathrm{e}}(\bar{x}=0)}{c_{e0}} + \ln c_{e0}
$$

$$
\bar{c}_{\mathrm{e}}(\bar{x}=0) \approx \mathrm{avg} (\bar{c}_{N E})
$$

$$
\ln \frac{\bar{c}_{\mathrm{e}}(\bar{x}=0)}{\bar{c}_{\mathrm{e}}(\bar{x}=1)} \approx  \frac{1}{c_{e0}}(\mathrm{avg} (\bar{c}_{N E}) - \mathrm{avg} (\bar{c}_{P E}))
$$

$$
\mathrm{avg}(f_{N E}) = \frac{1}{L N E}\int f_{N E} d x=\frac{1}{L N E}\int \bar{f}_{N E}(\bar{x}) \ \omega_{N E} d \bar{x} = \frac{N X}{N N E}\int \bar{f}_{N E}(\bar{x}) d \bar{x}
$$


According to [^3],
the terminal voltage $V$ can be written in the form

$$
V(t)=U_{\mathrm{eq}}-\eta_{\mathrm{r}}-\eta_{\mathrm{c}}-\Delta \phi_{\mathrm{e}}-\Delta \phi_{\mathrm{s}},
$$

where

$$
\eta_c = \left(1-t^{+}\right) \frac{2 R T}{F c_{\mathrm{e} 0}}\left(\frac{1}{L_{\mathrm{n}}} \int_0^{L_{\mathrm{n}}} c_{\mathrm{e}}(x, t) \mathrm{d} x-\frac{1}{L_{\mathrm{p}}} \int_{L-L_{\mathrm{p}}}^L c_{\mathrm{e}}(x, t) \mathrm{d} x\right) \qquad (2)
$$

which is equivalent to the first term in (1) with $\bar T = \dfrac{RT}{F}$ scaling.

*Code scaling follows equation (1), the averaging idea comes from equation (2).*

The other terms are defined as follows

$$
\eta_r = \frac{2 R T}{F}\left(\mathrm{arcsinh}\left(\frac{j_{\mathrm{n}}(t)}{j_{\mathrm{n} 0}(t)}\right)-\mathrm{arcsinh}\left(\frac{j_{\mathrm{p}}(t)}{j_{\mathrm{p} 0}(t)}\right)\right),
$$

$$
\Delta \phi_{\mathrm{e}} = \frac{i_{\mathrm{app}}}{\sigma_{\mathrm{e}}\left(c_{\mathrm{e} 0}\right)}\left(\frac{L_{\mathrm{n}}}{3 \mathcal{B}_{\mathrm{n}}}+\frac{L_{\mathrm{s}}}{\mathcal{B}_{\mathrm{s}}}+\frac{L_{\mathrm{p}}}{3 \mathcal{B}_{\mathrm{p}}}\right),
$$

$$
\Delta \phi_{\mathrm{s}} = \frac{i_{\text {app }}}{3}\left(\frac{L_{\mathrm{p}}}{\sigma_{\mathrm{p}}}+\frac{L_{\mathrm{n}}}{\sigma_{\mathrm{n}}}\right),
$$

$$
U_{\mathrm{eq}} = U_{\mathrm{p}}\left(\left.c_{\mathrm{p}}\right|_{r=R_{\mathrm{p}}}\right)-U_{\mathrm{n}}\left(\left.c_{\mathrm{n}}\right|_{r=R_{\mathrm{n}}}\right).
$$

## MFEM scaling

The problem is solved on a uniform grid $\bar{x} \in[0,1]$ with spacing $\Delta \bar{x}_{\mathrm{ref}}=1 / N X$.


In real physical space $x$, we have different normalized lengths $L N E, L S E P, L P E$ for different regions. Then the physical element sizes become:
- NE elements: $\Delta x_{N E}=L N E / N N E$
- SEP elements: $\Delta x_{S E P}=L S E P / N S E P$
- PE elements: $\Delta x_{P E}=L P E / N P E$

To address this, the physical coordinate is represented by a piecewise linear mapping $$x=F(\bar{x}),$$

$$
F(\bar{x})=
\begin{cases}
\frac{L N E \cdot N X}{N N E} \bar{x}, & \bar{x} \in\left[0,\frac{N N E}{N X}\right],\\
L N E+\frac{L S E P \cdot N X}{N S E P}\left(\bar{x}-\frac{N N E}{N X}\right), & \bar{x} \in\left[\frac{N N E}{N X},\frac{N N E+N S E P}{N X}\right],\\
L N E+L S E P+\frac{L P E \cdot N X}{N P E}\left(\bar{x}-\frac{N N E+N S E P}{N X}\right), & \bar{x} \in\left[\frac{N N E+N S E P}{N X},1\right],
\end{cases}
$$

with the piecewise constant Jacobian

$$
F'(\bar{x})=\omega_i=
\begin{cases}
\frac{L N E \cdot N X}{N N E}, &\quad i = N E,\\
\frac{L S E P \cdot N X}{N S E P}, &\quad i = S E P,\\
\frac{L P E \cdot N X}{N P E}, &\quad i = P E.
\end{cases}
$$

Let $\bar{\phi}=\phi \circ F = \phi(F(\cdot))$ and $\bar{v}=v \circ F$.
Additonally, $f_x$ denotes the derivative of $f$ with respect to $x$.
We then apply the change of variables to the integrals below

$$
\int \sigma \phi_x v_x d x=\int \sigma\left(\frac{1}{\omega_i} \bar\phi_{\bar{x}}\right)\left(\frac{1}{\omega_i} \bar v_{\bar{x}}\right) \omega_i d \bar{x}=\int \frac{\sigma}{\omega_i} \ \bar\phi_{\bar{x}} \bar v_{\bar{x}} d \bar{x},
$$

$$
\int f v(x) d x=\int f \bar{v}(\bar{x}) \ \omega_i d \bar{x}=\int f \omega_i \ \bar{v}(\bar{x}) d \bar{x},
$$

where $\sigma$ and $f$ are constants.

## References
[^1]: Ai, W., & Liu, Y. (2023). *Improving the convergence rate of Newman’s battery model using 2nd order finite element method*. Journal of Energy Storage, 67, 107512. https://doi.org/10.1016/j.est.2023.107512
[^2]: Chen et al., (2020). *Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models*. Journal of The Electrochemical Society, 167(8). https://iopscience.iop.org/article/10.1149/1945-7111/ab9050
[^3]: Brosa Planella et al., (2022). *A continuum of physics-based lithium-ion battery models reviewed*. Progress in Energy, 4(4), 042003. https://iopscience.iop.org/article/10.1088/2516-1083/ac7d31
