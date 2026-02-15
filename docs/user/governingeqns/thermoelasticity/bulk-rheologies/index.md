(sec-user-governingeqns-thermoelasticity-rheologies)=
# Bulk Rheologies

In this section we describe the mathematical formulations of bulk rheologies for thermoelasticity.
The thermal coupling introduces temperature-dependent material properties and thermal strain contributions to the constitutive relations.

## Linear Isotropic Thermoelasticity

The fundamental constitutive relation for linear isotropic thermoelasticity is:
%
\begin{equation}
\boldsymbol{\sigma} = \boldsymbol{C} : \boldsymbol{\epsilon} - \boldsymbol{C} : \alpha (T - T_0) \mathbf{I},
\end{equation}
%
where $\boldsymbol{C}$ is the elasticity tensor, $\alpha$ is the linear thermal expansion coefficient, and $T_0$ is a reference temperature.

For isotropic materials:
%
\begin{equation}
\boldsymbol{\sigma} = 2\mu \boldsymbol{\epsilon} + \lambda \mathop{\mathrm{Tr}}(\boldsymbol{\epsilon}) \mathbf{I} - (3\lambda + 2\mu) \alpha (T - T_0) \mathbf{I},
\end{equation}
%
where $\mu$ is the shear modulus and $\lambda$ is Lamé's first parameter.

## Temperature-Dependent Properties

Material properties may vary with temperature:
- Elastic moduli: $\mu(T)$, $\lambda(T)$
- Thermal expansion coefficient: $\alpha(T)$
- Thermal conductivity: $k(T)$
- Specific heat capacity: $c(T)$

For small temperature variations, linear approximations are often sufficient:
%
\begin{equation}
\mu(T) \approx \mu_0 \left[1 - \beta_\mu (T - T_0)\right],
\end{equation}
%
where $\beta_\mu$ is the temperature sensitivity coefficient.

## Thermal Strain

The thermal strain tensor for isotropic materials is:
%
\begin{equation}
\boldsymbol{\epsilon}^{\mathrm{th}} = \alpha (T - T_0) \mathbf{I}.
\end{equation}
%
The total strain is decomposed as:
%
\begin{equation}
\boldsymbol{\epsilon}^{\mathrm{total}} = \boldsymbol{\epsilon}^{\mathrm{elastic}} + \boldsymbol{\epsilon}^{\mathrm{th}}.
\end{equation}

## Thermoviscoelasticity

For viscoelastic materials with thermal coupling, temperature affects both elastic and viscous response:
- Elastic moduli depend on temperature
- Viscosity typically decreases with increasing temperature following Arrhenius-type relations:
%
\begin{equation}
\eta(T) = \eta_0 \exp\left(\frac{E_a}{R}\left(\frac{1}{T} - \frac{1}{T_0}\right)\right),
\end{equation}
%
where $E_a$ is the activation energy and $R$ is the gas constant.