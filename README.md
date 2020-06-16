# PHB model
Implementation of PHB kinetic model in Julia

The implementation here is for the model described in the paper <em>Fit-for-Purpose Models are Required for Competent Design of High-Volume Bioprocesses</em>. An excerpt from the paper describing the model is provided below.

## Model Description
The model is constructed with two underlying processes, growth and PHB production. The model includes substrate and ammonia inhibition, however, no analysis of the PHB consumption and biomass activity. The reactor contents simplified to substrate, ammonia, biomass and PHB are modelled. The pH, temperature, reactor headspace and fluid physical properties are not modelled.   The system considers the growth and accumulation phases. The mass stoichiometry for aerobic biomass production is modelled by:

<img src="https://latex.codecogs.com/svg.latex?&space;2C_{6}H_{12}O_{6}+0.15NH_{3}+0.72O_{2}\longrightarrow&space;Biomass+0.99CO_{2}+0.88H_{2}O" title=" 2C_{6}H_{12}O_{6}+0.15NH_{3}+0.72O_{2}\longrightarrow&space;Biomass+0.99CO_{2}+0.88H_{2}O" />

and the mass stoichiometry for PHB production from glucose by:

<img src="https://latex.codecogs.com/svg.latex?&space;1.67C_{6}H_{12}O_{6}+0.10O_{2}\longrightarrow&space;PHB+0.40CO_{2}+0.26H_{2}O" title=" 1.67C_{6}H_{12}O_{6}+0.10O_{2}\longrightarrow&space;PHB+0.40CO_{2}+0.26H_{2}O" />

The state variables modelled listed in the following table:

<table>
  <tr>
    <td>S<sub>s</sub></td>
    <td>Substrate (glucose) concentration</td>
    <td>(g/L)</td>
  </tr>
  <tr>
    <td>S<sub>NH</sub></td>
    <td>Ammonia concentration</td>
    <td>(g/L)</td>
  </tr>
  <tr>
    <td>X<sub>B</sub></td>
    <td>Biomass concentration</td>
    <td>(g/L)</td>
  </tr>
  <tr>
    <td>X<sub>PHB</sub></td>
    <td>PHB concentration</td>
    <td>(g/L)</td>
  </tr>
  <tr>
    <td>V</td>
    <td>Volume in reactor</td>
    <td>(L)</td>
  </tr>
 </table>

The process state is modelled according to:

<img src="https://latex.codecogs.com/svg.latex?&space;\frac{dS_s}{dt}=\frac{q_{in}}{V}(S_{s_{in}}-S_s)-2r_x-1.67r_p" />
<img src="https://latex.codecogs.com/svg.latex?&space;\frac{dS_{NH}}{dt}=\frac{q_{in}}{V}(S_{NH_{in}}-S_{NH})-0.15r_x" />
<img src="https://latex.codecogs.com/svg.latex?&space;\frac{dX_{B}}{dt}=\frac{q_{in}}{V}(X_{B_{in}}-X_B)+r_x" />
<img src="https://latex.codecogs.com/svg.latex?&space;\frac{dX_{PHB}}{dt}=\frac{q_{in}}{V}(X_{PHB_{in}}-X_{PHB})+r_p" />

where r<sub>x</sub> is the rate of biomass production and r<sub>p</sub> is the rate of PHB production. The volume in the reactor is modelled by:

<img src="https://latex.codecogs.com/svg.latex?&space;\frac{dV}{dt}=q_in-q_out" />

where q is the flow rate in L/d. The reaction rates r<sub>x</sub> and r<sub>p</sub> are calculated as follows:

<img src="https://latex.codecogs.com/svg.latex?&space;r_x=\mu_{max,g}X_B\frac{S_s}{K_{S,s,g}+S_s}\cdot\frac{S_{NH}}{K_{S,NH,g}+S_{NH}}" />
<img src="https://latex.codecogs.com/svg.latex?&space;r_x=\mu_{max,p}X_B\frac{S_s}{K_{S,s,g}+S_s}\cdot\frac{K_{I,NH,p}}{K_{I,NH,p}+S_{NH}}\cdot\frac{PHA_{max}}{PHA_{max}+X_{PHA}/(X_{PHA}+X_B)}\cdot\frac{K_{I,S}}{K_{I,S}+S_s}" />

The kinetic parameters applied in the model are representative of the microbial strain H. mediterranei and listed in the table below.

<table>
  <tr>
    <td>μ<sub>max,g</sub></td>
    <td>Maximum Biomass growth rate</td>
    <td>4.56</td>
    <td>d<sup>-1</sup></td>
  </tr>
  <tr>
    <td>μ<sub>max,p</sub></td>
    <td>Maximum PHA production rate</td>
    <td>3.6</td>
    <td>d<sup>-1</sup></td>
  </tr>
  <tr>
    <td>K<sub>S,s,g</sub></td>
    <td>Half saturation coefficient for growth</td>
    <td>2.98</td>
    <td>kg/m<sup>3</sup></td>
  </tr>
  <tr>
    <td>K<sub>S,s,p</sub></td>
    <td>Half saturation coefficient for PHA production</td>
    <td>9.3</td>
    <td>kg/m<sup>3</sup></td>
  </tr>
  <tr>
    <td>K<sub>S,NH,g</sub></td>
    <td>Ammonia saturation coefficient for growth</td>
    <td>0.021</td>
    <td>kg/m<sup>3</sup></td>
  </tr>
  <tr>
    <td>K<sub>I,NH,p</sub></td>
    <td>Ammonia inhibition coefficient for PHA uptake</td>
    <td>0.021</td>
    <td>kg/m<sup>3</sup></td>
  </tr>
  <tr>
    <td>K<sub>I,S</sub></td>
    <td>Substrate inhibition coefficient</td>
    <td>2.68</td>
    <td>kg/m<sup>3</sup></td>
  </tr>
  <tr>
    <td>PHA<sub>max</sub></td>
    <td>Max fraction of PHA content</td>
    <td>0.8</td>
    <td>-</td>
  </tr>
</table>
Data sourced from (Koller et al., 2006)

## Model Implementation
The model as described above is implemented using the Julia programming language. Julia is an innovative up and coming high-performance programming language with many features well suited to scientific and engineering computing applications. Julia’s state-of-the-art DifferentialEquations.jl library, which is both highly performant and feature rich, is particularly useful for creating kinetic models. While high-performance is not necessary for this implementation, it becomes important when applying the models to larger optimisation problems (i.e. optimising multiple parameters simultaneously) and applying the same model for online process optimisation. It was therefore interesting to explore this upcoming language here.

The Julia implementation showed an impressive performance advantage of up to 10x over the Python implementation (despite this being the first-time using Julia). The system of ordinary differential equations (ODEs) presented in above are stiff and are solved numerically using the Rodas4 method included in the DifferentialEquations.jl library. Rodas4 is a 5th order A-stable Rosenbrock method with a stiff-aware 3rd order interpolant (similar to ode15s in MATLAB).

### Reference
Koller, M., Horvat, P., Hesse, P., Bona, R., Kutschera, C., Atlić, A., Braunegg, G., 2006. Assessment of formal and low structured kinetic modeling of polyhydroxyalkanoate synthesis from complex substrates. Bioprocess Biosyst. Eng. 29, 367–377. https://doi.org/10.1007/s00449-006-0084-x
