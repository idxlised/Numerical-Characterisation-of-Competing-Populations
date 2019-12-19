# Numerical-Characterisation-of-Competing-Populations


**Problem Description and Importance**

The field of Synthetic Biology is an interdisciplinary branch that takes elements from Engineering, and Biology to both modify, and use existing biological systems to carry out an intended function.

The focal point of research within the field of _experimental_ Synthetic Biology is centred around the monitoring and characterisation of (multiple) modified bacterial species growing and interacting with one another in the same environment. These experiments allow scientists, and researchers alike to  gain a better understanding in the study of interspecies population dynamics. However, the study of such systems requires precise and non-invasive instrumentation and characterisation.

The scope of our capstone project focuses on the design and development of easy-to-use instrumentation that is able to not only accurately monitor and differentiate the cell densities of two different bacterial populations within a flask, but also characterise the results into usable data and parameters. It is this characterization that requires the use of predictive model systems that combine experiments with theoretical modelling to explain how population dynamics aries. It is as such that this project seeks to tackle the issue of characterising competitive interactions between two cell populations.

**Mathematical Formulation**

The problem outlined above can be modelled by the two variable Competitive Lotka-Volterra equation(s). The Competitive Lotka-Volterra equations are very similar to the Lotka-Volterra equations for predation in that the equation for each species has one term for self-interaction and one term for the  interaction with other species. This equation is shown as follows:

![img](https://latex.codecogs.com/gif.latex?%5Cfrac%7B%7Edx_%7B1%7E%7D%7D%7Bdt%7D%3Dr_%7B1%7Dx_%7B1%7D%20%28%201%20-%20%5Cfrac%7Bx_%7B1%7D&plus;%20%5Calpha%20_%7B12%7Dx_%7B2%7D%7E%7D%7BK_%7B1%7D%7D%20%29)

![img](https://latex.codecogs.com/gif.latex?%5Cfrac%7B%7Edx_%7B2%7E%7D%7D%7Bdt%7D%3Dr_%7B2%7Dx_%7B2%7D%20%28%201-%20%28%20%5Cfrac%7Bx_%7B2%7D&plus;%20%5Calpha%20_%7B21%7Dx_%7B1%7D%7E%7D%7BK_%7B2%7D%7D%29)

The equations above differ from the Lotka-Volterra equations by adding an additional term, ![img](https://latex.codecogs.com/gif.latex?%5Calpha%20_%7B12%7D%7Eand%7E%20%5Calpha%20_%7B21%7D), to account for the species&#39; interactions towards one another, respectively.

Here ![img](https://latex.codecogs.com/gif.latex?%5C%20x) is the size of the population at a given time, r is the inherent per-capita growth rate, and K is the population&#39;s carrying capacity. This model can be scaled up for systems with more than two species, where N is the total number of interacting species, as follows:

![img](https://latex.codecogs.com/gif.latex?%5Cfrac%7B%7Edx_%7Bi%7E%7D%7D%7Bdt%7D%3Dr_%7Bi%7Dx_%7Bi%7D%20%28%201-%20%5Cfrac%7B%20%5Csum%20_%7Bj%3D1%7D%5E%7BN%7D%20%5Calpha%20_%7Bij%7Dx_%7Bj%7D%7E%7D%7BK_%7Bi%7D%7D%29)

The general equation above can be constructed into an interaction matrix containing all ![img](https://latex.codecogs.com/gif.latex?%5Calpha) values. The coexisting equilibrium points can be found by isolating the point at which all derivatives are equal to zero, but that is not the origin. This can be found by inverting the interaction matrix and multiplying by the unit column vector. It is also important to note that there are always ![img](https://latex.codecogs.com/gif.latex?%5C2%5E%7BN%7D) equilibrium points, although the others often have at least one species&#39; population equal to zero.

**Problem Inputs and Outputs**

The inputs to the above model are characterised by ![img](https://latex.codecogs.com/gif.latex?%5C%20x_%7B1%7D%20%5C%20and%20%5C%20x_%7B2%7D), and the outputs are characterised by ![img](https://latex.codecogs.com/gif.latex?%5C%20r%2C%7E%20%5Calpha%20%2C) and ![img](https://latex.codecogs.com/gif.latex?%5C%20K) for both populations.

**Tentative Numerical Techniques**

In many cases, the process by which similar models have been used to predict the behaviour of important state variables to fit a solution of a differential equation to observed data has only been considered by few methods. It has predominantly been variations of the least squares regression method, though this is computationally exhaustive, especially for larger models. It is as such that we will be using an integration-based method, such as the one outlined by Kloppers et al. and Hodler et al.. This method will allow us to transform our system of ordinary differential equations to an algebraic system of equations that can be solved for the unknown parameters. This method is relatively robust and can be adapted to different and larger models.

**Numerical Experiments**

In order to verify that the numerical method above is capable of approximating the parameters of our model within the specified error limits, we will be running a series of iterative simulations in Simulink using known variables. The edge cases in which there are no interactions, or perhaps a large discrepancy in the initial concentrations, will be simulated and tested to further challenge the robustness of our model. If the model proves to be robust with a known dataset, then we can later continue to feed it with datasets embedded with simulated instrumentation errors. This will allow us to characterize the model&#39;s ability to handle external noise.
