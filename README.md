# comparment_model_identification
This work is intended to identify the compartment model (which are prone to interpretations), instead of differential system (which may be hard to interpret), given a bunch of data.
The data considered was infact epidemic number data in South Korea, taken from https://sites.google.com/view/snuaric/data-service/covid-19/covid-19-data?authuser=0.
Next, while for identification of non-linear differential system, reminiscent of methods like SINDy, we need to calculate the derivative.
Instead of doing that, we can use Fourier transform to nullify the treatment of any derivative, and produce a rather simple objective function with a parameters.
This objective function can be easily expandible into minimization of L1 norm of a vector x, subjected to Ax=B.
In this way, we can get the underlying differential system with simple components that can be back traced towards a compartment model.
For the case of demo, the data is divided into 50 days interval to model only confirmed and recovered case daily numbers, and in this way the corresponding dynamics of the pandemic (represented by compartment model)
evolve over time, contrary to a typical variable parameter system approach. 
The corresponding demo is shown in demo.gif

