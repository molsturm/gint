Provide wrappings for various integral libraries and a simple interface to choose between them

Idea is to have a default library for Gaussians and Sturmians and FEs and ...
and have ways to select different ones

The upstream integral library should make the decisision which matrix data type to use (if possible)
or else this wrapper libarary should be able to make the choice

In here we should define the interface of the integral contractor engine 
and provide implementations for the available types of integrals (sturmians, gaussians, ...)
