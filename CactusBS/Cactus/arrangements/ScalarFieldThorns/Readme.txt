This folder contains six thorns for the non-linear simulation of the scalar fields in a curved spacetime, where the curved spacetime could be bbh background, mixed binary background, binary neutron stars background, or single comapct star background.


The thorn EvolveScalarFields is to evolve the scalar fields in the curved spacetime.

The thorn IDScalar field is to provide the initial profile of the scalar field, it could set the initial values of the scalar field with a non-vanishing initial momentum or non-vanishing initial scalar fields.

The thorn Matter_EJP is to monitor the total energy, total momentum, and total angular momentum of the matter, where we use it to compute these values of the scalar field.


The thorn ScalarAnalysis is to compute the radiation of the matter in terms of the energy-momentum of the matter, where we use it to compute the scalar radiation. For example, it could compute the radiated linear momentum, the radiated energy, and the radiated angular momentum.


The thorn TwopuncturesMatterzyp is to generate the proper initial values of the variables, where we have added the scalar fields into the Hamitonian constraint and let the momentum constraints are satisfied automatically with a special initial setup of the scalar fields.

The thorn TwoPuncturesMatterzypbaum is used to generate the initial parameters for the par file. 
