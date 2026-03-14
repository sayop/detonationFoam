solve
(
    fvm::ddt(rhoU)
  + fvc::div(phiUp)
);

U.ref() = rhoU()/rho();
U.correctBoundaryConditions();
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

