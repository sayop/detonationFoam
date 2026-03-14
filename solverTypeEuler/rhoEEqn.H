solve
(
    fvm::ddt(rhoE)
  + fvc::div(phiEp)
  ==
    reaction->Qdot()
);

e = rhoE/rho - 0.5*magSqr(U);
e.correctBoundaryConditions();
thermo.correct();

rhoE.boundaryFieldRef() ==
    rho.boundaryField()*
    (
        e.boundaryField() + 0.5*magSqr(U.boundaryField())
    );

//- Compute pressure
p.ref() = rho()/psi();
p.correctBoundaryConditions();
rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();

Info<< "min/max(p) = "<< min(p).value() << ", " << max(p).value() << endl;
Info<< "min/max(T) = "<< min(T).value() << ", " << max(T).value() << endl;
