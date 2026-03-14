{
    reaction->correct();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (composition.solve(i))
        {
            volScalarField& rhoYi = rhoY[i];
            volScalarField& Yi = Y[i];

            surfaceScalarField rhoYi_pos(interpolate(rhoYi, pos, "Yi"));
            surfaceScalarField rhoYi_neg(interpolate(rhoYi, neg, "Yi"));

            surfaceScalarField phiYi(aphiv_pos*rhoYi_pos + aphiv_neg*rhoYi_neg);

            // --- Solve Yi
	    // operator splitting-1: Convective step
            solve(fvm::ddt(rhoYi) + fvc::div(phiYi));
            Yi.ref() = rhoYi()/rho();
            Yi.correctBoundaryConditions();
            rhoYi.boundaryFieldRef() == rho.boundaryField()*Yi.boundaryField();

	    // operator splitting-2: Diffusive/Reactive step
	    // only reactive step is solved here.
            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi) - fvc::ddt(rho, Yi)
                ==
                reaction->R(Yi)
            );
            YiEqn.solve("Yi");    
        }
    }
    composition.normalise();
}

