    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (!LTS)
        {
            #include "setDeltaT.H"

            // Update the mesh for topology change, mesh to mesh mapping
            mesh.update();

            runTime++;

            // Move the mesh
            mesh.move();
        }
	
        // --- Directed interpolation of primitive fields onto faces
        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        volScalarField rPsi("rPsi", 1.0/psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());

        // Make fluxes relative to mesh-motion
        if (mesh.moving())
        {
            phiv_pos -= mesh.phi();
            phiv_neg -= mesh.phi();
        }

	volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));
        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            interpolate(c, pos, T.name())*mesh.magSf()
        );
        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            interpolate(c, neg, T.name())*mesh.magSf()
        );

        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );
        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField a_pos("a_pos", ap/(ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf("aSf", am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        #include "centralCourantNo.H"

	if (LTS)
        {
            #include "setRDeltaT.H"
            runTime++;
        }

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

        surfaceVectorField phiUp
        (
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
          + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf()
        );

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg
        );

        // Make flux for pressure-work absolute
        if (mesh.moving())
        {
            phiEp += mesh.phi()*(a_pos*p_pos + a_neg*p_neg);
        }

        volScalarField muEff("muEff", rho*turbulence->nuEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        #include "rhoEqn.H"

        // --- Solve momentum
        #include "rhoUEqn.H"

        // --- Solve species
        #include "rhoYEqn.H"

        // --- Solve energy
        #include "rhoEEqn.H"

        turbulence->correct();
        thermophysicalTransport->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
