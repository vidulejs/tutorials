// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    Info<< "---[pimpleFoam] (pimple.momentumPredictor()) is true, solving momentum equation" << endl;
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);

    Info<< "---[pimpleFoam] fvOptions.correct(U)" << endl;
}
else
{
    Info<< "---[pimpleFoam] (pimple.momentumPredictor()) is false, skipping" << endl;
}
