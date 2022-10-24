reactions.correct(theta);

forAll(composition.Y(), speciesi)
{
    auto& C = composition.Y(speciesi);
    auto Deff = composition.Deff(speciesi, U, theta);
    volScalarField Rd{composition.retardation(speciesi, theta)};
    const auto& R = reactions.reactionTerm(speciesi);

    fvScalarMatrix CEqn
    {
            fvm::ddt(Rd*theta, C)
          + fvm::div(phi, C, "div(phi,C)")
          - fvm::laplacian(theta*Deff/Rd, C, "laplacian(Deff,C)")
        ==
            Rd*theta*R
          + theta*fvOptions(C)
    };

    fvOptions.constrain(CEqn);
    CEqn.solve(mesh.solver("C"));
    fvOptions.correct(C);
}
