#include "porousMedium.H"
#include "phaseFractionField.H"

inline Foam::tmp<Foam::volTensorField> Foam::Pmt::porousMixture::Deff
(
    label speciesi,
    const volVectorField& U,
    const phaseFractionField& frac
)
{
    return dispersionModels_[speciesi].Deff(U, frac);
}

inline Foam::tmp<Foam::volScalarField> Foam::Pmt::porousMixture::retardation
(
    label speciesi,
    const phaseFractionField& frac
) const
{
    const auto& Kd = Kd_[speciesi];

    if (Kd.value() == 0)
    {
        return volScalarField::New("Rd", medium_.mesh(), dimensionedScalar{dimless, One});
    }

    return volScalarField::New("Rd", One + medium_.rs()*Kd/frac);
}
