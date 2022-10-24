#ifndef PMT_porousMixture_H
#define PMT_porousMixture_H

#include "porousMediumFwd.H"
#include "fluidPhaseFwd.H"
#include "dispersionModel.H"
#include "phaseFractionFieldFwd.H"

#include <basicMultiComponentMixture.H>
#include <dictionary.H>
#include <tmp.H>
#include <volFields.H>
#include <labelFwd.H>
#include <PtrList.H>

namespace Foam
{
namespace Pmt
{

class porousMixture
:
    public basicMultiComponentMixture
{
public:
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Winconsistent-missing-override"
    TypeName("porousMixture");
    #pragma clang diagnostic pop

    porousMixture
    (
        const porousMedium& medium,
        const fluidPhase& phase,
        const dictionary& transportProperties
    );

    tmp<volTensorField> Deff
    (
        label speciesi,
        const volVectorField& U,
        const phaseFractionField& frac
    );

    tmp<volScalarField> retardation
    (
        label speciesi,
        const phaseFractionField& frac
    ) const;

private:
    const porousMedium& medium_;
    List<dimensionedScalar> Kd_;
    PtrList<dispersionModel> dispersionModels_;
};

}
}

#include "porousMixtureI.H"

#endif
