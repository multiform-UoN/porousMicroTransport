#include "RichardsModel.H"
#include "unsaturatedFlowModelNew.H"

Foam::autoPtr<Foam::Pmt::RichardsModel> Foam::Pmt::RichardsModel::New
(
    const porousMedium& medium,
    const fluidPhase& phase,
    const dictionary& transportProperties
)
{
    return Detail::unsaturatedFlowModelNew<RichardsModel>(medium, phase, transportProperties);
}
