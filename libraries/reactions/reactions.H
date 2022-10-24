#ifndef PMT_reactions_H
#define PMT_reactions_H

#include "speciesCoeffs.H"
#include "porousMixture.H"
#include "phaseFractionFieldFwd.H"

#include <volFieldsFwd.H>
#include <label.H>
#include <PtrList.H>
#include <Pair.H>
#include <List.H>
#include <scalar.H>
#include <dimensionSet.H>

namespace Foam
{
namespace Pmt
{

class reactions
{
public:
    reactions
    (
        const porousMixture& composition,
        const dictionary& transportProperties
    );

    void correct(const phaseFractionField& frac);

    const volScalarField& reactionTerm(label speciesi) const;

    reactions(const reactions&) = delete;
    void operator=(const reactions&) = delete;

private:
    struct reaction
    {
        List<speciesCoeffs> lhs;
        List<speciesCoeffs> rhs;
        scalar kf;
        scalar kr;
    };

    const porousMixture& composition_;
    PtrList<volScalarField> reactionTerms_;
    List<reaction> reactions_;

    void clearTerms();
    void setReactionRate
    (
        const List<speciesCoeffs>& lhs,
        const List<speciesCoeffs>& rhs, 
        scalar k,
        const phaseFractionField& frac
    );

    static scalar reactionOrder(const List<speciesCoeffs>& lhs);

    static void checkYDimensions
    (
        const List<speciesCoeffs>& s,
        const basicMultiComponentMixture& composition,
        const dimensionSet& expected
    );

    static const dimensionSet& YDimensions
    (
        const List<speciesCoeffs>& lhs,
        const List<speciesCoeffs>& rhs,
        const basicMultiComponentMixture& composition
    );

    static dimensionSet kDimensions
    (
        const List<speciesCoeffs>& lhs,
        const List<speciesCoeffs>& rhs,
        const basicMultiComponentMixture& composition
    );
};

}
}

#include "reactionsI.H"

#endif
