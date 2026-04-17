#ifndef modifiedKistlerAlphaContactAngleFvPatchScalarField_H
#define modifiedKistlerAlphaContactAngleFvPatchScalarField_H

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"

namespace Foam
{

class modifiedKistlerAlphaContactAngleFvPatchScalarField
:
    public alphaContactAngleTwoPhaseFvPatchScalarField
{
    // Required user inputs
    scalar thetaE_;        // equilibrium contact angle [deg]
    scalar thetaA_;        // advancing contact angle [deg]
    scalar thetaR_;        // receding  contact angle [deg]
    scalar muEffective_;   // dynamic viscosity [Pa.s]
    scalar sigma_;         // surface tension [N/m]

    // Helper functions
    static scalar degToRad(const scalar deg);
    static scalar radToDeg(const scalar rad);
    static scalar clampDeg(const scalar thetaDeg);

    scalar hoffmanFunctionRad(const scalar x) const;
    scalar inverseHoffmanRad(const scalar thetaRad) const;
    scalar modifiedKistlerThetaDeg(const scalar ucl) const;

public:

    TypeName("modifiedKistlerAlphaContactAngle");

    modifiedKistlerAlphaContactAngleFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    modifiedKistlerAlphaContactAngleFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    modifiedKistlerAlphaContactAngleFvPatchScalarField
    (
        const modifiedKistlerAlphaContactAngleFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    modifiedKistlerAlphaContactAngleFvPatchScalarField
    (
        const modifiedKistlerAlphaContactAngleFvPatchScalarField&
    );

    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>
        (
            new modifiedKistlerAlphaContactAngleFvPatchScalarField(*this)
        );
    }

    modifiedKistlerAlphaContactAngleFvPatchScalarField
    (
        const modifiedKistlerAlphaContactAngleFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new modifiedKistlerAlphaContactAngleFvPatchScalarField(*this, iF)
        );
    }

    virtual tmp<scalarField> theta
    (
        const fvPatchVectorField& Up,
        const fvsPatchVectorField& nHat
    ) const;

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
