function tformOut = normalizeTF(tformIn)

    if isa(tformIn, 'affine2d')
        tformOut = affinetform2d(tformIn.T');

    elseif isa(tformIn, 'affinetform2d') || ...
           isa(tformIn, 'images.geotrans.PolynomialTransformation2D') || ...
           isa(tformIn, 'images.geotrans.PiecewiseLinearTransformation2D') || ...
           isa(tformIn, 'images.geotrans.LocalWeightedMeanTransformation2D')
        tformOut = tformIn;

    else
        error('Unsupported transform type: %s', class(tformIn));
    end

end