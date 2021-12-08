pub enum FilterKind {
    LowPass,
    HighPass,
    BandPass,
    BandStop,
}

pub enum IIRFilterKind {
    Butterworth, // butter
    Chebyshev1,  // 'cheby1'
    Chebyshev2,  // 'cheby2'
    Elliptic,    //            - Cauer/elliptic: 'ellip'
    Bessel,      //Thomson: 'bessel'
}

pub enum FilterForm {
    Sos, //- second-order sections (recommended): 'sos'
    Ba,  //            - numerator/denominator (default)    : 'ba'
    Zpk, //            - pole-zero                          : 'zpk'
}
