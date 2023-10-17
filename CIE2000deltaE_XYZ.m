function E = CIE2000deltaE_XYZ( XYZ1, XYZ2, XYZ_w )

C = makecform( 'xyz2lab', 'WhitePoint', whitepoint( 'd65' ) );
Lab1 = applycform( XYZ1./repmat( XYZ_w, [size(XYZ1,1) 1] ), C );
Lab2 = applycform( XYZ2./repmat( XYZ_w, [size(XYZ1,1) 1] ), C );

E = CIE2000deltaE(Lab1,Lab2 );

end
