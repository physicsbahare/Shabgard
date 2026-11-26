# HST subgrid audit
Actual artifact 15136554839, run 45858906823, downloaded and inspected with Python. Fit local quadratic to SNR² within 100 mas of peak (grid spacing 40 mas), using design [1,x,y,x²,xy,y²] with x,y in grid units. Peak = grid maximum minus 40 H^-1 gradient; conditional covariance = -2*40² H^-1. These are approximate profile likelihood curvature errors, not calibrated astrometric errors.
| visit/filter | east,north residual from moving prior (mas) | conditional SD (mas) |
|---|---|---|
| 2010 F140W unregistered |101.54,43.20|6.99,6.61|
|2011 F125W|3.84,25.41|4.92,4.96|
|2011 F160W|-2.89,24.12|6.32,6.13|
|2012 F125W|-10.05,11.76|4.77,4.55|
|2012 F160W|-10.73,29.13|5.81,5.45|
2011 filter difference 6.73,1.30 mas; 2012 difference 0.68,-17.37 mas. The latter north difference exceeds twice the quadrature conditional curvature errors; shared registration and morphology not accounted for. Reference 2 in 2012 has east offsets 176–199 mas vs typical 85–132 mas for other refs; ref 16 north offsets -57 to +4 mas vs ref 11 +42 to +63 mas in same visit. These are evidence of reference-dependent offsets, not proof of motion. 2010 only refs 6 and 10, not sufficient for robust tie. Next remeasure with shared reference IDs and continuous centroid fitting, broaden stationary reference pool and test injection recovery; no parallax claim. No new FLT run launched in this audit.