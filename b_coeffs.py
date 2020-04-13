"""
Created by Ben Kaiser (UNC-Chapel Hill) 2020-03-31

All of the coefficients that will be needed to produce the b_n coefficients themselves for use in other equations.

All are taken from Appendix A of Hurley et al. 2000, I'm just putting them in a nice python format, or at least my 
conception of that...


b_coeffs is a dictionary with strings given like "b39" or for primed coefficients "b'41"


"""

b_coeffs={
    "b9":[
        2.751631e3,
        3.557098e2,
        0,
        0],
    "b10":[
        -3.820831e-2,
        5.872664e-2,
        0,
        0],
    "b'11":[
        1.071738e2,
        -8.970339e1,
        -3.949739e1,
        0],
    "b12":[
        7.348793e2,
        -1.531020e2,
        -3.793700e1,
        0],
    "b'13":[
        9.219293e0,
        -2.005865e0,
        -5.561309e-1,
        0],
    
    "b39":[
        1.314955e2,
        2.009258e1,
        -5.143082e-1,
        -1.379140e0],
    "b40_static":[
        1.823973e1,
        -3.074559e0,
        -4.307878e0,
        0],
    "b'41":[
        2.327037e0,
        2.403445e0,
        1.208407e0,
        2.087263e-1],
    "b42":[
        1.997378e0,
        -8.126205e-1,
        0,
        0],
    "b43":[
        1.079113e-1,
        1.762409e-2,
        1.096601e-2,
        3.058818e-3],
    "b'44":[
        2.327409e0,
        6.901582e-1,
        -2.158431e-1,
        -1.084117e-1]
    }
