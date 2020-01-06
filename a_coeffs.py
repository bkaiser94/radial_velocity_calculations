"""
Created by Ben Kaiser (UNC-Chapel Hill) 2020-01-06

All of the coefficients that will be needed to produce the a_n coefficients themselves for use in other equations.

All are taken from Appendix A of Hurley et al. 2000, I'm just putting them in a nice python format, or at least my 
conception of that...

a_coeffs[i][0]=alpha
a_coeffs[i][1]=beta
a_coeffs[i][2]=gamma
a_coeffs[i][3]=eta



"""


###original attempt without assuming scientific notation
##a_coeffs=[
    ##[],#0th index is empty because they index from 1 and I don't want to shift all of the numbers
    ##[
        ##1.593890,
        ##2.053038,
        ##1.231226,
        ##2.327785
        ##],#1
    ##[
        ##2.706708,
        ##1.483131,
        ##5.772723,
        ##7.411230
        ##],#2
    ##[
        ##1.466143,
        ##-1.048442,
        ##-6.795374,
        ##-1.391127
        ##],#3
    ##[
        ##4.141960,
        ##4.564888,
        ##2.958542,
        ##5.571483
        ##],#4
    ##[
        ##3.426349,
        ##0.,
        ##0.,
        ##0.
        ##],#5
    ##[
        ##1.949814,
        ##1.758178,
        ##-6.008212,
        ##-4.470533
        ##],#6
    ##[
        ##4.903830,
        ##0.,
        ##0.,
        ##0.
        ##],#7
    ##[
        ##5.212154,
        ##3.166411,
        ##-2.750074,
        ##-2.271549
        ##],#8
    ##[
        ##1.312179,
        ##-3.294936,
        ##9.231860,
        ##2.610989
        ##],#9
    ##[
        ##8.073972,
        ##0.,
        ##0.,
        ##0.
        ##]#10
    ##]
        
    
    
    
    
    
    
    
    
#new attempt assuming scientific notation
a_coeffs=[
    [],#0th index is empty because they index from 1 and I don't want to shift all of the numbers
    [
        1.593890e3,
        2.053038e3,
        1.231226e3,
        2.327785e2
        ],#1
    [
        2.706708e3,
        1.483131e3,
        5.772723e2,
        7.411230e1
        ],#2
    [
        1.466143e2,
        -1.048442e2,
        -6.795374e1,
        -1.391127e1
        ],#3
    [
        4.141960e-2,
        4.564888e-2,
        2.958542e-2,
        5.571483e-3
        ],#4
    [
        3.426349e-1,
        0.,
        0.,
        0.
        ],#5
    [
        1.949814e1,
        1.758178e0,
        -6.008212e0,
        -4.470533e0
        ],#6
    [
        4.903830e0,
        0.,
        0.,
        0.
        ],#7
    [
        5.212154e-2,
        3.166411e-2,
        -2.750074e-3,
        -2.271549e-3
        ],#8
    [
        1.312179e0,
        -3.294936e-1,
        9.231860e-2,
        2.610989e-2
        ],#9
    [
        8.073972e-1,
        0.,
        0.,
        0.
        ]#10
    ]
        
    
    
