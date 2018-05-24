#--------------------------- Functions --------------------------------------------

import numpy as np

def ang_sep(lon1, lat1, lon2, lat2):
    """
    Angular separation between two points on a sphere

    Parameters
    ----------
    lon1, lat1, lon2, lat2 : Angle, Quantity or float
        Longitude and latitude of the two points.  Quantities should be in
        angular units; floats in radians

    Returns
    -------
    angular separation : Quantity or float
        Type depends on input; Quantity in angular units, or float in radians

    Notes
    -----
    The angular separation is calculated using the Vincenty formula [1]_,
    which is slightly more complex and computationally expensive than
    some alternatives, but is stable at at all distances, including the
    poles and antipodes.

    .. [1] http://en.wikipedia.org/wiki/Great-circle_distance
    """
    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)
    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon
    sep = np.arctan2(np.sqrt(num1 ** 2 + num2 ** 2), denominator)
    return sep



def eq2p2(RA,Dec,RA_center,Dec_center): 
    
    #center all points in RA
    RAprime = RA-RA_center
    RA_center=0
    #convert to positive values of RA
    negative = (RAprime<0)
    RAprime[negative] = 2*np.pi+RAprime[negative]
    #define wich quadrant is each point
    Qd1 = (RAprime<np.pi)&(Dec-Dec_center>0)
    Qd2 = (RAprime<np.pi)&(Dec-Dec_center<0)
    Qd3 = (RAprime>np.pi)&(Dec-Dec_center<0)
    Qd4 = (RAprime>np.pi)&(Dec-Dec_center>0)
    
    #Calculate the distance between the center and object, and the azimuthal angle
    Sep = ang_sep(RAprime, Dec, RA_center,Dec_center)
        
    #build a triangle to calculate the spherical cosine law
    x = ang_sep(RA_center,Dec,RA_center,Dec_center)
    y = ang_sep(RAprime,Dec, RA_center,Dec)
    
    #Apply shperical cosine law
    cosT = (np.cos(y) - np.cos(x)*np.cos(Sep))/(np.sin(x)*np.sin(Sep))
    #Round cosines that went over because of rounding errors
    roundinghigh = (cosT >  1)
    roundinglow  = (cosT < -1)
    cosT[roundinghigh]=1
    cosT[roundinglow]=-1
    uTheta = np.arccos(cosT)
    #Correct the angle for quadrant (because the above law calculates the acute angle between the "horizontal-RA" direction
    #and the angular separation great circle between the center and the object
    Theta= np.zeros(len(uTheta))
    
    Theta[Qd1] = uTheta[Qd1]
    Theta[Qd2] = np.pi-uTheta[Qd2] 
    Theta[Qd3] = np.pi+uTheta[Qd3]
    Theta[Qd4] = 2*np.pi-uTheta[Qd4]
   
    return Sep, Theta, cosT,uTheta
