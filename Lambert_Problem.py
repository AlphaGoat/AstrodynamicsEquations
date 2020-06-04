import numpy as np
import Stumpff from Astro_Eqtns

#Implementation of example 5_02 from Curtis's 'Orbital Mechanics for 
#Engineering Students'

rad2deg = np.pi/180
mu = 398600


def Lambert_Solver(R1,R2,del_t,string):
	'''Takes input parameters for position of sat in 3-D plane and time
    between observations and solves for velocity of object'''
    #Calculate distance for R1 and R2:
    r1 = np.linalg.norm(R1,ord='fro')
    r2 = np.linalg.norm(R2,ord='fro')
    #Take cross product of two position vectors
    cross12 = np.cross(R1,R2)
    #Calculate the angle between the two position vectors
    theta = np.acos(np.dot(R1,R2)/(r1*r2))
    #Determine whether orbit is prograde or retrograde:
    if string == 'pro':
        if cross12[3] <= 0:
            theta = 2*np.pi - theta
    elif string == 'retro':
        if cross12[3] >= 0
            theta = 2*np.pi - theta
    else:
        string = 'pro'
        print('\n ** Prograde trajectory assumed.\n')
       
    #Solving for coefficients of lagrange polynomial
    A = np.sin(theta*np.sqrt(r1*r2/(1-np.cos(theta)))) 
    
    #Determine approximately which z forces F(z,t) to change sign,
    z = -100
    y = lambda x: r1 + r2 + A * (z * Stumpff.S(z) - 1)/np.sqrt(Stumpff.C(z))
    while F(z,del_t) < 0:
        z = z + 0.1
    
    #Set an error tolerance and a limit on the number of iterations
    tol = 1e-8
    nmax = 5000

    #Iterate using Newton's Method until z is determined within
    #given error tolerance, or until the max number of iterations
    #have been completed
    ratio = 1
    n = 0
    while abs(ratio) > tol) && (n <= nmax):
        n += 1
        ratio = F(z,del_t)/dFdz(z)
        z = z - ratio

    #Report if the max number of iterations was exceeded
    if n >= nmax:
        print('\n\n **number of iterations exceeded')
        print('\n%s' % nmax)

    #Solve for lagrange coefficients:
    f = 1 - y(z)/r1
    g = A*np.sqrt(y(z)/mu)
    gdot = 1 - y(z)/r2

    #Solve for the first component of velocity
    v1 = 1/g*(R2 - f*R1)

    #Solve for the second component of velocity
    v2 = 1/g*(gdot*R2 - R1)

	#Compile State Vectors
	SV1 = np.array([R1,v1])
	SV2 = np.array([R2,v2])

    return SV1, SV2
     

def F(z,del_t):
    '''Function formed by Newton's method to solve for z'''
    global mu
    term1 = (y(z)/StumpffC(z))^1.5*StumpffS(z)
    term2 = A*np.sqrt(y(z))
    term3 = np.sqrt(mu)*del_t
    return term1 + term2 - term3 

def dFdz(z):
    '''Derivative of F wrt z'''
    if z == 0:
        term1 = (np.sqrt(2)/40)*y(0)^(3/2)
        term2 = (A/8)*(np.sqrt(y(0))+A*np.sqrt(1/(2*y(0))))
        return term1 + term2
    else:
        term1 = (y(z)/StumpffC(z))^(3/2)
        term2 = (1/(2*z))*(StumpffC(z)-(3/2)*(StumpffS(z)/StumpffC(z)))
        term3 = (3/4)*(StumpffS(z)^2/StumpffC(z))
        term4 = 3*StumpffS(z)/StumpffC(z)))*np.sqrt(y(z)
        term5 = A*np.sqrt(StumpffC(z)/y(z))
        return term1*(term2+term3)+(A/8)*(term4+term5)

