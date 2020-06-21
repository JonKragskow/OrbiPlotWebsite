import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt


def radial_p(Orb,r):
    # !!! Radial Wavefunctions of p orbitals
    if Orb=='2p': return 1./(2.*np.sqrt(6.))
 
    if Orb=='3p': return (4. - 2.*r/3.)/(9.*np.sqrt(6.))

    if Orb=='4p': return (20. - 5.*r + (r/2.)**2.)/(32.*np.sqrt(15.))

    if Orb=='5p': return (120. - 90.*(2.*r/5.) + 18.*(2*r/5.)**2. - (2.*r/5.)**3.)/(150.*np.sqrt(30.))

    if Orb=='6p': return (840. - 840.*r/3. + 252.*(r/3.)**2. - 28.*(r/3.)**3. + (r/3.)**4.)/(432.*np.sqrt(210.))

def radial_d(Orb,r): 
    # !!! Radial Wavefunctions of d orbitals
    if Orb=='3d': out = 1./(9.*np.sqrt(30.))

    if Orb=='4d': out = (6.-r/2.)/(96.*np.sqrt(5.))

    if Orb=='5d': out = (42. - 28.*r/5. + (2.*r/5.)**2.)/(150.*np.sqrt(70.))

    if Orb=='6d': out = (336. - 168.*(r/3.) + 24.*(r/3.)**2. - (r/3.)**3.)/(864.*np.sqrt(105.))
    return out

def p_z(n,r,f,c):
    return n*np.exp(r/n)*c*np.sqrt(np.pi/3.)/abs(f)

def dz2_z(n,r,f,c):
    return np.sqrt(r**2./3. + n**2 * c * np.exp(r/n) * np.sqrt(np.pi/5) * 1/(3*abs(f)))


def bound_fun_dz2_extra(r):

    Orb = '4d'
    c = 0.003
    n = float(Orb[0])
    f = radial_d(Orb,r)

    return r**2./3. - n**2 * c * np.exp(r/n) * np.sqrt(np.pi/5) * 1/(3*abs(f))

def bound_fun_p(r):

    Orb = '5p'
    c = 0.0003
    n = float(Orb[0])
    
    return r**2 - p_z(n,r,radial_p(Orb,r),c)**2

def bound_fun_dz2(r):

    Orb = '4d'
    c = 0.003
    n = float(Orb[0])

    return r**2 - dz2_z(n,r,radial_d(Orb,r),c)**2 


#p orbital r domains

#val = 25.8375464341635
#
#sol = optimize.fsolve(bound_fun_p,val,xtol = 0.00000000001)
#
#print('%15.13f '% (sol))
#print(bound_fun_p(25.8546178917617))


# d orbital r domains excluding dz2

val = 10
sol = optimize.fsolve(bound_fun_dz2_extra,val,xtol = 0.000000001)

print('%15.13f '% (sol))
print(bound_fun_dz2_extra(sol))
print(bound_fun_dz2_extra(10.0090911475380))


fig, (ax) = plt.subplots( 1, 1, sharex=False, sharey='all', figsize=(8,6))

r = np.linspace(0,20, 1000)

ax.plot(r,r**2 - dz2_z(3,r,radial_d('3d',r),0.003)**2 )
ax.plot(r,np.zeros(1000))

#plt.show()

# dx2 orbital r domains



