import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import subprocess


def radial_p(Orb,r):
    # !!! Radial Wavefunctions of p orbitals
    if Orb=='2p': return 1./(2.*np.sqrt(6.))
 
    if Orb=='3p': return (4. - 2.*r/3.)/(9.*np.sqrt(6.))

    if Orb=='4p': return (20. - 5.*r + (r/2.)**2.)/(32.*np.sqrt(15.))

    if Orb=='5p': return (120. - 90.*(2.*r/5.) + 18.*(2*r/5.)**2. - (2.*r/5.)**3.)/(150.*np.sqrt(30.))

    if Orb=='6p': return (840. - 840.*r/3. + 252.*(r/3.)**2. - 28.*(r/3.)**3. + (r/3.)**4.)/(432.*np.sqrt(210.))

def radial_d(Orb,r): 
    # !!! Radial Wavefunctions of d orbitals
    if Orb == '3d_z2': 
        return 1./(9.*np.sqrt(30.))

    if Orb == '4d_z2':
        return (6.-r/2.)/(96.*np.sqrt(5.))

    if Orb == '5d_z2':
        return (42. - 28.*r/5. + (2.*r/5.)**2.)/(150.*np.sqrt(70.))

    if Orb == '6d_z2':
        return (336. - 168.*(r/3.) + 24.*(r/3.)**2. - (r/3.)**3.)/(864.*np.sqrt(105.))

def d_z_ax_pos(n,r,f,c):
    return np.sqrt(r**2/3 + (n**2 * c * np.exp(r/n) * np.sqrt(np.pi/5) * 1./(3.*abs(f))))

def d_z_ax_neg(n,r,f,c):
    return np.sqrt(r**2/3 - (n**2 * c * np.exp(r/n) * np.sqrt(np.pi/5) * 1./(3.*abs(f))))

def p_z(n,r,f,c):
    return n*np.exp(r/n)*c*np.sqrt(np.pi/3.)/abs(f)

def s_func(n,r,f,c)
    return n*np.exp(r/n)*c*np.sqrt(np.pi/3.)/abs(f)

def bound_fun_s(r, Orb, c):

    n = float(Orb[0])

    return r**2 - s_func(n,r,radial_s(Orb,r),c)**2

def bound_fun_p(r):

    Orb = '5p'
    c = 0.0003
    n = float(Orb[0])
    
    return r**2 - p_z(n,r,radial_p(Orb,r),c)**2

def bound_fun_dz2_pos(r, Orb, c):

    n = float(Orb[0])

    return r**2 - d_z_ax_pos(n,r,radial_d(Orb,r),c)**2

def bound_fun_dz2_neg(r, Orb, c):

    n = float(Orb[0])

    return np.nan_to_num(r**2 - d_z_ax_neg(n,r,radial_d(Orb,r),c)**2, 200000000000.)


fig, (ax) = plt.subplots( 1, 1)

# dz2 orbital r domains

c = 0.0003
orb = '1s'
n = float(Orb[0])
tolerance = 0.000000000000001

val_1 = 0.5
sol_1 = optimize.root_scalar(bound_fun_s, x0 = val_1, x1 = val_1+0.01, args = (orb, c), xtol = tolerance)

# print(sol_1, '\n')
# print(sol_2)
r = np.arange(58.47034306, 58.47034307, 0.000000000001)
 
cv = 58.470343060020994
print(cv)
print(d_z_ax_neg(n,cv,radial_d(orb,cv),c))
#Calculate z(r)
z_of_r_neg = d_z_ax_neg(n,r,radial_d(orb,r),c)

subprocess.run(["rm", "outfile"])

np.savetxt('outfile', np.transpose(np.vstack((r,z_of_r_neg))))

exit()

output = np.nan_to_num(r**2 - d_z_ax_neg(n,r,radial_d(orb,r),c)**2, 20000.)

ax.plot(r, np.zeros(np.size(r)))
ax.plot(r, output)
plt.show()

print(output.min())
