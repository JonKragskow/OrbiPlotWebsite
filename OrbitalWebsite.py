import numpy as np
import dash
import math
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import dash.dependencies as ddep
from subprocess import call
import json
from plotly.subplots import make_subplots
import dash_defer_js_import as dji
from scipy.special import sph_harm
from plotly.subplots import make_subplots


####################################################################################################
####################################### Radial Functions ###########################################
####################################################################################################

def p_z_ax(n,r,f,c):
    return n*np.exp(r/n)*c*np.sqrt(np.pi/3.)/abs(f)

def d_z_ax_pos(n,r,f,c): 
    return np.sqrt(r**2 * 1./3. + n**2 * c * np.exp(r/n) * np.sqrt(np.pi/5) * 1./(3.*np.abs(f)))

def d_z_ax_neg(n,r,f,c):
    return np.sqrt(r**2 * 1./3. - n**2 * c * np.exp(r/n) * np.sqrt(np.pi/5) * 1./(3.*np.abs(f)))

def radial_s(Orb,r):
    # !!! Radial Wavefunction of s orbitals
    if Orb=='1s': 
        return 2.*np.exp(-r)
 
def radial_p(Orb,r):
    # !!! Radial Wavefunctions of p orbitals
    if Orb=='2p': 
        return 1./(2.*np.sqrt(6.))
 
    if Orb=='3p': 
        return (4. - 2.*r/3.)/(9*np.sqrt(6.))

    if Orb=='4p': 
        return (20. - 5.*r + (r/2.)**2.)/(32.*np.sqrt(15.))

    if Orb=='5p': 
        return (120. - 90.*(2*r/5.) + 18.*(2*r/5.)**2. - (2.*r/5.)**3.)/(150.*np.sqrt(30.))

    if Orb=='6p': 
        return (840. - 840.*r/3. + 252.*(r/3.)**2. - 28.*(r/3.)**3. + (r/3.)**4.)/(432.*np.sqrt(210.))

############################ d Orbital functions #########################

def radial_d(Orb,r): 
    # !!! Radial Wavefunctions of d orbitals
    if Orb.find('3d') == 0: 
        return 1./(9.*np.sqrt(30.))

    if Orb.find('4d') == 0: 
        return (6.-r/2.)/(96.*np.sqrt(5.))

    if Orb.find('5d') == 0: 
        return (42. - 28.*r/5. + (2.*r/5.)**2.)/(150.*np.sqrt(70.))

    if Orb.find('6d') == 0: 
        return (336. - 168.*(r/3.) + 24.*(r/3.)**2. - (r/3.)**3.)/(864.*np.sqrt(105.))

############################f Orbital Wavefunction functions#########################

#def radial_f(Orb,r): 
#    if Orb=='4f': 

#
#    if Orb=='5f': 

#
#    if Orb=='6f': 

#    return out


####################################################################################################
########################################## Lobe bounds #############################################
####################################################################################################

def dz_domain(Orb):
    # !!! Returns an array containing the bounds
    # of r for each lobe of the requested dz2 orbital

    if Orb.find('3d') == 0:

        num_lobes = 2

        lobe_1_lower =  0.8347013549589228 
        lobe_1_upper = 19.84715220152567 

        ring_1_lower = 1.269074408199
        ring_1_upper = 16.749590243795

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper]])
        ring_domains = np.array([[ring_1_lower, ring_1_upper]])

    if Orb.find('4d') == 0:

        num_lobes = 4

        lobe_1_lower = 0.9717206679891478 
        lobe_1_upper = 10.947850494510153

        lobe_2_lower = 13.281145273906876  
        lobe_2_upper = 30.71438185359193 

        ring_1_lower = 1.506089675248620585
        ring_1_upper = 10.009091147538

        ring_2_lower = 15.13702469934
        ring_2_upper = 24.553806957999903

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[lobe_2_lower, lobe_2_upper]])
        ring_domains = np.array([[ring_1_lower, ring_1_upper],[ring_2_lower, ring_2_upper]])

    if Orb.find('5d') == 0:

        num_lobes = 6

        lobe_1_lower = 0.31412614081764606 
        lobe_1_upper = 10.755951991517627 

        lobe_2_lower = 11.018457454488443 
        lobe_2_upper = 23.74309091422983 

        lobe_3_lower = 24.49747212710316 
        lobe_3_upper = 64.21693289152138 

        ring_1_lower = 0.4549182747066
        ring_1_upper = 10.629142563719

        ring_2_lower = 11.15480678351499755
        ring_2_upper = 23.38146805779990

        ring_3_lower = 24.8948878699568
        ring_3_upper = 58.470343060020994

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],
                                 [lobe_2_lower, lobe_2_upper],
                                 [lobe_3_lower, lobe_3_upper]])
        ring_domains = np.array([[ring_1_lower, ring_1_upper],
                                 [ring_2_lower, ring_2_upper],
                                 [ring_3_lower, ring_3_upper]])

    if Orb.find('6d') == 0:

        num_lobes = 8

        lobe_1_lower = 0.3565338974235
        lobe_1_upper = 10.272518865906342 

        lobe_2_lower = 10.590710660732219 
        lobe_2_upper = 21.562394639355578 

        lobe_3_lower = 22.37410716607108
        lobe_3_upper = 38.69673874491944 

        lobe_4_lower = 40.56767789682293 
        lobe_4_upper = 85.63716821662403 

        ring_1_lower = 0.518213089425
        ring_1_upper = 10.120097300417532

        ring_2_lower = 10.75781337004409
        ring_2_upper = 21.175221605007522

        ring_3_lower = 22.8072828487972
        ring_3_upper = 37.8107854483668

        ring_4_lower = 41.59113348635129
        ring_4_upper = 77.715204335756

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],
                                 [lobe_2_lower, lobe_2_upper],
                                 [lobe_3_lower, lobe_3_upper],
                                 [lobe_4_lower, lobe_4_upper]])
        ring_domains = np.array([[ring_1_lower, ring_1_upper],
                                 [ring_2_lower, ring_2_upper],
                                 [ring_3_lower, ring_3_upper],
                                 [ring_4_lower, ring_4_upper]])

    return lobe_domains, ring_domains, num_lobes

def dxy_domain(Orb):

    if Orb.find('3dxy') == 0:

        num_lobes = 2

        lobe_1_lower = 0.25762972206931206
        lobe_1_upper = 30.2964334304816 

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper]])

    return lobe_domains, num_lobes


def r_domain(Orb):
    # !!! Returns an array containing the bounds of r for each lobe of the requested orbital

    if Orb == '1s':
        num_lobes = 2
        lobe_1_lower = 0.01
        lobe_1_upper = 0.05

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper]])

#####2p orbital#####
    if Orb == '2p':
        num_lobes = 2
        lobe_1_lower = 0.030542415534
        lobe_1_upper = 11.973154175074
       
        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper]])

#####3p orbital#####
    if Orb == '3p':
        num_lobes = 4
        lobe_1_lower = 0.0521008719026 
        lobe_1_upper =  5.6457958052515
        
        lobe_2_lower = 6.4019092339593
        lobe_2_upper = 20.7317063563592

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper]])
    
#####4p orbital#####

    if Orb=='4p':
        num_lobes = 6

        lobe_1_lower = 0.0791782710855
        lobe_1_upper = 5.0739651680262
        
        lobe_2_lower = 6.0725425741786
        lobe_2_upper = 12.8640838577464
        
        lobe_3_lower = 16.5661209026166
        lobe_3_upper = 30.3121031661214

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper],[ lobe_3_lower, lobe_3_upper]])
    
#####5p orbital#####

    if Orb=='5p':
        num_lobes = 8
    
        lobe_1_lower = 0.0105650461441 
        lobe_1_upper = 5.2893751552301
        
        lobe_2_lower = 5.4182579720594
        lobe_2_upper = 13.087803450358
        
        lobe_3_lower = 13.4960576912696
        lobe_3_upper = 25.8375464341634
        
        lobe_4_lower = 26.9002312729155
        lobe_4_upper = 64.0314459536619

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper],[ lobe_3_lower, lobe_3_upper],[ lobe_4_lower, lobe_4_upper]])
    
#####6p orbital#####
    
    if Orb=='6p':
        num_lobes = 10

        lobe_1_lower = 0.0138230698659
        lobe_1_upper = 5.1856979119147
        
        lobe_2_lower = 5.3499953596572
        lobe_2_upper = 12.5530614690789
    
        lobe_3_lower = 13.0502450037780
        lobe_3_upper = 23.6116572475096

        lobe_4_lower = 24.7635128466900 
        lobe_4_upper = 40.5262363941473
        
        lobe_5_lower = 43.0775334433834
        lobe_5_upper = 84.65464 # Fix me!!!
    
        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper],[ lobe_3_lower, lobe_3_upper],[ lobe_4_lower, lobe_4_upper],[ lobe_5_lower, lobe_5_upper]])

#####3d orbital#####

    if Orb == '3d':

        num_lobes = 2

        lobe_1_lower = 0.0305424
        lobe_1_upper = 11.973154-0.03054255
       
        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper]])


#####4d orbitals other than z2#####

    if Orb == '4d':

        lobe_1_lower = 0.9717206679892
        lobe_1_upper = 10.9478504945101
        
        lobe_1_lower = 13.2811452739069
        lobe_2_upper = 30.7143818535919

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper]])
    
#####5d orbitals other than z2#####

    if Orb=='5d':

        lobe_1_lower = 0.6156
        lobe_1_upper = 10.3787
        
        lobe_2_lower = 11.3617
        lobe_2_upper = 12.1429
        
        lobe_3_lower = 25.5091
        lobe_3_upper = 33.2447

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper],[ lobe_3_lower, lobe_3_upper]])
    
#####4f orbital#####

    if Orb == '4f':
        lobe_1_lower = 0.0305424
        lobe_1_upper = 11.973154-0.03054255
       
        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper]])

#####5f orbital#####

    if Orb == '5f':

        lobe_1_lower = 0.0521075737439532  
        lobe_1_upper =  5.660956689880348059366305 - 0.0521075737439532
        
        lobe_2_lower = 6.382478125927591589829610
        lobe_2_upper = 20.94330871272566767441322 - 6.382478125927591589829610

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper]])
    
#####6f orbital#####

    if Orb=='6f':

        lobe_1_lower = 0.0791834876755957
        lobe_1_upper =  4.99478167864062
        
        lobe_2_lower = 6.07254257647759
        lobe_2_upper =  6.79154128088496
        
        lobe_3_lower = 16.5661209046895
        lobe_3_upper =  13.745982260283

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper],[ lobe_3_lower, lobe_3_upper]])


    return lobe_domains, num_lobes

####################################################################################################
###################################### Orbital Mathematics #########################################
####################################################################################################

def set_3d_colour(colour_name):

    if colour_name == 'go':

        colours = [[0, 'rgb(48,148,33)'], [1, 'rgb(255,185,0)']]

    elif colour_name == 'yp':

        colours = [[0, 'rgb(109,0,157)'], [1, 'rgb(255,204,51)']]

    elif colour_name == 'rb':

        colours = [[0, 'rgb(12,116,235)'], [1, 'rgb(255,0,0)']]

    return colours

def get_orb_name(orb):

    n = int(orb[0])
    l = orb[1]

    return n,l

def OrbCalc(orbital_input, colour_name, fig, cutaway):

    # Get colours of lobes
    colours = set_3d_colour(colour_name)
    
    # Get orbital n value and name
    n, l = get_orb_name(orbital_input)

    print(cutaway, flush=True)

    if l == 's':

        data = plot_s_orb(n, orbital_input, colours, fig, cutaway)

    elif l == 'p':

        data = plot_p_orb(n, orbital_input, colours, fig, cutaway)

    elif l == 'd':  

        print(orbital_input, orbital_input.find('xy'), flush=True)

        if orbital_input.find('xy') != -1:

            data = plot_dxy_orb(n, orbital_input,colours, fig, cutaway)

        else:

            data = plot_dz_orb(n, orbital_input,colours, fig, cutaway)

    elif l == 'f':

        data = plot_f_orb(n, colour_1, colour_2)

    return data

def plot_s_orb(n, orbital_input,  colours, fig):

    # Set contour level
    if n < 5:
        c = 0.003
    else:
        c = 0.0003

    # Set num steps for angle and r
    r_mini_steps = 20
    r_steps = 4*r_mini_steps
    angle_steps = 80

    #Array of angle values
    ang = np.linspace(-np.pi/2, np.pi, num = angle_steps)

    # Get bounds of r for each lobe of orbital
    orb_r_bounds, num_lobes = r_domain(orbital_input)

    data = []

    flag = True

    #Calculate coordinates of isosurface
    # loop over each lobe of the orbital
    
    num_shells = n

    for shell in np.arange(num_shells):

        low_bound = orb_r_bounds[shell,0]

        gap = np.abs(orb_r_bounds[shell,0] - orb_r_bounds[shell,1])

        x, y, z = calc_s_orb(n, c, orbital_input, gap, ang, orb_r_bounds, num_lobes, angle_steps, r_steps, r_mini_steps, low_bound)

        fig.add_trace(go.Surface(x=x, y=y, z=z[0,:,:], colorscale= colours, showscale=False), 1, 1)
        fig.add_trace(go.Surface(x=x, y=y, z=z[1,:,:], colorscale= colours, showscale=False), 1, 1)


    return fig

def calc_s_orb(n, c, orbital_input, gap, ang, orb_r_bounds, num_lobes, angle_steps, r_steps, r_mini_steps, low_bound):

    #Arrays for x,y,z coordinates
    # z has two separate arrays as a p orbital has two halves separated by a plane
    x = np.zeros([angle_steps, r_steps])
    y = np.zeros([angle_steps, r_steps])
    x_y_no_angle = np.zeros([r_steps])
    z = np.zeros([num_lobes, angle_steps, r_steps])

    #Array of r values for each lobe
    # Sample more frequently closer to the pole of the lobe
    r = np.linspace(low_bound, low_bound + gap*0.025, r_mini_steps)
    r = np.append(r,np.linspace(low_bound+ gap*0.025, low_bound + gap*0.5, r_mini_steps))
    r = np.append(r,np.linspace(low_bound+ gap*0.5, low_bound + gap*0.975, r_mini_steps))
    r = np.append(r,np.linspace(low_bound+ gap*0.975, low_bound + gap, r_mini_steps))

    #Calculate z(r)
    z_of_r = r

    #Calculate angle independent term of x and y
    x_y_no_angle = np.sqrt(r**2. - z_of_r**2.)

    #Calculate x and y by multiplying by angular factor
    #And set z coordinate for every angle and r value
    for r_it in np.arange(r_steps):
        for ang_it in np.arange(angle_steps):
            x[ang_it,r_it]   = z_of_r[r_it]*np.cos(ang[ang_it])
            y[ang_it,r_it]   = z_of_r[r_it]*np.sin(ang[ang_it])
            z[0,ang_it,r_it] = z_of_r[r_it]
            z[1,ang_it,r_it] = -z_of_r[r_it]

    return x, y, z


def plot_p_orb(n, orbital_input,  colours, fig, cutaway):

    # Set contour level
    if n < 5:
        c = 0.003
    else:
        c = 0.0003

    # Set num steps for angle and r
    r_mini_steps = 20
    r_steps = 4*r_mini_steps
    angle_steps = 50

    #Array of angle values
    ang = np.linspace(-np.pi * cutaway, np.pi, num = angle_steps)

    # Get bounds of r for each lobe of orbital
    orb_r_bounds, num_lobes = r_domain(orbital_input)

    data = []

    flag = True

    #Calculate coordinates of isosurface
    # loop over each lobe of the orbital

    num_shells = n-1
    sw = 0

    for shell in np.arange(num_shells):

        xt, yt, zt = calc_p_orb(n, c, orbital_input, orb_r_bounds[shell, :], ang, num_lobes, angle_steps, 
                             r_steps, r_mini_steps)

        if sw == 0 :
            cot1 = np.zeros(np.shape(zt))
            cot2 = np.ones(np.shape(zt))
            sw = 1
        else:
            cot2 = np.zeros(np.shape(zt))
            cot1 = np.ones(np.shape(zt))
            sw = 0

        if shell == 0:
            x = np.copy(xt)
            y = np.copy(yt)
            z = np.copy(zt)
            co = np.copy(cot1)
        else:
            x = np.append(x, xt, axis=0)
            y = np.append(y, yt, axis=0)
            z = np.append(z, zt, axis=0)
            co = np.append(co, cot1, axis=0)

        x = np.append(x, xt, axis=0)
        y = np.append(y, yt, axis=0)
        z = np.append(z, -zt, axis=0)
        co = np.append(co, cot2, axis=0)

    fig.add_trace(go.Surface(x=x, y=y, z=z, surfacecolor = co, colorscale=colours, showscale=False), 1, 1)

    return fig

def calc_p_orb(n, c, orbital_input, bounds, ang, num_lobes, angle_steps, 
               r_steps, r_mini_steps):

    print(bounds, flush=True)

    gap = np.abs(bounds[1] - bounds[0])

    #Array of r values for each lobe
    # Sample more frequently closer to the pole of the lobe
    r = np.linspace(bounds[0], bounds[0] + gap*0.015, r_mini_steps)
    r = np.append(r,np.linspace(bounds[0]+ gap*0.015, bounds[0] + gap*0.5, r_mini_steps))
    r = np.append(r,np.linspace(bounds[0]+ gap*0.5, bounds[0] + gap*0.975, r_mini_steps))
    r = np.append(r,np.linspace(bounds[0]+ gap*0.975, bounds[0] + gap, r_mini_steps))

    angm, rm = np.meshgrid(ang, r)

    zor = p_z_ax(n,rm,radial_p(orbital_input,rm),c)

    x = np.sqrt(rm**2. - zor**2.)*np.cos(angm)
    y = np.sqrt(rm**2. - zor**2.)*np.sin(angm)
    z = zor

    return x, y, z

def plot_dz_orb(n, orbital_input, colours, fig, cutaway):


    if n == 3:
        x,y,z = np.mgrid[-30:30:60j, -30:30*cutaway:60j, -30:30:60j]
    elif n == 4:
        x,y,z = np.mgrid[-42:42:60j, -42:42*cutaway:60j, -42:42:60j]
    elif n == 5:
        x,y,z = np.mgrid[-65:65:60j, -65:65*cutaway:60j, -65:65:60j]
    elif n == 6:
        x,y,z = np.mgrid[-90:90:110j, -90:90*cutaway:110j, -90:90:110j]

    r = np.sqrt(x**2 + y**2 + z**2)

    rho = 2*r/n

    if n == 3:
        Rad = 1./(9.*np.sqrt(30.))*rho**2.*np.exp(-rho/2.)
    elif n == 4:
        Rad = 1./(96.*np.sqrt(5.))*(6.-rho)*rho**2.*np.exp(-rho/2.)
    elif n == 5:
        Rad = 1./(150.*np.sqrt(70.))*(42.-14.*rho+rho**2)*rho**2.*np.exp(-rho/2.)
    elif n == 6:
        Rad = 1./(864.*np.sqrt(105.))*(336.-168.*rho+24.*rho**2.-rho**3.)*rho**2.*np.exp(-rho/2.)

    Ang = 2*z**2-x**2-y**2

    wav = Rad*np.sqrt(1/(4*np.pi)*np.sqrt(10)/r**2)*Ang

    print(np.shape(x), flush=True)

    if n == 3:
        ival=0.04
    elif n == 4:
        ival = 0.04
    elif n == 5:
        ival = 0.04
    elif n == 6:
        ival = 0.04

    print(ival, flush=True)


    fig.add_trace(go.Isosurface(
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        value=wav.flatten(),
        isomin=-ival,
        isomax=ival,
        caps=dict(x_show=False, y_show=False, z_show=False),
        showscale=False,
        colorscale=colours
        ))
    return fig


def plot_dxy_orb(n, orbital_input, colours, fig, cutaway):


    if n == 3:
        x,y,z = np.mgrid[-25:25:45j, -25:25:45j, -25:25*cutaway:45j]
    elif n == 4:
        x,y,z = np.mgrid[-42:42:50j, -42:42:50j, -42:42*cutaway:50j]
    elif n == 5:
        x,y,z = np.mgrid[-60:60:50j, -60:60:50j, -60:60*cutaway:50j]
    elif n == 6:
        x,y,z = np.mgrid[-77:77:60j, -77:77:60j, -77:77*cutaway:60j]

    r = np.sqrt(x**2 + y**2 + z**2)

    rho = 2*r/n

    if n == 3:
        Rad = 1./(9.*np.sqrt(30.))*rho**2.*np.exp(-rho/2.)
    elif n == 4:
        Rad = 1./(96.*np.sqrt(5.))*(6.-rho)*rho**2.*np.exp(-rho/2.)
    elif n == 5:
        Rad = 1./(150.*np.sqrt(70.))*(42.-14.*rho+rho**2)*rho**2.*np.exp(-rho/2.)
    elif n == 6:
        Rad = 1./(864.*np.sqrt(105.))*(336.-168.*rho+24.*rho**2.-rho**3.)*rho**2.*np.exp(-rho/2.)

    Ang = x*y

    wav = Rad*np.sqrt(1/(4*np.pi)*np.sqrt(10)/r**2)*Ang

    print(x)

    if n == 3:
        ival=0.005
    elif n == 4:
        ival = 0.01
    elif n == 5:
        ival = 0.01
    elif n == 6:
        ival = 0.01
        
    fig.add_trace(go.Isosurface(
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        value=wav.flatten(),
        isomin=-ival,
        isomax=ival,
        caps=dict(x_show=False, y_show=False, z_show=False),
        showscale=False,
        colorscale=colours
        ))
    return fig


def calc_dz_orb(n, c, orbital_input, bounds, ang, num_lobes, angle_steps, r_steps,
                r_mini_steps, pos):

    gap = np.abs(bounds[1] - bounds[0])

    #Array of r values for each lobe
    # Sample more frequently closer to the pole of the lobe
    r =             np.linspace(bounds[0], bounds[0] + gap*0.025, r_mini_steps)
    r = np.append(r,np.linspace(bounds[0]+ gap*0.025, bounds[0] + gap*0.5, r_mini_steps))
    r = np.append(r,np.linspace(bounds[0]+ gap*0.5, bounds[0] + gap*0.975, r_mini_steps))
    r = np.append(r,np.linspace(bounds[0]+ gap*0.975, bounds[0] + gap, r_mini_steps))

    angm, rm = np.meshgrid(ang, r)

    #Calculate z(r)
    if pos :
        zor = d_z_ax_pos(n,rm,radial_d(orbital_input,rm),c)
    else :
        zor = d_z_ax_neg(n,rm,radial_d(orbital_input,rm),c)

    x = np.sqrt(rm**2. - zor**2.)*np.cos(angm)
    y = np.sqrt(rm**2. - zor**2.)*np.sin(angm)
    z = zor

    return x, y, z

def calc_dxy_orb(n, c, orbital_input, bounds, ang, num_lobes, angle_steps, r_steps,
                r_mini_steps):

    gap = np.abs(bounds[1] - bounds[0])

    #Array of r values for each lobe
    # Sample more frequently closer to the pole of the lobe
    r =             np.linspace(bounds[0], bounds[0] + gap*0.025, r_mini_steps)
    r = np.append(r,np.linspace(bounds[0]+ gap*0.025, bounds[0] + gap*0.5, r_mini_steps))
    r = np.append(r,np.linspace(bounds[0]+ gap*0.5, bounds[0] + gap*0.975, r_mini_steps))
    r = np.append(r,np.linspace(bounds[0]+ gap*0.975, bounds[1], r_mini_steps))

    angm, rm = np.meshgrid(ang, r)

    #Calculate g(r)
    g_r = rm**2. - c * float(n)**2. * np.exp(rm/float(n)) * np.sqrt(np.pi/15.) / np.abs(radial_d(orbital_input,rm))

    x = np.sqrt(g_r/2.)*np.cos(angm)
    z = np.sqrt(g_r)*np.sin(angm)
    y = np.sqrt(rm**2. -x**2. -z**2.)
    
    return x, y, z


def swap(colour_1, colour_2):

    return colour_2, colour_1

###Functions for Radial Wave Functions and probability density functions

def s_1(r,mode, zeff):
    rho = np.copy(2.*r/1.  * zeff)
    O1s = 2.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O1s**2.
    if mode==2:
        return O1s

def s_2(r,mode, zeff):
    rho = np.copy(2.*r/2. * zeff)
    O2s = 1./(2.*np.sqrt(2.))*(2.-rho)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O2s**2.
    if mode==2:
        return O2s

def s_3(r,mode, zeff):
    rho = np.copy(2.*r/3. * zeff)
    O3s = 1./(9.*np.sqrt(3.))*(6.-6.*rho+rho**2.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O3s**2.
    if mode==2:
        return O3s

def s_4(r,mode, zeff):
    rho = np.copy(2.*r/4. * zeff)
    O4s = (1./96.)*(24.-36.*rho+12.*rho**2.-rho**3.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O4s**2.
    if mode==2:
        return O4s

def s_5(r,mode, zeff):
    rho = np.copy(2.*r/5. * zeff)
    O5s = (1./(300.*np.sqrt(5.)))*(120.-240.*rho+120.*rho**2.-20.*rho**3.+rho**4.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O5s**2.
    if mode==2:
        return O5s

def s_6(r,mode, zeff):
    rho = np.copy(2.*r/6. * zeff)
    O6s = (1./(2160.*np.sqrt(6.)))*(720.-1800.*rho+1200.*rho**2.-300.*rho**3.+30.*rho**4.-rho**5.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O6s**2.
    if mode==2:
        return O6s

def p_2(r,mode, zeff):
    rho = np.copy(2.*r/2. * zeff)
    O2p = 1./(2.*np.sqrt(6.))*rho*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O2p**2.
    if mode==2:
        return O2p

def p_3(r,mode, zeff):
    rho = np.copy(2.*r/3. * zeff)
    O3p = 1./(9.*np.sqrt(6.))*rho*(4.-rho)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O3p**2.
    if mode==2:
        return O3p

def p_4(r,mode, zeff):
    rho = np.copy(2.*r/4. * zeff)
    O4p = 1./(32.*np.sqrt(15.))*rho*(20.-10.*rho+rho**2.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O4p**2.
    if mode==2:
        return O4p

def p_5(r,mode, zeff):
    rho = np.copy(2.*r/5. * zeff)
    O5p = 1./(150.*np.sqrt(30.))*rho*(120.-90.*rho+18.*rho**2.-rho**3.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O5p**2.
    if mode==2:
        return O5p

def p_6(r,mode, zeff):
    rho = np.copy(2.*r/6. * zeff)
    O6p = 1./(432.*np.sqrt(210.))*rho*(840.-840.*rho+252.*rho**2.-28.*rho**3.+rho**4.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O6p**2.
    if mode==2:
        return O6p

def d_3(r,mode, zeff):
    rho =np.copy(2.*r/3. * zeff)
    O3d = 1./(9.*np.sqrt(30.))*rho**2.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O3d**2.
    if mode==2:
        return O3d

def d_4(r,mode, zeff):
    rho = np.copy(2.*r/4. * zeff)
    O4d = 1./(96.*np.sqrt(5.))*(6.-rho)*rho**2.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O4d**2.
    if mode==2:
        return O4d

def d_5(r,mode, zeff):
    rho = np.copy(2.*r/5. * zeff)
    O5d = 1./(150.*np.sqrt(70.))*(42.-14.*rho+rho**2)*rho**2.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O5d**2.
    if mode==2:
        return O5d

def d_6(r,mode, zeff):
    rho = np.copy(2.*r/6. * zeff)
    O6d = 1./(864.*np.sqrt(105.))*(336.-168.*rho+24.*rho**2.-rho**3.)*rho**2.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O6d**2.
    if mode==2:
        return O6d

def f_4(r,mode, zeff):
    rho = np.copy(2.*r/4. * zeff)
    O4f = 1./(96.*np.sqrt(35.))*rho**3.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O4f**2.
    if mode==2:
        return O4f

def f_5(r,mode, zeff):
    rho = np.copy(2.*r/5. * zeff)
    O5f = 1./(300.*np.sqrt(70.))*(8.-rho)*rho**3.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O5f**2.
    if mode==2:
        return O5f

def f_6(r,mode, zeff):
    rho = np.copy(2.*r/5. * zeff)
    O6f = 1./(2592.*np.sqrt(35.))*(rho**2.-18.*rho+72.)*rho**3.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O6f**2.
    if mode==2:
        return O6f


functiondict = {'1s': s_1,
                '2s': s_2, 
                '3s': s_3, 
                '4s': s_4, 
                '5s': s_5, 
                '6s': s_6, 
                '2p': p_2,
                '3p': p_3,
                '4p': p_4,
                '5p': p_5,
                '6p': p_6, 
                '3d': d_3, 
                '4d': d_4, 
                '5d': d_5, 
                '6d': d_6, 
                '4f': f_4, 
                '5f': f_5, 
                '6f': f_6}


##################################################################################################################################
############################################################ Webpage layout ######################################################
##################################################################################################################################

#Build the webpage layout

mathjax_script = dji.Import(src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/latest.js?config=TeX-AMS-MML_SVG")
refresh_plots =  dji.Import(src="https://codepen.io/chrisvoncsefalvay/pen/ExPJjWP.js")

app = dash.Dash(__name__)


# app.css.append_css({"external_url": "https://max.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css"})

# Change header 
# %thing% are defined by dash and are replaced at runtime
app.index_string = r'''
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>Orbiplot</title>
        {%favicon%}
        {%css%}
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            <script type="text/x-mathjax-config">
            MathJax.Hub.Config({
                tex2jax: {
                inlineMath: [ ['$','$'],],
                processEscapes: true
                }
            });
            </script>
            {%renderer%}
        </footer>
    </body>
</html>
'''

orbital_plot_options = [html.Div(className = "container", 
         style = {
                  'display' : 'grid',
                  'grid-template-columns': r'50% 50%',
                  'grid-template-rows' : r'50% 50%'
                  },
         children=[
                    html.Div(className = "item", 
                        style = {
                                 'grid-column-start': '1',
                                 'grid-column-end': '2',
                                 'grid-row-start': '1',
                                 'grid-row-end': '2'},
                        children=[
                                  html.H2(style = {
                                                  'textAlign' : 'center', 
                                                  'fontFamily' : 'sans-serif'
                                                 },
                                         children = 'Orbital'),
                                 ]
                            ),
                    html.Div(className = "item", 
                        style = {
                                 'grid-column-start': '1',
                                 'grid-column-end': '2',
                                 'grid-row-start': '2',
                                 'grid-row-end': '3'},
                        children=[
                                  dcc.Checklist(id = 'OrbCheck', 
                                                style = {
                                                         'textAlign' : 'left',
                                                         'fontFamily' : 'sans-serif'
                                                        },
                                                options=[
                                                         {'label': '1s', 'value': '1s'},
                                                         {'label': '2s', 'value': '2s'},
                                                         {'label': '3s', 'value': '3s'},
                                                         {'label': '4s', 'value': '4s'},
                                                         {'label': '5s', 'value': '5s'},
                                                         {'label': '6s', 'value': '6s'},
                                                         {'label': '2p', 'value': '2p'},
                                                         {'label': '3p', 'value': '3p'},
                                                         {'label': '4p', 'value': '4p'},
                                                         {'label': '5p', 'value': '5p'},
                                                         {'label': '6p', 'value': '6p'},
                                                         {'label': '3d', 'value': '3d'},
                                                         {'label': '4d', 'value': '4d'},
                                                         {'label': '5d', 'value': '5d'},
                                                         {'label': '6d', 'value': '6d'},
                                                         {'label': '4f', 'value': '4f'},
                                                         {'label': '5f', 'value': '5f'},
                                                         {'label': '6f', 'value': '6f'},
                                                        ],
                                                 value=['2p', '3p', '4p'],
                                                labelStyle={
                                                            'maxwidth' : '20px',
                                                            'display': 'inline-block'
                                                           }
                                               ),
                                 ]
                            ),
                    html.Div(className = "item", 
                        style = {
                                 'grid-column-start': '2',
                                 'grid-column-end': '3',
                                 'grid-row-start': '1',
                                 'grid-row-end': '2'},
                        children=[
                                  html.H2(style = {
                                                  'textAlign' : 'center', 
                                                  'fontFamily' : 'sans-serif'
                                                 },
                                         children = 'Function'),
                                 ]
                            ),
                    html.Div(className = "item", 
                        style = {
                                 'grid-column-start': '2',
                                 'grid-column-end': '3',
                                 'grid-row-start': '2',
                                 'grid-row-end': '3'},
                        children=[
                                  dcc.RadioItems(id = 'FuncType', 
                                                 style = {
                                                          'textAlign' : 'center',
                                                          'fontFamily' : 'sans-serif',
                                                          'display': 'block'
                                                         },
                                                 options=[
                                                          {
                                                           'label': 'Radial Distribution Function',
                                                           'value': 'RDF'
                                                          },
                                                          {
                                                           'label': 'Radial Wave Function',
                                                           'value': 'RWF'
                                                          },
                                                          {
                                                           'label': '3D', 
                                                           'value': '3DWF'
                                                          }
                                                         ],
                                                 value='RDF',
                                                 labelStyle={
                                                             'float':'left'
                                                             }
                                                 )
                                 ]
                            )
                  ]
        )]


plot_options_2d = [
                
            html.H2(style = {
                             'textAlign' : 'center', 
                             'fontFamily' : 'sans-serif'
                            },
                    children = 'Plot Options'),

            html.Div(className = "container", 
                     style = {'display' : 'grid',
                              'grid-template-columns': r'33% 33% 33%',
                              'grid-template-rows' : r'25% 25% 25% 25%'},
                     children=[
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '1',
                                             'grid-column-end': '2',
                                             'grid-row-start': '1',
                                             'grid-row-end': '2'},
                                    children=[
                                              html.P(style = {
                                                              'textAlign' : 'center', 
                                                              'fontFamily' : 'sans-serif'
                                                             },
                                                     children = 'Gridlines'),
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '2',
                                             'grid-column-end': '3',
                                             'grid-row-start': '1',
                                             'grid-row-end': '2'},
                                    children=[
                                              html.P(style = {
                                                              'textAlign' : 'center', 
                                                              'fontFamily' : 'sans-serif'
                                                             },
                                                     children = 'Lower x limit'),
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '3',
                                             'grid-column-end': '4',
                                             'grid-row-start': '1',
                                             'grid-row-end': '2'},
                                    children=[
                                              html.P(style = {
                                                              'textAlign' : 'center', 
                                                              'fontFamily' : 'sans-serif'
                                                             },
                                                     children = 'Upper x limit'),
                                             ]
                                        ),
                                html.Div(className = "container", 
                                    style = {'display' : 'grid',
                                             'grid-template-columns': r'100%',
                                             'grid-template-rows' : r'50% 50%'
                                            },
                                    children=[
                                              html.Div(className = "item", 
                                                       style = {
                                                                'grid-column-start': '1',
                                                                'grid-column-end': '2',
                                                                'grid-row-start': '1',
                                                                'grid-row-end': '2'
                                                               },
                                                       children=[
                                                                 dcc.Checklist(id = 'xgridcheck',
                                                                               style = {
                                                                                        'textAlign' : 'center', 
                                                                                        'fontFamily' : 'sans-serif'
                                                                                        },
                                                                               options=[
                                                                                        {
                                                                                         'label': ' x-Axis ', 
                                                                                         'value': 'xgridval'
                                                                                        }
                                                                               ],
                                                                               value=[],
                                                                               )
                                                                ]
                                                           ),
                                              html.Div(className = "item", 
                                                       style = {
                                                                'grid-column-start': '1',
                                                                'grid-column-end': '2',
                                                                'grid-row-start': '2',
                                                                'grid-row-end': '3'},
                                                       children=[
                                                                 dcc.Checklist(id = 'ygridcheck',
                                                                               style = {
                                                                                        'textAlign' : 'center', 
                                                                                        'fontFamily' : 'sans-serif'
                                                                                        },
                                                                               options=[
                                                                                        {
                                                                                         'label': ' y-Axis ', 
                                                                                         'value': 'ygridval'
                                                                                        }
                                                                               ],
                                                                               value=[],
                                                                               )
                                                                ]
                                                           )
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '2',
                                             'grid-column-end': '3',
                                             'grid-row-start': '2',
                                             'grid-row-end': '3'},
                                    children=[
                                              dcc.Input(
                                                  id = 'LowerxInput',
                                                  placeholder = 0,
                                                  type = 'number',
                                                  min = -10,
                                                  max = 100,
                                                  value = 0,
                                                  style = {'width':'40%'}
                                              ),
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '3',
                                             'grid-column-end': '4',
                                             'grid-row-start': '2',
                                             'grid-row-end': '3'},
                                    children=[
                                              dcc.Input(
                                                  id = 'UpperxInput',
                                                  placeholder = 80,
                                                  type='number',
                                                  max = 100,
                                                  value = 80,
                                                  style = {'width':'40%'}
                                              ),
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '1',
                                             'grid-column-end': '2',
                                             'grid-row-start': '3',
                                             'grid-row-end': '4'},
                                    children=[
                                              html.P(style = {
                                                              'textAlign' : 'center', 
                                                              'fontFamily' : 'sans-serif'
                                                             },
                                                     children = 'Line Width')
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '1',
                                             'grid-column-end': '2',
                                             'grid-row-start': '4',
                                             'grid-row-end': '5'},
                                    children=[
                                              dcc.Slider(
                                                         id = 'LineWidthSlider',
                                                         min=1,
                                                         max=10,
                                                         step=0.5,
                                                         value=5,
                                                        )
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '2',
                                             'grid-column-end': '3',
                                             'grid-row-start': '3',
                                             'grid-row-end': '4'},
                                    children=[
                                              html.P(style = {
                                                              'textAlign' : 'center', 
                                                              'fontFamily' : 'sans-serif'
                                                             },
                                                     children = 'Text Size')
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '2',
                                             'grid-column-end': '3',
                                             'grid-row-start': '4',
                                             'grid-row-end': '5'},
                                    children=[
                                              dcc.Slider(
                                                  id = 'TextSizeSlider',
                                                  min=15,
                                                  max=25,
                                                  step=0.5,
                                                  value=19,
                                                        )
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '3',
                                             'grid-column-end': '4',
                                             'grid-row-start': '3',
                                             'grid-row-end': '4'},
                                    children=[
                                              html.P(style = {
                                                              'textAlign' : 'center', 
                                                              'fontFamily' : 'sans-serif'
                                                             },
                                                     children = 'Plot Colours')
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '3',
                                             'grid-column-end': '4',
                                             'grid-row-start': '4',
                                             'grid-row-end': '5',
                                             'justify-self':'center',
                                             'align-self':'center',
                                             'padding-right' : '40%'
                                             },
                                    children=[
                                              dcc.Dropdown(
                                                             id = 'Colours_2d',
                                                             options=[
                                                                     { 
                                                                      'label': 'normal', 
                                                                      'value': 'normal'
                                                                     },
                                                                     {
                                                                      'label': 'deuteranopia',
                                                                      'value': 'deut'
                                                                     },
                                                                     {
                                                                      'label': 'protanopia', 
                                                                      'value': 'prot'
                                                                     }
                                                                     ],
                                                             style = {
                                                                      'width' : '150%'
                                                                     }, 
                                                             value='normal',
                                                             searchable=False,
                                                             clearable=False
                                                            )
                                             ]
                                        )

                              ]
                              )]



plot_options_3d = [
                
            html.H2(style = {
                             'textAlign' : 'center', 
                             'fontFamily' : 'sans-serif'
                            },
                    children = 'Plot Options'),

            html.Div(className = "container", 
                     style = {'display' : 'grid',
                              'grid-template-columns': r'50% 50%',
                              'grid-template-rows' : r'50% 50%',
                              'justify-items' : 'center',
                              'align-items' : 'center'
                             },
                     children=[
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '1',
                                             'grid-column-end': '2',
                                             'grid-row-start': '1',
                                             'grid-row-end': '2'},
                                    children=[
                                              html.P(style = {
                                                              'textAlign' : 'center', 
                                                              'fontFamily' : 'sans-serif'
                                                             },
                                                     children = 'Lobe Colours')
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '2',
                                             'grid-column-end': '3',
                                             'grid-row-start': '1',
                                             'grid-row-end': '2'},
                                    children=[
                                              html.P(style = {
                                                              'textAlign' : 'center', 
                                                              'fontFamily' : 'sans-serif'
                                                             },
                                                     children = 'Cutaway')
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '1',
                                             'grid-column-end': '2',
                                             'grid-row-start': '2',
                                             'grid-row-end': '3',
                                             'justify-self': 'stretch'
                                             },
                                    children=[
                                              dcc.Dropdown(
                                                             id = 'Colours_3d',
                                                             options=[
                                                                     { 
                                                                      'label': 'yellow/purple', 
                                                                      'value': 'yp'
                                                                     },
                                                                     {
                                                                      'label': 'red/blue',
                                                                      'value': 'rb'
                                                                     },
                                                                     {
                                                                      'label': 'green/orange', 
                                                                      'value': 'go'
                                                                     }
                                                                     ],
                                                             style = {
                                                                     }, 
                                                             value='yp',
                                                             searchable=False,
                                                             clearable=False
                                                            )
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '2',
                                             'grid-column-end': '3',
                                             'grid-row-start': '2',
                                             'grid-row-end': '3'},
                                    children=[
                                              dcc.RadioItems(id = 'Cutaway', 
                                                             style = {
                                                                      'textAlign' : 'center',
                                                                      'fontFamily' : 'sans-serif',
                                                                      'display': 'block'
                                                                     },
                                                             options=[
                                                                      {
                                                                       'label': '0', 
                                                                       'value': 1.0
                                                                      },
                                                                      # {
                                                                      #  'label': '1/4',
                                                                      #  'value': 0.5
                                                                      # },
                                                                      {
                                                                       'label': '1/2',
                                                                       'value': 0.
                                                                      },
                                                                     ],
                                                             value=1.0,
                                                             labelStyle={
                                                                         'float':'left'
                                                                         }
                                                             )
                                              ]
                                        )

                              ]
                    )]

plot_save_options = [

            html.H2(style = {
                             'textAlign' : 'center', 
                             'fontFamily' : 'sans-serif'
                            },
                    children = 'Save Options'
                    ),
                
            html.Div(className = "container", 
                     style = {
                              'display' : 'grid',
                              'grid-template-columns': r'33% 33% 33%',
                              'grid-template-rows' : r'50% 50%',
                               'justify-items' : 'center',
                               'align-items' : 'center'
                              },
                     children=[
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '1',
                                             'grid-column-end': '2',
                                             'grid-row-start': '1',
                                             'grid-row-end': '2'
                                             },
                                    children=[
                                              html.P(style = {
                                                              'textAlign' : 'center', 
                                                              'fontFamily' : 'sans-serif'
                                                             },
                                                     children = 'Output Height'),
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '2',
                                             'grid-column-end': '3',
                                             'grid-row-start': '1',
                                             'grid-row-end': '2'
                                             },
                                    children=[
                                              html.P(style = {
                                                              'textAlign' : 'center', 
                                                              'fontFamily' : 'sans-serif'
                                                             },
                                                     children = 'Output Width'),
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '1',
                                             'grid-column-end': '2',
                                             'grid-row-start': '2',
                                             'grid-row-end': '3',
                                             'padding-left': '40%',
                                             },
                                    children=[
                                              dcc.Input(
                                                  id = 'PlotHeightInput',
                                                  placeholder=500,
                                                  type='number',
                                                  value=500,
                                                  style = {
                                                           'width':'35%'
                                                           }
                                                        )
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '2',
                                             'grid-column-end': '3',
                                             'grid-row-start': '2',
                                             'grid-row-end': '3',
                                             'padding-left': '40%',
                                             },
                                    children=[
                                              dcc.Input(
                                                  id = 'PlotWidthInput',
                                                  placeholder=700,
                                                  type='number',
                                                  value=700,
                                                  style = {
                                                           'width':'35%'
                                                           }
                                                        )
                                             ]
                                        ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '3',
                                             'grid-column-end': '4',
                                             'grid-row-start': '2',
                                             'grid-row-end': '3',
                                             'justify-self':'center'
                                             },
                                    children=[                                 
                                                dcc.Dropdown(
                                                             id = 'PlotFormatDropdown',
                                                             options=[
                                                                     { 
                                                                      'label': 'svg', 
                                                                      'value': 'svg'
                                                                     },
                                                                     {
                                                                      'label': 'png',
                                                                      'value': 'png'
                                                                     },
                                                                     {
                                                                      'label': 'jpeg', 
                                                                      'value': 'jpeg'
                                                                     }
                                                                     ],
                                                             style = {
                                                                      'width' : '130%'
                                                                     }, 
                                                             value='png',
                                                             searchable=False,
                                                             clearable=False
                                                            )
                                              ]
                                         ),
                                html.Div(className = "item", 
                                    style = {
                                             'grid-column-start': '3',
                                             'grid-column-end': '4',
                                             'grid-row-start': '1',
                                             'grid-row-end': '2'
                                            },
                                    children=[                                 
                                              html.P(style = {
                                                              'textAlign' : 'center', 
                                                              'fontFamily' : 'sans-serif'
                                                             },
                                                     children = 'Output Format'),
                                              ]
                                         )

                              ]
                              )]

# Layout of webpage
app.layout = html.Div(children=[

html.Title(children='orbiplot'),

html.Div(style = {'background-color': '#3977AF'},
            children = [

                        html.H1(style = {
                                         'textAlign'  : 'center',
                                         'fontFamily' : 'sans-serif',
                                         'color'      : 'white'
                                         },
                                children = 'Orbiplot')
                        ]
            ),

html.Div(className = "container", 
         style = { 'display' : 'grid',
                   'grid-template-columns': r'50% 50%',
                   'grid-template-rows' : r'100%',
                   'margin': 'auto',
                   'width' : '95%'

                  },
         children=[
                    html.Div(className = "item", 
                        style = {
                                 'grid-column-start': '1',
                                 'grid-column-end': '2',
                                 'grid-row-start': '1',
                                 'grid-row-end': '2',
                                 'justify-self': 'stretch',
                                 'align-self': 'stretch'
                                 },
                        children=[
                                  dcc.Graph(id='RDF_Graph', 
                                            style = {
                                                     'responsive' : 'true',
                                                     'height' : '580px',
                                                     'automargin' : 'true'
                                                     }
                                            )
                                 ]
                            ),
                    html.Div(className = "item", 
                        style = {
                                 'grid-column-start': '2',
                                 'grid-column-end': '3',
                                 'grid-row-start': '1',
                                 'grid-row-end': '2'
                                 },
                        children=[
                                  html.Div(className = "container", 
                                           style = { 'display' : 'grid',
                                                     'grid-template-columns': r'100%',
                                                     'margin': 'auto',
                                                     'justify-content':'stretch',
                                                     'align-self':'center'
                                                    },
                                           children=[
                                                     html.Div(className = "item", 
                                                         style = {
                                                                  'grid-column-start': '1',
                                                                  'grid-column-end': '2',
                                                                  'grid-row-start': '1',
                                                                  'grid-row-end': '2',
                                                                  },
                                                         children=[
                                                                    html.Div(id        = 'orbitals box', 
                                                                             style     = {}, 
                                                                             children  = orbital_plot_options
                                                                             ),
                                                                  ]
                                                             ),
                                                     html.Div(className = "item", 
                                                         style = {
                                                                  'grid-column-start': '1',
                                                                  'grid-column-end': '2',
                                                                  'grid-row-start': '2',
                                                                  'grid-row-end': '3',
                                                                  },
                                                         children=[
                                                                    html.Div(id        = 'plot_options_box_2d', 
                                                                             style     = {}, 
                                                                             children  = plot_options_2d
                                                                             ),
                                                                  ]
                                                             ),
                                                     html.Div(className = "item", 
                                                         style = {
                                                                  'grid-column-start': '1',
                                                                  'grid-column-end': '2',
                                                                  'grid-row-start': '2',
                                                                  'grid-row-end': '3',
                                                                  },
                                                         children=[
                                                                    html.Div(id        = 'plot_options_box_3d', 
                                                                             style     = {}, 
                                                                             children  = plot_options_3d
                                                                             ),
                                                                  ]
                                                             ),
                                                     html.Div(className = "item", 
                                                         style = {
                                                                  'grid-column-start': '1',
                                                                  'grid-column-end': '2',
                                                                  'grid-row-start': '3',
                                                                  'grid-row-end': '4',
                                                                  },
                                                         children=[
                                                                    html.Div(className = 'Save Options Box', 
                                                                             style = {}, 
                                                                             children  = plot_save_options
                                                                            )

                                                                  ]
                                                             ),
                                                     ]
                                          )

                                 ]
                            )
                 ]
        ),
html.Footer(style = {
                     'textAlign'       : 'center', 
                     'fontFamily'      : 'sans-serif',
                     'font-size'       : 'smaller',
                     'color'           : 'white',
                     'background-color': '#3977AF'
                    }, 
            children=[
                      html.P(children = ['Author: Jon Kragskow']),
                      html.A(href = 'www.kragskow.com',
                             style = {'color':'white'},
                             children = 'www.kragskow.com')
                     ]
           ),
refresh_plots,
mathjax_script
])

def toggle_pob(on_off) :

    if on_off == 'on' :
        return   {'display' : 'inline-block' , 
                'padding-left':'5%', 
                'padding-bottom':'5%' , 
                'padding-right':'5%', 
                'width' : '90%', 
                'borderStyle' : 'solid', 
                'borderWidth' : '0px',
                'textAlign' : 'center'}
    else :
        return {'display' : 'none' , 
                'padding-left':'5%', 
                'padding-bottom':'5%' , 
                'padding-right':'5%', 
                'width' : '90%', 
                'borderStyle' : 'solid', 
                'borderWidth' : '0px',
                'textAlign' : 'center'}

def orb_checklist(dimension):
    if dimension == '2d':
        return[
               {'label': '1s', 'value': '1s'},
               {'label': '2s', 'value': '2s'},
               {'label': '3s', 'value': '3s'},
               {'label': '4s', 'value': '4s'},
               {'label': '5s', 'value': '5s'},
               {'label': '6s', 'value': '6s'},
               {'label': '2p', 'value': '2p'},
               {'label': '3p', 'value': '3p'},
               {'label': '4p', 'value': '4p'},
               {'label': '5p', 'value': '5p'},
               {'label': '6p', 'value': '6p'},
               {'label': '3d', 'value': '3d'},
               {'label': '4d', 'value': '4d'},
               {'label': '5d', 'value': '5d'},
               {'label': '6d', 'value': '6d'},
               {'label': '4f', 'value': '4f'},
               {'label': '5f', 'value': '5f'},
               {'label': '6f', 'value': '6f'},
              ]

    elif dimension == '3d':
        return[
               {'label': '2p', 'value': '2p'},
               {'label': '3p', 'value': '3p'},
               {'label': '4p', 'value': '4p'},
               {'label': '5p', 'value': '5p'},
               {'label': '6p', 'value': '6p'},
               {'label': '3dz²', 'value': '3d_z2'},
               {'label': '4dz²', 'value': '4d_z2'},
               {'label': '5dz²', 'value': '5d_z2'},
               {'label': '6dz²', 'value': '6d_z2'},
               {'label': '3dxy', 'value': '3dxy'},
               {'label': '4dxy', 'value': '4dxy'},
               {'label': '5dxy', 'value': '5dxy'},
               {'label': '6dxy', 'value': '6dxy'},
              ]
##################################################################################################################################
########################################################### Callbacks ############################################################
##################################################################################################################################


#Callback which defines what changes ('figure') and what causes the change (Checkbox being pressed)
@app.callback([ddep.Output('RDF_Graph', 'figure'),
               ddep.Output('RDF_Graph', 'config'),
               ddep.Output('plot_options_box_2d', 'style'),
               ddep.Output('plot_options_box_3d', 'style'),
               ddep.Output('OrbCheck', 'options')],
              [ddep.Input('OrbCheck', 'value'),
               ddep.Input('FuncType', 'value'), 
               ddep.Input('LineWidthSlider','value'), 
               ddep.Input('TextSizeSlider','value'), 
               ddep.Input('xgridcheck','value'), 
               ddep.Input('ygridcheck','value'), 
               ddep.Input('UpperxInput', 'value'),
               ddep.Input('LowerxInput', 'value'),
               ddep.Input('PlotFormatDropdown', 'value'), 
               ddep.Input('PlotHeightInput', 'value'),
               ddep.Input('PlotWidthInput', 'value'),
               ddep.Input('Colours_2d','value'),
               ddep.Input('Colours_3d','value'),
               ddep.Input('Cutaway','value')])

#Function which is called after the element is pressed
def UpdatePlot(Orbitals, FuncType, Thickness, TextSize, xgridinput,
               ygridinput, upperxlim, lowerxlim,PlotFormat,
               PlotHeight, PlotWidth,colour_name_2d, colour_name_3d, cutaway):

    # If xlims nonetype then set to default
    if lowerxlim == None:
        lowerxlim = 0

    if upperxlim == None:
        upperxlim = lowerxlim + 80

    #Read type of wavefunction being requested and set axis names
    if FuncType == 'RDF':
        WFFlag = 1
        WFName = r'$\text{Radial Distribution Function}  \ \ \ 4\pi r^2 R(r)^2$'
        file_name = 'radial_distribution_functions'
    if FuncType == 'RWF':
        WFFlag = 2
        WFName = r'$\text{Radial Wavefunction}  \ \ \ R(r) $\n'
        file_name = 'radial_wavefunctions'
    if FuncType == '3DWF':
        WFFlag = 3
        WFName = ''

    # Plot radial wavefunction or radial distribution function
    if WFFlag == 2 or WFFlag == 1:
        return {
                'data': get_2d_plots(Orbitals, max(lowerxlim,0), max(upperxlim,0), WFFlag, Thickness),
                'layout': ax_config_2d(xgridinput, ygridinput, TextSize, WFName, lowerxlim, upperxlim)
                }, modebar_config(PlotFormat, PlotHeight, PlotWidth, file_name), toggle_pob('on'), toggle_pob('off'), orb_checklist('2d')

    # Plot wavefunction isosurface
    # use first selection
    if WFFlag == 3 :
        if len(Orbitals) > 0:

            fig = make_subplots(rows=1, cols=1,
                                specs=[[{'is_3d': True}]]
                               )

            fig.update_layout(ax_config_3d())

            fig = OrbCalc(Orbitals[0], colour_name_3d, fig, np.float64(cutaway))

            return  fig,modebar_config(PlotFormat, PlotHeight, PlotWidth, Orbitals[0]+' Orbital'), toggle_pob('off'), toggle_pob('on'), orb_checklist('3d')
        else:
            return {
                    'data'   : None,
                    'layout' : ax_config_3d()
                    }, modebar_config(PlotFormat, PlotHeight, PlotWidth, ' Orbital'), toggle_pob('off'), toggle_pob('on'), orb_checklist('3d')


def get_2d_plots(Orbitals, lowerxlim, upperxlim, WFFlag, Thickness):

    # Calculate plot for each orbital and add to list
    traces = []
    curr_ymax = 0.
    if len(Orbitals) > 0:
        for n in range(0,len(Orbitals)):
            traces.append(go.Scatter(x = np.linspace(lowerxlim,upperxlim,1000), 
                                     y = functiondict[Orbitals[n]](np.linspace(lowerxlim,upperxlim,1000), WFFlag, 1.),
                                     line = dict(width = Thickness),
                                     name = Orbitals[n], 
                                     hoverinfo = 'none')
                                    )
    return traces


def ax_config_3d():
    return go.Layout(
                     hovermode=False,
                     dragmode="orbit",
                     scene_aspectmode='cube',
                     scene=dict(
                                xaxis=dict(
                                           gridcolor='rgb(255, 255, 255)',
                                           zerolinecolor='rgb(255, 255, 255)',
                                           showbackground=False,
                                           showgrid=False,
                                           zeroline=False,
                                           title='',
                                           showline=False,
                                           ticks='',
                                           showticklabels=False,
                                           backgroundcolor='rgb(255, 255,255)'
                                          ),
                                yaxis=dict(
                                           gridcolor='rgb(255, 255, 255)',
                                           zerolinecolor='rgb(255, 255, 255)',
                                           showgrid=False,
                                           zeroline=False,
                                           title='',
                                           showline=False,
                                           ticks='',
                                           showticklabels=False,
                                           backgroundcolor='rgb(255, 255,255)'
                                          ),
                                zaxis=dict(
                                           gridcolor='rgb(255, 255, 255)',
                                           zerolinecolor='rgb(255, 255, 255)',
                                           showgrid=False,
                                           title='',
                                           zeroline=False,
                                           showline=False,
                                           ticks='',
                                           showticklabels=False,
                                           backgroundcolor='rgb(255, 255,255)'
                                          ),
                                aspectratio=dict(
                                                 x=1,
                                                 y=1,
                                                 z=1
                                                ),
                                
                                ),
                     margin=dict(l=20, r=30, t=30, b=20),
                    )

def ax_config_2d(xgridinput, ygridinput, TextSize, WFName, xlow, xup):

    #Turn on x axis gridlines
    if xgridinput == ['xgridval']:
        xgrid = True
    else:
        xgrid = False

    #Turn on y axis gridlines
    if ygridinput == ['ygridval']:
        ygrid = True
    else:
        ygrid = False

    return go.Layout(
                     xaxis = {
                              'autorange' : True,
                              'showgrid'  : xgrid,
                              'zeroline'  : False,
                              'showline'  : True,
                              'range'     : [xlow, xup],
                              'title' : {
                                         'text' : "$(Distance a_0)$",
                                         'font' :
                                                 {
                                                  'size' : TextSize
                                                 } 
                                        },
                              'ticks' :'outside',
                              'tickfont' : 
                                           {
                                            'size' : TextSize
                                           },
                              'showticklabels' : True
                             },
                     yaxis = {
                              'autorange' :True,
                              'showgrid' :ygrid,
                              'zeroline' :False,
                              'fixedrange' : True,
                              'title' : 
                                        {
                                         'text' : WFName,
                                         'font' : 
                                                  {
                                                   'size' : TextSize
                                                  }
                                        },
                              'title_standoff' : 100,
                              'showline' :True,
                              'ticks' :'outside',
                              'tickfont' : 
                                          {
                                           'size' : TextSize
                                          },
                              'showticklabels' :True
                             },
                     legend = {
                               'x' : 0.8,
                               'y' : 1,
                               'font' : 
                                       {
                                        'size' : TextSize - 3
                                       }
                              },
                     margin=dict(l=90, r=30, t=30, b=60),
    )

def modebar_config(PlotFormat, PlotHeight, PlotWidth, file_name):

    return  {
            'toImageButtonOptions': 
              {
              'format': PlotFormat, 
              'filename': file_name,
              'height': PlotHeight,
              'width': PlotWidth,
              }, 
            'modeBarButtonsToRemove':
                [
                'sendDataToCloud',
                'autoScale2d',
                'resetScale2d',
                'hoverClosestCartesian',
                'toggleSpikelines',
                'zoom2d',
                'zoom3d',
                'pan3d',
                'pan2d',
                'select2d',
                'zoomIn2d',
                'zoomOut2d',
                'hovermode',
                'resetCameraLastSave3d', 
                'hoverClosest3d',
                'hoverCompareCartesian',
                'resetViewMapbox',
                'orbitRotation', 
                'tableRotation',
                'resetCameraDefault3d'
                ],
            'displaylogo': False,
            'displayModeBar' : True,
            }



if __name__ == '__main__':

    app.run_server(debug=True, port=8053)
