#! /usr/bin/env python3

import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash.dependencies as ddep
import dash_defer_js_import as dji
import plotly.colors as pc
import plotly.graph_objs as go
import numpy as np
import os
import time

from app import app

def p_z_ax(n,r,f,c):
    return n*np.exp(r/n)*c*np.sqrt(np.pi/3.)/abs(f)

############################ Radial functions #########################
def radial_s(n, rho):
    """
    Calculates radial Wavefunction of s orbital
    for the specified principal quantum number

    Input:
        n (int)  ::  principal quantum number
        rho (numpy float array) :: 2.*r/n, where r^2 = x^2+y^2+z^2

    Returns:
        rad (numpy float array) :: radial wavefunction
    """
    if n == 1:
        rad = 2.*np.exp(-rho/2.)
    if n == 2:
        rad = 1./(2.*np.sqrt(2.))*(2.-rho)*np.exp(-rho/2.)
    if n == 3:
        rad = 1./(9.*np.sqrt(3.))*(6.-6.*rho+rho**2.)*np.exp(-rho/2.)
    if n == 4:
        rad = (1./96.)*(24.-36.*rho+12.*rho**2.-rho**3.)*np.exp(-rho/2.)
    if n == 5:
        rad = (1./(300.*np.sqrt(5.)))*(120.-240.*rho+120.*rho**2.-20.*rho**3.+rho**4.)*np.exp(-rho/2.)
    if n == 6:
        rad = (1./(2160.*np.sqrt(6.)))*(720.-1800.*rho+1200.*rho**2.-300.*rho**3.+30.*rho**4.-rho**5.)*np.exp(-rho/2.)

    return rad

def radial_p(n, rho):
    """
    Calculates radial Wavefunction of p orbital
    for the specified principal quantum number

    Input:
        n (int)  ::  principal quantum number
        rho (numpy float array) :: 2.*r/n, where r^2 = x^2+y^2+z^2

    Returns:
        rad (numpy float array) :: radial wavefunction
    """

    if n == 2:
        rad = 1./(2.*np.sqrt(6.))*rho*np.exp(-rho/2.)
    elif n == 3:
        rad = 1./(9.*np.sqrt(6.))*rho*(4.-rho)*np.exp(-rho/2.)
    elif n == 4:
        rad = 1./(32.*np.sqrt(15.))*rho*(20.-10.*rho+rho**2.)*np.exp(-rho/2.)
    elif n == 5:
        rad = 1./(150.*np.sqrt(30.))*rho*(120.-90.*rho+18.*rho**2.-rho**3.)*np.exp(-rho/2.)
    elif n == 6:
        rad = 1./(432.*np.sqrt(210.))*rho*(840.-840.*rho+252.*rho**2.-28.*rho**3.+rho**4.)*np.exp(-rho/2.)
    return rad

def radial_p_mod(n, r):
    """
    Calculates radial Wavefunction of p orbital
    for the specified principal quantum number

    Input:
        n (int)  ::  principal quantum number
        r (numpy float array) :: r^2 = x^2+y^2+z^2

    Returns:
        rad (numpy float array) :: radial wavefunction
    """

    if n == 2:
        rad = 1./(2.*np.sqrt(6.))
    elif n == 3:
        rad = (4. - 2.*r/3.)/(9*np.sqrt(6.))
    elif n == 4:
        rad = (20. - 5.*r + (r/2.)**2.)/(32.*np.sqrt(15.))
    elif n == 5:
        rad = (120. - 90.*(2*r/5.) + 18.*(2*r/5.)**2. - (2.*r/5.)**3.)/(150.*np.sqrt(30.))
    elif n == 6:
        rad = (840. - 840.*r/3. + 252.*(r/3.)**2. - 28.*(r/3.)**3. + (r/3.)**4.)/(432.*np.sqrt(210.))
    return rad

def radial_d(n,rho): 
    """
    Calculates radial Wavefunction of d orbital
    for the specified principal quantum number

    Input:
        n (int)  ::  principal quantum number
        rho (numpy float array) :: 2.*r/n, where r^2 = x^2+y^2+z^2

    Returns:
        rad (numpy float array) :: radial wavefunction
    """

    if n == 3:
        rad = 1./(9.*np.sqrt(30.))*rho**2.*np.exp(-rho/2.)
    elif n == 4:
        rad = 1./(96.*np.sqrt(5.))*(6.-rho)*rho**2.*np.exp(-rho/2.)
    elif n == 5:
        rad = 1./(150.*np.sqrt(70.))*(42.-14.*rho+rho**2)*rho**2.*np.exp(-rho/2.)
    elif n == 6:
        rad = 1./(864.*np.sqrt(105.))*(336.-168.*rho+24.*rho**2.-rho**3.)*rho**2.*np.exp(-rho/2.)

    return rad

def radial_f(n,rho): 
    """
    Calculates radial Wavefunction of f orbital
    for the specified principal quantum number

    Input:
        n (int)  ::  principal quantum number
        rho (numpy float array) :: 2.*r/n, where r^2 = x^2+y^2+z^2

    Returns:
        rad (numpy float array) :: radial wavefunction
    """     
    if n == 4:
        rad = 1./(96.*np.sqrt(35.))*rho**3.*np.exp(-rho/2.)
    elif n == 5:
        rad = 1./(300.*np.sqrt(70.))*(8.-rho)*rho**3.*np.exp(-rho/2.)
    elif n == 6:
        rad = 1./(2592.*np.sqrt(35.))*(rho**2.-18.*rho+72.)*rho**3.*np.exp(-rho/2.)

    return rad

####################################################################################################
########################################## Lobe bounds #############################################
####################################################################################################

def p_domain(n):
    # !!! Returns an array containing the bounds of r for each lobe of the requested p orbital

    #####2p orbital#####
    if n == 2:
        num_lobes = 2
        lobe_1_lower = 0.030542415534
        lobe_1_upper = 11.973154175074
       
        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper]])

    #####3p orbital#####
    if n == 3:
        num_lobes = 4
        lobe_1_lower = 0.0521008719026 
        lobe_1_upper =  5.6457958052515
        
        lobe_2_lower = 6.4019092339593
        lobe_2_upper = 20.7317063563592

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper]])
    
    #####4p orbital#####

    if n == 4:
        num_lobes = 6

        lobe_1_lower = 0.0791782710855
        lobe_1_upper = 5.0739651680262
        
        lobe_2_lower = 6.0725425741786
        lobe_2_upper = 12.8640838577464
        
        lobe_3_lower = 16.5661209026166
        lobe_3_upper = 30.3121031661214

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper],[ lobe_3_lower, lobe_3_upper]])
    
    #####5p orbital#####

    if n == 5:
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
    
    if n == 6:
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

        colours = [[0, 'rgb(0,90,181)'], [1, 'rgb(220,50,32)']]

    return colours

def name_to_qn(orb):
    """
    Separates string containing orbital name into n and l quantum numbers

    Input:
        orb (string) :: name of orbital
    Returns:
        n (integer)  :: principal quantum number
        l (integer)  :: azimuthal quantum number
    """
    n = int(orb[0])
    l = orb[1]

    return n,l

def calc_s_orb(n, cutaway):
    """
    Calculates s orbital wavefunction on a grid

    Input:
        n (int)       :: prinipal quantum number of orbital
        cutaway (int) :: number used to split orbital in half 
    Returns:
        x (np.mgrid)   :: x values
        y (np.mgrid)   :: y values
        z (np.mgrid)   :: z values
        wav (np.mgrid) :: wavefunction values at x, y, z
        upper (float)  :: max value of axes
        lower (float)  :: min value of axes
        ival, (float)  :: isoval to plot orbital at
    """


    if n == 1:
        upper =  10
        lower = -10
        step  =  50j
    if n == 2:
        upper =  17
        lower = -17
        step  =  50j
    if n == 3:
        upper =  30
        lower = -30
        step  =  50j
    elif n == 4:
        upper =  45
        lower = -45
        step  =  50j
    elif n == 5:
        upper =  58
        lower = -58
        step  =  60j
    elif n ==6:
        upper =  75
        lower = -75
        step  =  70j

    x,y,z = np.mgrid[upper:lower:step, upper:lower*cutaway:step, upper:lower:step]

    r = np.sqrt(x**2 + y**2 + z**2)

    rho = 2*r/n

    rad = radial_s(n, 2*r/n)

    ang = 0.5/np.sqrt(np.pi)

    wav = ang*rad

    if n == 1:
        ival=0.0005
    if n == 2:
        ival=0.0005
    if n == 3:
        ival=0.0005
    elif n == 4:
        ival = 0.0005
    elif n == 5:
        ival =0.0005
    elif n == 6:
        ival =0.0005

    return x, y, z, wav, upper, lower, ival

def plot_p_orb(n, cutaway, colours):

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
    orb_r_bounds, num_lobes = p_domain(n)

    data = []

    flag = True

    #Calculate coordinates of isosurface
    # loop over each lobe of the orbital

    num_shells = n-1
    sw = 0

    for shell in np.arange(num_shells):

        xt, yt, zt = calc_p_orb(n, c, orb_r_bounds[shell, :], ang, num_lobes, angle_steps, 
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

        x  = np.append(x, xt, axis=0)
        y  = np.append(y, yt, axis=0)
        z  = np.append(z, -zt, axis=0)
        co = np.append(co, cot2, axis=0)

    data = go.Surface(x=x, y=y, z=z, surfacecolor = co, colorscale=colours, showscale=False)

    return data, -orb_r_bounds[-1,1], orb_r_bounds[-1,1]

def calc_p_orb(n, c, bounds, ang, num_lobes, angle_steps, r_steps, r_mini_steps):

    gap = np.abs(bounds[1] - bounds[0])

    #Array of r values for each lobe
    # Sample more frequently closer to the pole of the lobe
    r = np.linspace(bounds[0], bounds[0] + gap*0.015, r_mini_steps)
    r = np.append(r,np.linspace(bounds[0]+ gap*0.015, bounds[0] + gap*0.5, r_mini_steps))
    r = np.append(r,np.linspace(bounds[0]+ gap*0.5, bounds[0] + gap*0.975, r_mini_steps))
    r = np.append(r,np.linspace(bounds[0]+ gap*0.975, bounds[0] + gap, r_mini_steps))

    angm, rm = np.meshgrid(ang, r)

    zor = p_z_ax(n,rm,radial_p_mod(n, rm),c)

    x = np.sqrt(rm**2. - zor**2.)*np.cos(angm)
    y = np.sqrt(rm**2. - zor**2.)*np.sin(angm)
    z = zor

    return x, y, z

def calc_dz_orb(n, cutaway):
    """
    Calculates dz2 orbital wavefunction on a grid

    Input:
        n (int)       :: prinipal quantum number of orbital
        cutaway (int) :: number used to split orbital in half 
    Returns:
        x (np.mgrid)   :: x values
        y (np.mgrid)   :: y values
        z (np.mgrid)   :: z values
        wav (np.mgrid) :: wavefunction values at x, y, z
        upper (float)  :: max value of axes
        lower (float)  :: min value of axes
        ival, (float)  :: isoval to plot orbital at
    """
    if n == 3:
        upper = 40
        lower = -40
        step  = 60j
    elif n == 4:
        upper = 70
        lower = -70
        step  = 70j
    elif n == 5:
        upper = 98
        lower = -98
        step  = 80j
    elif n ==6:
        upper = 135
        lower = -135
        step  = 90j

    x,y,z = np.mgrid[upper:lower:step, upper:lower*cutaway:step, upper:lower:step]

    r = np.sqrt(x**2 + y**2 + z**2)

    rho = 2*r/n

    rad = radial_d(n, 2*r/n)

    ang = 2*z**2-x**2-y**2

    wav = rad*ang

    if n == 3:
        ival=0.08
    elif n == 4:
        ival = 0.08
    elif n == 5:
        ival = 0.08
    elif n == 6:
        ival = 0.08

    return x, y, z, wav, upper, lower, ival


def calc_dxy_orb(n, cutaway):
    """
    Calculates dxy orbital wavefunction on a grid

    Input:
        n (int)       :: prinipal quantum number of orbital
        cutaway (int) :: number used to split orbital in half 
    Returns:
        x (np.mgrid)   :: x values
        y (np.mgrid)   :: y values
        z (np.mgrid)   :: z values
        wav (np.mgrid) :: wavefunction values at x, y, z
        upper (float)  :: max value of axes
        lower (float)  :: min value of axes
        ival, (float)  :: isoval to plot orbital at
    """

    if n == 3:
        upper =  45
        lower = -45
        step  =  60j
    elif n == 4:
        upper =  70
        lower = -70
        step  =  70j
    elif n == 5:
        upper =  98
        lower = -98
        step  =  80j
    elif n ==6:
        upper =  135
        lower = -135
        step  =  90j

    x,y,z = np.mgrid[upper:lower:step, upper:lower:step, upper:lower*cutaway:step]

    r = np.sqrt(x**2 + y**2 + z**2)

    rad = radial_d(n, 2*r/n)

    ang = x*y

    wav = rad*ang

    if n == 3:
        ival=0.005
    elif n == 4:
        ival = 0.01
    elif n == 5:
        ival = 0.01
    elif n == 6:
        ival = 0.01
        
    return x, y, z, wav, upper, lower, ival


def calc_fz_orb(n, cutaway):
    """
    Calculates fz3 orbital wavefunction on a grid

    Input:
        n (int)       :: prinipal quantum number of orbital
        cutaway (int) :: number used to split orbital in half 
    Returns:
        x (np.mgrid)   :: x values
        y (np.mgrid)   :: y values
        z (np.mgrid)   :: z values
        wav (np.mgrid) :: wavefunction values at x, y, z
        upper (float)  :: max value of axes
        lower (float)  :: min value of axes
        ival, (float)  :: isoval to plot orbital at
    """

    if n == 4:
        upper = 70
        lower = -70
        step  = 60j
    elif n == 5:
        upper = 100
        lower = -100
        step  = 75j
    elif n ==6:
        upper = 130
        lower = -130
        step  = 85j

    x,y,z = np.mgrid[upper:lower*cutaway:step, upper:lower:step, upper:lower:step]

    r = np.sqrt(x**2 + y**2 + z**2)

    rad = radial_f(n, 2*r/n)

    ang = 0.25 * np.sqrt(7/np.pi) * z*(2*z**2-3*x**2-3*y**2)/(r**3)

    wav = rad*ang

    if n == 4:
        ival = 0.000005
    elif n == 5:
        ival = 0.000005
    elif n == 6:
        ival = 0.000005
  
    return x, y, z, wav, upper, lower, ival

def calc_fxyz_orb(n, cutaway):
    """
    Calculates fxyz orbital wavefunction on a grid

    Input:
        n (int)       :: prinipal quantum number of orbital
        cutaway (int) :: number used to split orbital in half 
    Returns:
        x (np.mgrid)   :: x values
        y (np.mgrid)   :: y values
        z (np.mgrid)   :: z values
        wav (np.mgrid) :: wavefunction values at x, y, z
        upper (float)  :: max value of axes
        lower (float)  :: min value of axes
        ival, (float)  :: isoval to plot orbital at
    """

    if n == 4:
        upper = 60
        lower = -60
        step  = 60j
    elif n == 5:
        upper = 90
        lower = -90
        step  = 70j
    elif n ==6:
        upper = 115
        lower = -115
        step  = 80j

    x,y,z = np.mgrid[upper:lower*cutaway:step, upper:lower:step, upper:lower:step]

    r = np.sqrt(x**2 + y**2 + z**2)

    rad = radial_f(n, 2*r/n)

    ang = 0.5* np.sqrt(105/np.pi) * x*y*z/(r**3)

    wav = rad*ang

    if n == 4:
        ival = 0.000005
    elif n == 5:
        ival = 0.000005
    elif n == 6:
        ival = 0.000005
  
    return x, y, z, wav, upper, lower, ival

def calc_fyz2_orb(n, cutaway):
    """
    Calculates fyz2 orbital wavefunction on a grid

    Input:
        n (int)       :: prinipal quantum number of orbital
        cutaway (int) :: number used to split orbital in half 
    Returns:
        x (np.mgrid)   :: x values
        y (np.mgrid)   :: y values
        z (np.mgrid)   :: z values
        wav (np.mgrid) :: wavefunction values at x, y, z
        upper (float)  :: max value of axes
        lower (float)  :: min value of axes
        ival, (float)  :: isoval to plot orbital at
    """

    if n == 4:
        upper = 65
        lower = -65
        step  = 60j
    elif n == 5:
        upper = 90
        lower = -90
        step  = 70j
    elif n ==6:
        upper = 120
        lower = -120
        step  = 80j

    x,y,z = np.mgrid[upper:lower*cutaway:step, upper:lower:step, upper:lower:step]

    r = np.sqrt(x**2 + y**2 + z**2)

    rad = radial_f(n, 2*r/n)

    ang = 0.25* np.sqrt(35/(2*np.pi)) * (3*x**2-y**2)*y/r**3

    wav = rad*ang

    if n == 4:
        ival = 0.000005
    elif n == 5:
        ival = 0.000005
    elif n == 6:
        ival = 0.000005
  
    return x, y, z, wav, upper, lower, ival

def swap(colour_1, colour_2):
    return colour_2, colour_1

def calc_radial_s(n, r, wf_type):
    """
    Calculates s orbital rwf or rdf

    Input:
        n (int)      :: prinipal quantum number of orbital
        r (np.array) :: r values
    Returns:
        y (np.array) :: y values corresponding to rwf or rdf
    """
    if "RDF" in wf_type:
        return r**2.* radial_s(n, 2.*r/n)**2
    if "RWF" in wf_type:
        return radial_s(n, 2.*r/n)

def calc_radial_p(n, r, wf_type):
    """
    Calculates p orbital rwf or rdf

    Input:
        n (int)      :: prinipal quantum number of orbital
        r (np.array) :: r values
    Returns:
        y (np.array) :: y values corresponding to rwf or rdf
    """
    if "RDF" in wf_type:
        return r**2.* radial_p(n, 2.*r/n)**2
    if "RWF" in wf_type:
        return radial_p(n, 2.*r/n)

def calc_radial_d(n, r, wf_type):
    """
    Calculates d orbital rwf or rdf

    Input:
        n (int)      :: prinipal quantum number of orbital
        r (np.array) :: r values
    Returns:
        y (np.array) :: y values corresponding to rwf or rdf
    """
    if "RDF" in wf_type:
        return r**2.* radial_d(n, 2.*r/n)**2
    if "RWF" in wf_type:
        return radial_d(n, 2.*r/n)

def calc_radial_f(n, r, wf_type):
    """
    Calculates f orbital rwf or rdf

    Input:
        n (int)      :: prinipal quantum number of orbital
        r (np.array) :: r values
    Returns:
        y (np.array) :: y values corresponding to rwf or rdf
    """
    if "RDF" in wf_type:
        return r**2.* radial_f(n, 2.*r/n)**2
    if "RWF" in wf_type:
        return radial_f(n, 2.*r/n)

##################################################################
########################## Webpage layout ########################
##################################################################

#Build the webpage layout

# Change header 
# %thing% are defined by dash and are replaced at runtime
app.index_string = r'''
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>Waveplot</title>
        <meta name="description" content="Online atomic orbital viewer">
        <meta name="keywords" content="Online atomic orbital viewer">
        <meta name="author" content="Jon Kragskow">
        {%favicon%}
        {%css%}
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            <script src="assets/3Dmol-min.js" async></script>
            {%renderer%}
        </footer>
    </body>
</html>
'''

navbar = dbc.NavbarSimple(
    id = "navbar",
    children=[
        dbc.NavItem(dbc.NavLink(id= "orb_tab", children = "Orbitals", href="/apps/app1", active=True)),
        dbc.NavItem(dbc.NavLink(id= "vib_tab", children = "Vibrations", href="/apps/app2")),
        dbc.NavItem(dbc.NavLink(id= "trans_tab", children = "Translations", href="/apps/app3")),
    ],
    brand="Waveplot",
    brand_href="/apps/app1",
    color="primary",
    dark=True,
)

orb_graph = dcc.Graph(
                    id='plot_area', 
                    style = {
                        'responsive' : 'true',
                        'height' : '580px',
                        'automargin' : 'true'
                    }
                )

orb_select = [dbc.Row(className="h-33", children =[
        dbc.Col(
            html.Div(
                children=html.H4(
                    style = {
                        'textAlign' : 'center', 
                        },
                    children = 'Orbital'
                )
            )
        ),
        dbc.Col(
            html.Div(
                children=html.H4(
                    style = {
                        'textAlign' : 'center', 
                        },
                    children = 'Function'
                )
            )
        ),
    ]),
dbc.Row([
    dbc.Col(
        html.Div(
            children=[
                dcc.Dropdown(id = 'orb_checklist', 
                    style = {
                        'textAlign' : 'left'
                    },
                    options=[
                        {"label": '1s', "value": '1s'},
                        {"label": '2s', "value": '2s'},
                        {"label": '3s', "value": '3s'},
                        {"label": '4s', "value": '4s'},
                        {"label": '5s', "value": '5s'},
                        {"label": '6s', "value": '6s'},
                        {"label": '2p', "value": '2p'},
                        {"label": '3p', "value": '3p'},
                        {"label": '4p', "value": '4p'},
                        {"label": '5p', "value": '5p'},
                        {"label": '6p', "value": '6p'},
                        {"label": '3d', "value": '3d'},
                        {"label": '4d', "value": '4d'},
                        {"label": '5d', "value": '5d'},
                        {"label": '6d', "value": '6d'},
                        {"label": '4f', "value": '4f'},
                        {"label": '5f', "value": '5f'},
                        {"label": '6f', "value": '6f'},
                    ],
                    value=['2p'],
                    multi=True, # browser autocomplete needs to be killed here, when they implement it
                    placeholder = "Orbital..."
                )
            ]
        )
    ),
    dbc.Col(
        html.Div(
            children=[
                dcc.Dropdown(
                    id = 'function_type', 
                    style = {
                        'textAlign' : 'center',
                        'display': 'block'
                    },
                    options=[
                        {
                         "label": 'Radial Distribution Function',
                         "value": 'RDF'
                        },
                        {
                         "label": 'Radial Wave Function',
                         "value": 'RWF'
                        },
                        {
                         "label": '3D Surface', 
                         "value": '3DWF'
                        }
                    ],
                    value='RDF',
                    searchable=False,
                    clearable=False
                )
            ]
        )
    )
])]

orb_save = [
    dbc.Row([
        dbc.Col([
            html.H4(
                style = {
                    'textAlign' : 'center', 
                },
                children = 'Save Options',
                id = "save_options_header"
            ),
            dbc.Tooltip(
                children="Use the camera button in the top right of the plot to save",
                target="save_options_header",
                style = {
                    'textAlign' : 'center', 
                },
            )
        ])
    ]),
    dbc.Row(
        [
            dbc.Col(
                        html.P(
                            style = {
                                'textAlign' : 'center', 
                            },
                            children = 'Output height'
                        )
                    ),
            dbc.Col(
                        html.P(
                            style = {
                                'textAlign' : 'center', 
                            },
                            children = 'Output width'
                        )
                    ),
            dbc.Col(
                        html.P(
                            style = {
                                'textAlign' : 'center', 
                            },
                            children = 'Save format'
                        )
                    ),
        ]
    ),
    dbc.Row(
        [
            dbc.Col(
                dbc.Input(
                    id = 'save_height_in',
                    placeholder=500,
                    type='number',
                    value=500
                )
            ),
            dbc.Col(
                dbc.Input(
                    id = 'save_width_in',
                    placeholder=500,
                    type='number',
                    value=500
                )
            ),
            dbc.Col(
                dcc.Dropdown(
                    id = 'save_format',
                    options=[
                        { 
                         "label": 'svg', 
                         "value": 'svg'
                        },
                        {
                         "label": 'png',
                         "value": 'png'
                        },
                        {
                         "label": 'jpeg', 
                         "value": 'jpeg'
                        }
                    ],
                    value='svg',
                    searchable=False,
                    clearable=False
                )
            ),
        ]
    ),
    ]

orb_customise_2d = [
    dbc.Row(
        dbc.Col(
                html.H4(
                    style = {'textAlign' : 'center'},
                    children = 'Plot Options')
                )
    ),
    dbc.Row(
        [
            dbc.Col(
                html.P(
                    style = {
                        'textAlign' : 'center', 
                    },
                    children = 'Gridlines'
                )
            ),
            dbc.Col(
                html.P(
                    style = {
                        'textAlign' : 'center', 
                    },
                    children = 'Lower x limit'
                )
            ),
            dbc.Col(
                html.P(
                    style = {
                        'textAlign' : 'center', 
                    },
                    children = 'Upper x limit'
                )
            )
        ]
    ),
    dbc.Row(children=
        [
            dbc.Col(
                dcc.Checklist(
                    id = 'gridlines',
                    style = {
                        'textAlign' : 'center', 
                        'verticalAlign': 'middle',
                        'horizontalAlign': 'middle'
                    },
                    options=[
                        {"label": ' x', "value": 'x'},
                        {"label": ' y', "value": 'y'}
                    ],
                    labelStyle={"display":"block"},
                    value=[],
                )
            ),
            dbc.Col(
                dbc.Input(
                    id = 'lower_x_in',
                    placeholder = 0,
                    type = 'number',
                    min = -10,
                    max = 100,
                    value = 0,
                    style = {
                        # 'width' : '50%',
                        'textAlign' : 'center', 
                        'verticalAlign': 'middle',
                        'horizontalAlign': 'middle'
                    }
                )
            ),
            dbc.Col(
                dbc.Input(
                    id = 'upper_x_in',
                    placeholder = 100,
                    type = 'number',
                    max = 100,
                    value = 100,
                    style = {
                        'textAlign' : 'center', 
                        'verticalAlign': 'middle',
                        'horizontalAlign': 'middle'
                    }
                )
            )
        ]
    ),
    dbc.Row(
        [
            dbc.Col(
                html.P(
                    style = {
                        'textAlign' : 'center', 
                    },
                    children = 'Line width'
                )
            ),
            dbc.Col(
                html.P(
                    style = {
                        'textAlign' : 'center', 
                    },
                    children = 'Text size'
                )
            ),
            dbc.Col(
                html.P(
                    style = {
                        'textAlign' : 'center', 
                    },
                    children = 'Plot colours'
                )
            )
        ]
    ),
    dbc.Row(
        [
            dbc.Col(
                dcc.Slider(
                    id = 'linewidth_slider',
                    min=1,
                    max=10,
                    step=0.5,
                    value=5,
                )
            ),
            dbc.Col(
                dcc.Slider(
                    id = 'text_size_slider',
                    min=15,
                    max=25,
                    step=0.5,
                    value=19,
                )
            ),
            dbc.Col(
                dcc.Dropdown(
                    id = 'colours_2d',
                    options=[
                            { 
                             "label": 'Standard', 
                             "value": 'normal'
                            },
                            {
                             "label": 'Tol',
                             "value": 'tol'
                            },
                            {
                             "label": 'Wong', 
                             "value": 'wong'
                            }
                    ],
                    value='normal',
                    searchable=False,
                    clearable=False
                )
            ),
        ]
    )
]

orb_customise_3d = [ 
    dbc.Row([
        dbc.Col([
            html.H4(
                style = {
                    'textAlign' : 'center', 
                },
                children = 'Plot Options'
            )
        ])
    ]),
    dbc.Row(
        [dbc.Col(
            html.P(
                style = {
                    'textAlign' : 'center', 
                },
                children = 'Colours'
            )
        ),
        dbc.Col(
            html.P(
                style = {
                    'textAlign' : 'center', 
                },
                children = 'Cutaway'
            )
        )
        ]
    ),
    dbc.Row(
        [dbc.Col(
            dcc.Dropdown(
                id = 'colours_3d',
                options=[
                        { 
                         "label": 'Yellow-Purple', 
                         "value": 'yp'
                        },
                        {
                         "label": 'Red-Blue',
                         "value": 'rb'
                        },
                        {
                         "label": 'Green-Orange', 
                         "value": 'go'
                        }
                ],
                style = {}, 
                value='yp',
                searchable=False,
                clearable=False
            )
        ),
        dbc.Col(
            dcc.RadioItems(id = 'cutaway_in', 
                style = {
                    'textAlign' : 'center', 
                    'verticalAlign': 'middle',
                    'horizontalAlign': 'middle'
                },
                options=[
                    {
                     "label": '  None', 
                     "value": 1.0
                    },
                    # {
                    #  "label": '1/4',
                    #  "value": 0.5
                    # },
                    {
                     "label": '  1/2',
                     "value": 0.
                    },
                ],
                value=1.0,
                labelStyle={"display":"block"}
            )
        )
    ])
]

orb_customise = [
    html.Div(
        id = "orb_customise_2d",
        children = orb_customise_2d,
        style = {}
    ),
    html.Div(
        id = "orb_customise_3d",
        children = orb_customise_3d,
        style = {}
    )]
            

orb_options = orb_select + orb_customise + orb_save

orb_body = dbc.Container(
    dbc.Row(
        [
            dbc.Col(orb_graph, style={"height" : "100%"},className="h-100"),
            dbc.Col(orb_options, style={"height" : "100%"},className="h-100"),
        ],
        className="h-100",
    ),
    style = {"height" : "100vh"}

)

footer = html.Footer(
    style = {
        'textAlign'       : 'center', 
        'font-size'       : 'smaller',
        'color'           : 'white',
        'background-color': '#307cf6'
    }, 
    children=[
        html.P(
            children = [
            'Jon Kragskow'
            ]
        ),
        html.A(
            href = 'https://www.kragskow.com/',
            style = {
                'color':'white'
            },
            children = 'https://www.kragskow.com/'
        )
    ]
)

##################################################################
########################## Webpage Main ##########################
##################################################################

# Layout of webpage
layout = html.Div(
    children=[
        navbar,
        orb_body,
        footer,
    ]
)

def toggle_pob(wf_type) :

    if on_off == 'on' :
        return {
            'display' : 'inline-block' , 
            'padding-left':'5%', 
            'padding-bottom':'5%' , 
            'padding-right':'5%', 
            'width' : '90%', 
            'borderStyle' : 'solid', 
            'borderWidth' : '0px',
            'textAlign' : 'center'
        }
    elif on_off == None:
        return {
            'display' : 'none' , 
        }     
    else :
        return {
            'display' : 'none' , 
            'padding-left':'5%', 
            'padding-bottom':'5%' , 
            'padding-right':'5%', 
            'width' : '90%', 
            'borderStyle' : 'solid', 
            'borderWidth' : '0px',
            'textAlign' : 'center'
        }


def orb_customise_2d(wf_type) :

    #3D
    if "3" in wf_type:
        style = {'display' : 'none' , }

    # 2D
    else:
        style = {}


    return style

def orb_customise_3d(wf_type) :

    #3D
    if "3" in wf_type:
        style = {}

    # 2D
    else:
        style = {'display' : 'none' , }


    return style

def orb_checklist(wf_type):


    if "3" in wf_type:
        checklist = [
            {"label": "1s", "value": "1s"},
            {"label": "2s", "value": "2s"},
            {"label": "3s", "value": "3s"},
            {"label": "4s", "value": "4s"},
            {"label": "5s", "value": "5s"},
            {"label": "6s", "value": "6s"},
            {"label": "2p", "value": "2p"},
            {"label": "3p", "value": "3p"},
            {"label": "4p", "value": "4p"},
            {"label": "5p", "value": "5p"},
            {"label": "6p", "value": "6p"},
            {"label": "3dz²", "value": "3d_z2"},
            {"label": "4dz²", "value": "4d_z2"},
            {"label": "5dz²", "value": "5d_z2"},
            {"label": "6dz²", "value": "6d_z2"},
            {"label": "3dxy", "value": "3dxy"},
            {"label": "4dxy", "value": "4dxy"},
            {"label": "5dxy", "value": "5dxy"},
            {"label": "6dxy", "value": "6dxy"},
            {"label": "4fz³", "value": "4fz3"},
            {"label": "5fz³", "value": "5fz3"},
            {"label": "6fz³", "value": "6fz3"},
            {"label": "4fxyz", "value": "4fxyz"},
            {"label": "5fxyz", "value": "5fxyz"},
            {"label": "6fxyz", "value": "6fxyz"},
            {"label": "4fyz²", "value": "4fyz2"},
            {"label": "5fyz²", "value": "5fyz2"},
            {"label": "6fyz²", "value": "6fyz2"},
        ]
    else:
        checklist = [
            {"label": "1s", "value": "1s"},
            {"label": "2s", "value": "2s"},
            {"label": "3s", "value": "3s"},
            {"label": "4s", "value": "4s"},
            {"label": "5s", "value": "5s"},
            {"label": "6s", "value": "6s"},
            {"label": "2p", "value": "2p"},
            {"label": "3p", "value": "3p"},
            {"label": "4p", "value": "4p"},
            {"label": "5p", "value": "5p"},
            {"label": "6p", "value": "6p"},
            {"label": "3d", "value": "3d"},
            {"label": "4d", "value": "4d"},
            {"label": "5d", "value": "5d"},
            {"label": "6d", "value": "6d"},
            {"label": "4f", "value": "4f"},
            {"label": "5f", "value": "5f"},
            {"label": "6f", "value": "6f"},
        ]

    return checklist


##################################################################
############################ Callback ############################
##################################################################

#Callback which defines what changes (e.g. the plot) and what causes 
# the change (e.g. a checkbox being pressed)

orbitals_input = [ddep.Input('orb_checklist', "value"),
               ddep.Input('function_type', "value"), 
               ddep.Input('linewidth_slider',"value"), 
               ddep.Input('text_size_slider',"value"), 
               ddep.Input('gridlines',"value"), 
               ddep.Input('upper_x_in', "value"),
               ddep.Input('lower_x_in', "value"),
               ddep.Input('save_format', "value"), 
               ddep.Input('save_height_in', "value"),
               ddep.Input('save_width_in', "value"),
               ddep.Input('colours_2d',"value"),
               ddep.Input('colours_3d',"value"),
               ddep.Input('cutaway_in',"value")]

@app.callback([ddep.Output('plot_area', 'figure'),
               ddep.Output('plot_area', 'config'),
               ddep.Output('orb_customise_2d', 'style'),
               ddep.Output('orb_customise_3d', 'style'),
               ddep.Output('orb_checklist', 'options')],
               orbitals_input
              )

def update_app(orbitals, wf_type, linewidth, text_size, gridlines,
               x_up, x_low, save_format, save_height, save_width, 
               colours_2d, colours_3d, cutaway):
    """
    Updates the app, given the current state of the UI
    All inputs correspond (in the same order) to the list 
    of ddep.Input a few lines above^^^

    Input:
        orbitals (list, string)   :: names of orbitals
        wf_type (string)          :: type of wf 
        linewidth (float)         :: linewidth for 2d plot
        text_size (float)         :: plot label text size for 2d plot
        gridlines (list, string)  :: yes or no to gridlines on either axis
        x_up (float)              :: upper x limit for 2d plot
        x_low (float)             :: lower x limit for 2d plot
        save_format (string)      :: save format for plot
        save_height (float)       :: height of saved image
        save_width  (float)       :: width of saved image
        colour_2d (string)        :: colour style for 2d plot
        colour_3d (string)        :: colour style for 3d plot
        cutaway (float)           :: controls 3d slicing of orbitals

    Returns:
        figure (dict)             :: dictionary item for list of go.XX objects
        config (dict)             :: dictionary item for go.layout object, e.g. linewidth, colormap...
        orb_customise_2d (dict)   :: dictionary of style options for 2d plot UI input area
        orb_customise_3d (dict)   :: dictionary of style options for 2d plot UI input area
        orb_checklist  (dict)     :: dictionary of dictionaries containing label value pairs for orbitals checkbox

    """

    return [orb_fig(orbitals, x_up, x_low, wf_type, linewidth, colours_2d, colours_3d, cutaway, gridlines, text_size),
    orb_modebar(save_format, save_height, save_width, wf_type, orbitals),
    orb_customise_2d(wf_type), 
    orb_customise_3d(wf_type), 
    orb_checklist(wf_type)]

def orb_fig(orbitals, x_up, x_low, wf_type, linewidth, colours_2d, colours_3d, cutaway, gridlines, text_size):

    # Nothing to plot - exit
    if len(orbitals) == 0:
        output = {
            'data'   : None,
            'layout' : ax_null()
            }
        return output

    # Set limits to default if needed
    if x_low == None or x_low < 0:
        x_low = 0
    if x_up == None:
        x_up = x_low + 100

    #Turn on x axis gridlines
    if 'x' in gridlines:
        x_grid = True
    else:
        x_grid = False

    #Turn on y axis gridlines
    if 'y' in gridlines:
        y_grid = True
    else:
        y_grid = False

    y_labels = {
        "RDF" : "Radial distribution function 4𝜋𝑟²𝑅(𝑟)²",
        "RWF" : "Radial wavefunction 𝑅(𝑟)"
    }

    if "3" in wf_type:
        data, x_up, x_low = orb_plot_3d(orbitals[0], colours_3d, cutaway)
        layout = orb_ax_3d(x_up, x_low)
    else:
        layout = orb_ax_2d(y_labels[wf_type], text_size, x_grid, y_grid, x_up, x_low)
        data = orb_plot_2d(orbitals, x_up, x_low, wf_type, linewidth, colours_2d)

    output = {
        "data" : data,
        "layout" : layout
    }

    return output

def orb_plot_2d(orbitals, x_up, x_low, wf_type, linewidth, colours_2d):

    # Define colour lists
    # Paul Tol list of colourblindness friendly colours
    # https://personal.sron.nl/~pault/
    tol_cols = [
        'rgb(0  , 0  , 0)',
        'rgb(230, 159, 0)',
        'rgb(86 , 180, 233)',
        'rgb(0  , 158, 115)',
        'rgb(240, 228, 66)',
        'rgb(0  , 114, 178)',
        'rgb(213, 94 , 0)',
        'rgb(204, 121, 167)'
    ]
    # Bang wong list of colourblindness friendly colours
    # https://www.nature.com/articles/nmeth.1618
    wong_cols = [
        'rgb(51 , 34 , 136)',
        'rgb(17 , 119, 51)',
        'rgb(68 , 170, 153)',
        'rgb(136, 204, 238)',
        'rgb(221, 204, 119)',
        'rgb(204, 102, 119)',
        'rgb(170, 68 , 153)',
        'rgb(136, 34 , 85)'
    ]
    # Default list of colours is plotly's safe colourlist
    def_cols = pc.qualitative.Safe

    if colours_2d == 'tol':
        cols = tol_cols + wong_cols + def_cols
    elif colours_2d == 'wong':
        cols = wong_cols + def_cols + tol_cols
    else:
        cols = def_cols + tol_cols + wong_cols

    traces = []

    curr_ymax = 0.
    x = np.linspace(x_low,x_up,1000)

    # Plot each requested function
    for it, orbital in enumerate(orbitals):
        # Get orbital n value and name
        n, l = name_to_qn(orbital)

        if l == 's':
            y = calc_radial_s(n, x, wf_type)

        elif l == 'p':
            y = calc_radial_p(n, x, wf_type)

        elif l == 'd':  
            y = calc_radial_d(n, x, wf_type)

        elif l == 'f':
            y = calc_radial_f(n, x, wf_type)

        traces.append(
            go.Scatter(
                x = x,
                y = y,
                line = dict(width = linewidth),
                name = orbital,
                hoverinfo = 'none',
                marker={"color":cols[it]}
            )
        )

    return traces


def orb_plot_3d(orb_name, colour_name, cutaway):
    """
    Plots isosurfaces of atomic orbitals

    Input:
        orb_name (string)       :: name of orbital
        colour_name (string)    :: name of colour palette
        cutaway (float)         :: % of orbital to keep/plot
    Returns:
        data (plotly go object) :: plotly fig object with new graph added
        upper (float)           :: upper bound of all axis limits 
        lower (float)           :: lower bound of all axis limits 
    """

    # Get colours of lobes
    colours = set_3d_colour(colour_name)
    
    # Get orbital n value and name
    n, l = name_to_qn(orb_name)

    if l == 's':

        x, y, z, wav, upper, lower, ival = calc_s_orb(n, cutaway)

    elif l == 'p':

        data, upper, lower= plot_p_orb(n, cutaway, colours)

    elif l == 'd':  

        if 'xy' in orb_name:

            x, y, z, wav, upper, lower, ival = calc_dxy_orb(n, cutaway)

        else:

            x, y, z, wav, upper, lower, ival = calc_dz_orb(n, cutaway)

    elif l == 'f':

        if 'xyz' in orb_name:

            x, y, z, wav, upper, lower, ival = calc_fxyz_orb(n, cutaway)


        elif 'yz2' in orb_name:

            x, y, z, wav, upper, lower, ival = calc_fyz2_orb(n, cutaway)

        else:

            x, y, z, wav, upper, lower, ival = calc_fz_orb(n, cutaway)

    if l != "p":

        data = go.Isosurface(
            x=x.flatten(),
            y=y.flatten(),
            z=z.flatten(),
            value=wav.flatten(),
            isomin=-ival,
            isomax=ival,
            flatshading=False,
            caps=dict(x_show=False, y_show=False, z_show=False),
            showscale=False,
            colorscale=colours
            )

    return [data], upper, lower

def orb_ax_2d(y_label, text_size, x_grid, y_grid, x_up, x_low):

    layout = go.Layout(
                xaxis = {
                    'autorange' : True,
                    'showgrid'  : x_grid,
                    'zeroline'  : False,
                    'showline'  : True,
                    'range'     : [x_low, x_up],
                    'title' : {
                        'text' : "Distance (𝑎₀)",
                        'font' :{'size' : text_size} 
                    },
                    'ticks' :'outside',
                    'tickfont' :{'size' : text_size},
                    'showticklabels' : True
                },
                yaxis = {
                    'autorange'  : True,
                    'showgrid'   : y_grid,
                    'zeroline'   : False,
                    'fixedrange' : True,
                    'title' :{
                        'text' : y_label,
                        'font' :{
                            'size' : text_size
                        }
                    },
                    'title_standoff' : 100,
                    'showline' :True,
                    'ticks' :'outside',
                    'tickfont' :{'size' : text_size},
                    'showticklabels' :True
                },
                legend = {
                    'x' : 0.8,
                    'y' : 1,
                    'font' :{
                        'size' : text_size - 3
                    }
                },
                margin=dict(l=90, r=30, t=30, b=60),
    )

    return layout

def orb_ax_3d(upper=0, lower=0):

    layout = go.Layout(
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
                backgroundcolor='rgb(255, 255,255)',
                range=[lower, upper]
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
                backgroundcolor='rgb(255, 255,255)',
                range=[lower, upper]
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
                       backgroundcolor='rgb(255, 255,255)',
                       range=[lower, upper]
                      ),
            aspectratio=dict(
                x=1,
                y=1,
                z=1
            ),
        ),
        margin=dict(l=20, r=30, t=30, b=20),
    )

    return layout

def ax_null():

    #Turn off axis gridlines

    return go.Layout(
                     xaxis = {
                              'autorange' : True,
                              'showgrid'  : False,
                              'zeroline'  : False,
                              'showline'  : False,
                              'range'     : [0, 1],
                              'showticklabels' : False
                             },
                     yaxis = {
                              'autorange' :True,
                              'showgrid' :False,
                              'zeroline' :False,
                              'showline' :False,
                              'showticklabels' :False
                             },
                     margin=dict(l=90, r=30, t=30, b=60),
    )


def orb_modebar(save_format, save_height, save_width, wf_type, orbitals):

    # No plot, no modebar
    if len(orbitals) == 0:
        options = {"displayModeBar" : False}
        return options

    if "RDF" in wf_type:
        file_name = "radial_distribution_function"
    elif "RWF" in wf_type:
        file_name = "radial_wavefunction"
    elif "3DWF" in wf_type:
        file_name = "{}_orbital".format(orbitals[0])

    options = {
        'toImageButtonOptions':{
            'format': save_format, 
            'filename': file_name,
            'height': save_height,
            'width': save_width,
        }, 
        'modeBarButtonsToRemove':[
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

    return options