import numpy as np
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import dash.dependencies as ddep
from subprocess import call
import json

#####################################################################################  Experimental  #####################################################################################

#####################################################################################
################################ Orbital Mathematics ################################
#####################################################################################


############################ p Orbital functions #########################

def PrincipleAxisFuncp(n,r,f,c):
    return n*np.exp(r/n)*c*np.sqrt(np.pi/3.)/abs(f)

def radial_p(Orb,r):
    # !!! Radial Wavefunctions of p orbitals
    if Orb=='2p': return 1./(2.*np.sqrt(6.))
 
    if Orb=='3p': return (4. - 2.*r/3.)/(9*np.sqrt(6.))

    if Orb=='4p': return (20. - 5.*r + (r/2.)**2.)/(32.*np.sqrt(15.))

    if Orb=='5p': return (120. - 90.*(2*r/5.) + 18.*(2*r/5.)**2. - (2.*r/5.)**3.)/(150.*np.sqrt(30.))

    if Orb=='6p': return (840. - 840.*r/3. + 252.*(r/3.)**2. - 28.*(r/3.)**3. + (r/3.)**4.)/(432.*np.sqrt(210.))

############################ d Orbital functions #########################

def radial_d(Orb,r): 
    # !!! Radial Wavefunctions of d orbitals
    if Orb=='3d': out = 1./(9.*np.sqrt(30.))

    if Orb=='4d': out = (6.-r/2.)/(96.*np.sqrt(5.))

    if Orb=='5d': out = (42. - 28.*r/5. + (2.*r/5.)**2.)/(150.*np.sqrt(70.))

    if Orb=='6d': out = (336. - 168.*(r/3.) + 24.*(r/3.)**2. - (r/3.)**3.)/(864.*np.sqrt(105.))
    return out

############################f Orbital Wavefunction functions#########################

#def radial_f(Orb,r): 
#    if Orb=='4f': out = 
#
#    if Orb=='5f': out = 
#
#    if Orb=='6f': out = 
#    return out



def r_domain(Orb):
    # !!! Returns an array containing the bounds of r for each lobe of the requested orbital

#####2p orbital#####
    if Orb == '2p':
        lobe_1_lower = 0.030542415534
        lobe_1_upper = 11.973154175074
       
        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper]])

#####3p orbital#####
    if Orb == '3p':

        lobe_1_lower = 0.0521008719026 
        lobe_1_upper =  5.6457958052515
        
        lobe_2_lower = 6.4019092339593
        lobe_2_upper = np.float128(20.7317063563592)

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper]])
    
#####4p orbital#####

    if Orb=='4p':

        lobe_1_lower = 0.0791782710855
        lobe_1_upper = 5.0739651680262
        
        lobe_2_lower = 6.0725425741786
        lobe_2_upper = 12.8640838577464
        
        lobe_3_lower = 16.5661209026166
        lobe_3_upper = 30.3121031661214

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper],[ lobe_3_lower, lobe_3_upper]])
    
#####5p orbital#####

    if Orb=='5p':
    
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


        lobe_1_lower = 0.0138230698659
        lobe_1_upper = 5.1856979119147
        
        lobe_2_lower = 5.3499953596572
        lobe_2_upper = 12.5530614690789
    
        lobe_3_lower = 13.0502450037780
        lobe_3_upper = 23.6116572475096

        lobe_4_lower = 24.7635128466900 
        lobe_4_upper = 40.5262363941473
        
        lobe_5_lower = 43.0775334433834
        lobe_5_upper = np.float128(84.65464)
    
        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper],[ lobe_3_lower, lobe_3_upper],[ lobe_4_lower, lobe_4_upper],[ lobe_5_lower, lobe_5_upper]])

#####3d orbital#####

    if Orb == '3d':
        lobe_1_lower = 0.0305424
        lobe_1_upper = 11.973154-0.03054255
       
        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper]])

#####4d orbitals other than z2#####

    if Orb == '4d_nz2':

        lobe_1_lower = 1.0606
        lobe_1_upper = 12.3325
        
        lobe_1_lower = 13.515
        lobe_2_upper =  17.9724

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper]])
    
#####4d_z2 orbital#####

    if Orb == '4d_z2':

        lobe_1_lower = 0.412143781754179
        lobe_1_upper = 11.4749904407017
        
        lobe_2_lower = 0.39411
        lobe_2_upper = 11.382175
        
        lobe_3_lower = 12.23238573
        lobe_3_upper = 29.08904
        
        lobe_4_lower = 12.11503
        lobe_4_upper = 33.23954

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper],[ lobe_3_lower, lobe_3_upper],[ lobe_4_lower, lobe_4_upper]])

#####5d orbitals other than z2#####

    if Orb=='5d_nz2':

        lobe_1_lower = 0.6156
        lobe_1_upper = 10.3787
        
        lobe_2_lower = 11.3617
        lobe_2_upper = 12.1429
        
        lobe_3_lower = 25.5091
        lobe_3_upper = 33.2447

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper],[ lobe_3_lower, lobe_3_upper]])
    
#####5d_z2 orbital#####

    if Orb=='5d_z2':
    
        lobe_1_lower = 0.330 
        lobe_1_upper = 10.42595

        lobe_2_lower = 0.5
        lobe_2_upper = 10.129142

        lobe_3_lower = 11.15481
        lobe_3_upper = 12.22665

        lobe_4_lower = 11.0185
        lobe_4_upper = 13.4314

        lobe_5_lower = 24.4976
        lobe_5_upper = 39.718

        lobe_6_lower = 24.8949
        lobe_6_upper = 33.57544

        lobe_domains = np.array([[lobe_1_lower, lobe_1_upper],[ lobe_2_lower, lobe_2_upper],[ lobe_3_lower, lobe_3_upper],[ lobe_4_lower, lobe_4_upper],[ lobe_5_lower, lobe_5_upper],[ lobe_6_lower, lobe_6_upper]])
    

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
    return lobe_domains

############################Orbital Wavefunction Calculation#########################


def OrbCalc(Orb):

    letter = Orb.strip('Orb')
    n = float(letter.strip('p'))

    if n < 5:
        c = 0.003
    else:
        c = 0.0003

    #Number of steps for angle and r
    r_steps = 80
    angle_steps = 50

    #Arrays for x,y,z coordinates
    # z has two separate arrays as it has two lobes separated by z axis
    x = np.zeros([angle_steps, r_steps])
    y = np.zeros([angle_steps, r_steps])
    x_y_no_angle = np.zeros([r_steps])
    zp = np.zeros([angle_steps, r_steps])
    zm = np.zeros([angle_steps, r_steps])

    data = []

    #Plot different radial portions of orbital

    flag = True

    #Get bounds of r for each lobe of orbital
    orb_r_bounds = r_domain(Orb)

    #Calculate coordinates of isosurface
    # loop over each lobe of the orbital
    
    for lobe_it in np.arange(np.shape(orb_r_bounds)[0]):

        #Array of r values for each lobe
        r = np.linspace(orb_r_bounds[lobe_it,0], orb_r_bounds[lobe_it,1], num = r_steps)


        #Calculate z(r)
        z_of_r = PrincipleAxisFuncp(n,r,radial_p(Orb,r),c)

        #Calculate angle independent term of x and y
        x_y_no_angle = np.sqrt(r**2. - z_of_r**2.)
        
        #Array of angle values
        ang = np.linspace(0., np.pi*2., num = angle_steps)

        #Calculate x and y by multiplying by angular factor
        #And set z coordinate for every angle and r value
        for r_it in np.arange(r_steps):
            for ang_it in np.arange(angle_steps):
                x[ang_it,r_it] = x_y_no_angle[r_it]*np.cos(ang[ang_it])
                y[ang_it,r_it] = x_y_no_angle[r_it]*np.sin(ang[ang_it])
                zp[ang_it,r_it] = z_of_r[r_it]

        zm = -zp

        #Switch lobe colours when plotting
        flag = not flag

        if flag == True:
            data0 = go.Surface(x=x, y=y, z=zp, colorscale= [[0, 'rgb(109,0,157)'], [1,'rgb(109,0,157)']],showscale=False)
            data1 = go.Surface(x=x, y=y, z=zm, colorscale= [[0, 'rgb(255,204,51)'], [1,'rgb(255,204,51)']],showscale=False)

        elif flag == False:
            data0 = go.Surface(x=x, y=y, z=zm, colorscale= [[0, 'rgb(109,0,157)'], [1,'rgb(109,0,157)']],showscale=False)
            data1 = go.Surface(x=x, y=y, z=zp, colorscale= [[0, 'rgb(255,204,51)'], [1,'rgb(255,204,51)']],showscale=False)
        
        data.append(data0)
        data.append(data1)     

    return data


###Functions for Radial Wave Functions

def s_1(r,mode, zeff):
    rho = 2.*r/1.  * zeff
    O1s = 2.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O1s**2.
    if mode==2:
        return O1s

def s_2(r,mode, zeff):
    rho = 2.*r/2. * zeff
    O2s = 1./(2.*np.sqrt(2.))*(2.-rho)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O2s**2.
    if mode==2:
        return O2s

def s_3(r,mode, zeff):
    rho = 2.*r/3. * zeff
    O3s = 1./(9.*np.sqrt(3.))*(6.-6.*rho+rho**2.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O3s**2.
    if mode==2:
        return O3s

def s_4(r,mode, zeff):
    rho = 2.*r/4. * zeff
    O4s = (1./96.)*(24.-36.*rho+12.*rho**2.-rho**3.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O4s**2.
    if mode==2:
        return O4s

def s_5(r,mode, zeff):
    rho = 2.*r/5. * zeff
    O5s = (1./(300.*np.sqrt(5.)))*(120.-240.*rho+120.*rho**2.-20.*rho**3.+rho**4.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O5s**2.
    if mode==2:
        return O5s

def s_6(r,mode, zeff):
    rho = 2.*r/6. * zeff
    O6s = (1./(2160.*np.sqrt(6.)))*(720.-1800.*rho+1200.*rho**2.-300.*rho**3.+30.*rho**4.-rho**5.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O6s**2.
    if mode==2:
        return O6s

def p_2(r,mode, zeff):
    rho = 2.*r/2. * zeff
    O2p = 1./(2.*np.sqrt(6.))*rho*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O2p**2.
    if mode==2:
        return O2p

def p_3(r,mode, zeff):
    rho = 2.*r/3. * zeff
    O3p = 1./(9.*np.sqrt(6.))*rho*(4.-rho)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O3p**2.
    if mode==2:
        return O3p

def p_4(r,mode, zeff):
    rho = 2.*r/4. * zeff
    O4p = 1./(32.*np.sqrt(15.))*rho*(20.-10.*rho+rho**2.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O4p**2.
    if mode==2:
        return O4p

def p_5(r,mode, zeff):
    rho = 2.*r/5. * zeff
    O5p = 1./(150.*np.sqrt(30.))*rho*(120.-90.*rho+18.*rho**2.-rho**3.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O5p**2.
    if mode==2:
        return O5p

def p_6(r,mode, zeff):
    rho = 2.*r/6. * zeff
    O6p = 1./(432.*np.sqrt(210.))*rho*(840.-840.*rho+252.*rho**2.-28.*rho**3.+rho**4.)*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O6p**2.
    if mode==2:
        return O6p

def d_3(r,mode, zeff):
    rho = 2.*r/3. * zeff
    O3d = 1./(9.*np.sqrt(30.))*rho**2.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O3d**2.
    if mode==2:
        return O3d

def d_4(r,mode, zeff):
    rho = 2.*r/4. * zeff
    O4d = 1./(96.*np.sqrt(5.))*(6.-rho)*rho**2.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O4d**2.
    if mode==2:
        return O4d

def d_5(r,mode, zeff):
    rho = 2.*r/5. * zeff
    O5d = 1./(150.*np.sqrt(70.))*(42.-14.*rho+rho**2)*rho**2.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O5d**2.
    if mode==2:
        return O5d

def d_6(r,mode, zeff):
    rho = 2.*r/6. * zeff
    O6d = 1./(864.*np.sqrt(105.))*(336.-168.*rho+24.*rho**2.-rho**3.)*rho**2.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O6d**2.
    if mode==2:
        return O6d

def f_4(r,mode, zeff):
    rho = 2.*r/4. * zeff
    O4f = 1./(96.*np.sqrt(35.))*rho**3.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O4f**2.
    if mode==2:
        return O4f

def f_5(r,mode, zeff):
    rho = 2.*r/5. * zeff
    O5f = 1./(300.*np.sqrt(70.))*(8.-rho)*rho**3.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O5f**2.
    if mode==2:
        return O5f

def f_6(r,mode, zeff):
    rho = 2.*r/5. * zeff
    O6f = 1./(2592.*np.sqrt(35.))*(rho**2.-18.*rho+72.)*rho**3.*np.exp(-rho/2.) * zeff**1.5
    if mode==1:
        return r**2.*O6f**2.
    if mode==2:
        return O6f


functiondict = {'1s': s_1, '2s': s_2, '3s': s_3, '4s': s_4, '5s': s_5, '6s': s_6, '2p': p_2,'3p': p_3,'4p': p_4,'5p': p_5,
    '6p': p_6, '3d': d_3, '4d': d_4, '5d': d_5, '6d': d_6, '4f': f_4, '5f': f_5, '6f': f_6}


##################################################################################################################################
##################################################################################################################################
##################################################################################################################################



#Build the webpage layout
app = dash.Dash(__name__)

app.css.append_css({"external_url": "https://max.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css"})

app.css.append_css({
   'external_url': (
       'assets/external.css'
   )
})

colors = {
    'background': '#c1c1c1',
    'Header': '#c1c1c1',
    'text': '#7FDBFF'
}


app.layout = html.Div(children=[

    html.H1(style = {'textAlign' : 'center', 'fontFamily' : 'sans-serif'},children = 'Orbiplot'),
    html.P(style = {'textAlign' : 'center', 'fontFamily' : 'sans-serif'},children = 'Jon Kragskow 2019'),
    html.P(style = {'textAlign' : 'center', 'fontFamily' : 'sans-serif'},children = 'Please Cite OrbiPlot if you use it!'),
    html.P(style = {'textAlign' : 'center', 'fontFamily' : 'sans-serif'},children = 'www.kragskow.com'),

    html.Div(className = 'below_title', style = {'width' : '100%', 'height' : '100%'}, children= [

        html.Div(className = 'graph box', style = {'float' : 'left', 'width' : '59%'}, children= [
            dcc.Graph(id='RDF_Graph', style = {'autosize' : 'true', 'width' : '100%', 'height' : '600px'}),
        ]),
                
        html.Div(className = 'navbar', style = {'autosize' : 'true', 'float' : 'right','display' : 'inline-block' , 'width' : '39%','textAlign' : 'left', 'borderStyle' : 'solid', 'borderWidth' : '0px'}, children=[

            html.Div(className = 'top_box', style = {'display' : 'inline-block' , 'padding-left':'5%', 'padding-right':'5%', 'padding-bottom':'5%','width' : '90%', 'borderwidth' : '0px','textAlign' : 'center'},children = [

                html.Div(className = 'orbitals box', style = {'width' : '47%' , 'padding-right' : '2%' ,'float' : 'left', 'display' : 'inline-block'}, children = [

                    html.H2(children = 'Orbitals'),

                    dcc.Checklist(id = 'OrbCheck',
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
                        value=['6p'],
                        labelStyle={'maxwidth' : '20px'}
                    ),

                ]),

                html.Div(className = 'WF box', style = {'width' : '47%' , 'padding-left' : '2%' ,'float' : 'right', 'display' : 'inline-block'}, children = [

                    html.H2(children = 'Orbital Options'),

                    dcc.RadioItems(id = 'WFType',
                        options=[
                            {'label': 'Radial Distribution Function', 'value': 'RDF'},
                            {'label': 'Radial Wave Function', 'value': 'RWF'},
                            {'label': 'Experimental 3D Wave Function', 'value': '3DWF'}
                        ],
                        value='RWF',
                        labelStyle={'width' : '10px'}
                    ),


                    html.P(children = 'Z Effective'),
                    dcc.Slider(
                            id = 'ZeffSlider',
                            min=1,
                            max=100,
                            step=1,
                            value=1,
                    ),
        
                    ]),



                ]),

            html.P(children = ''),

            html.Div(className = 'Save Options Box', style = {'textAlign' : 'center', 'padding-bottom':'5%'}, children=[
                
                    
            html.H2(children = 'Save Options'),

            html.Div(className = 'leftsavebox', style = {'float' : 'left' ,'display' : 'inline-block', 'width' : '45%' ,'padding-right' : '3%'}, children = [
                dcc.Dropdown(
                    id = 'PlotFormatDropdown',
                    options=[
                        {'label': 'svg', 'value': 'svg'},
                        {'label': 'png', 'value': 'png'},
                        {'label': 'jpeg', 'value': 'jpeg'}
                    ],
                value='svg'
                )
            ]),

            html.Div(className = 'rightsavebox', style = {'float' : 'right' ,'display' : 'inline-block', 'width' : '45%', 'padding-left' : '3%'}, children = [
                html.P('Output Height'),
                dcc.Input(
                    id = 'PlotHeightInput',
                    placeholder='Height',
                    type='number',
                    value=500
                ),
                html.P('Output Width'),
                dcc.Input(
                    id = 'PlotWidthInput',
                    placeholder='Width',
                    type='number',
                    value=700
                )
            ])
        ]),

        html.Div(className = 'bottom_box', style = {'display' : 'inline-block' , 'width' : '100%', 'borderStyle' : 'solid', 'borderWidth' : '0px','textAlign' : 'center'}, children=[

            html.H2(children = 'Plot Options'),
            html.P(children = ''),
            
            html.Div(className = 'bottom_box_below_title', style = {'display' : 'inline-block' , 'padding-left':'5%', 'padding-bottom':'5%' , 'padding-right':'5%', 'width' : '90%', 'borderStyle' : 'solid', 'borderWidth' : '0px','textAlign' : 'center'}, children=[


                html.Div(className = 'gridlines box', style = {'width' : '29%' , 'padding-right' : '2%' ,'float' : 'left', 'display' : 'inline-block'}, children = [

                    html.P(children = 'Gridlines'),

                    dcc.Checklist(id = 'xgridcheck',
                        options=[
                            {'label': ' x-Axis ', 'value': 'xgridval'}
                        ],
                        value=[],
                        ),

                    dcc.Checklist(id = 'ygridcheck',
                        options=[
                            {'label': ' y-Axis ', 'value': 'ygridval'}
                        ],
                        value=[],
                        ),

                    ]),
                html.Div(className = 'limits box', style = {'width' : '29%' , 'padding-right' : '2%' ,'float' : 'center', 'display' : 'inline-block'}, children = [

                    html.P(children = 'Lower x Limit'),
                    dcc.Input(
                        id = 'LowerxInput',
                        placeholder='lower x limit',
                        type='number',
                        value=0
                    ),
    
                    html.P(children = 'Upper x Limit'),
                    dcc.Input(
                        id = 'UpperxInput',
                        placeholder='Upper x limit',
                        type='number',
                        value=40
                    ),

                ]),

                html.Div(className = 'sliders box', style = {'width' : '29%' , 'padding-right' : '2%' ,'float' : 'right','display' : 'inline-block'}, children = [

                    html.P(children = 'Line Width'),
                    dcc.Slider(
                        id = 'LineWidthSlider',
                        min=1,
                        max=10,
                        step=0.5,
                        value=5,
                    ),
                    html.P(children = 'Text Size '),
                    dcc.Slider(
                        id = 'TextSizeSlider',
                        min=5,
                        max=25,
                        step=0.5,
                        value=20,
                    )
                ]),
            ]),
        ])
    ])  
])
])




#Callback which defines what changes ('figure') and what causes the change (Checkbox being pressed)
@app.callback(ddep.Output('RDF_Graph', 'figure'),[ddep.Input('OrbCheck', 'value'),ddep.Input('WFType', 'value'), 
    ddep.Input('LineWidthSlider','value'), ddep.Input('TextSizeSlider','value'), ddep.Input('xgridcheck','value'), 
    ddep.Input('ygridcheck','value'), ddep.Input('ZeffSlider','value'),ddep.Input('UpperxInput', 'value'),ddep.Input('LowerxInput', 'value')])

#Function which is called after the element is pressed
def UpdatePlot(Orbitals, WFType, Thickness, TextSize, xgridinput, ygridinput, zeff, upperxlim, lowerxlim):


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


    #Read type of wavefunction being requested and set axis names
    if WFType == 'RDF':
        WFFlag = 1
        WFName = 'Radial Distribution Function' #' 4\pi r R(r)'
    if WFType == 'RWF':
        WFFlag = 2
        WFName = 'Radial Wavefunction'#'$R(r)$' 
    if WFType == '3DWF':
        WFFlag = 3
        WFName = 'Radial Wavefunction'#'$R(r)$' 
    

    #Set layout options common to all 3 wavefunctions
    layout = go.Layout(
        hovermode=False,
        dragmode="orbit",
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
                backgroundcolor='rgb(230, 230,230)'
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
                backgroundcolor='rgb(230, 230,230)'
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
                backgroundcolor='rgb(230, 230,230)'
            ),aspectratio=dict(
                x=1,
                y=1,
                z=1
            )
        )
    )


    #Plot radial wavefunction or radial distribution function
    if WFFlag == 2 or WFFlag == 1:
        traces = []
        if len(Orbitals) > 0:
            for n in range(0,len(Orbitals)):
                id = Orbitals[n]
                traces.append(go.Scatter(x = np.linspace(lowerxlim,upperxlim,1000), y = functiondict[id](np.linspace(lowerxlim,upperxlim,1000), WFFlag, zeff), line = dict(width = Thickness), name = id, hoverinfo = 'none'))
    
        return {'data': traces,
            'layout': go.Layout(
                xaxis = {
                    'autorange' : True,
                    'showgrid' : xgrid,
                    'zeroline' : False,
                    'showline' : True,
                    'title' : {'text' : 'Distance / a.u.', 'font' : {'size' : TextSize} },
                    'ticks' :'outside',
                    'tickfont' : {'size' : TextSize},
                    'showticklabels' : True
                },
                yaxis = {
                    'autorange' :True,
                    'showgrid' :ygrid,
                    'zeroline' :False,
                    'fixedrange' : True,
                    'title' : {'text' : WFName, 'font' : {'size' : TextSize} },
                    'showline' :True,
                    'ticks' :'outside',
                    'tickfont' : {'size' : TextSize},
                    'showticklabels' :True
                    },
                legend = {
                    'font' : {'size' : TextSize}
                }
            )
            }

    #Plot wavefunction isosurface
    if WFFlag == 3 :
        if len(Orbitals) > 0:
            return  {'data': OrbCalc(Orbitals[0]),'layout': layout}
        else:
            return {'data' : None, 'layout' : layout}


#Callback which defines what changes ('config') and what causes the change 
#Callback for save settings box
@app.callback(ddep.Output('RDF_Graph', 'config'),[ ddep.Input('PlotFormatDropdown', 'value'), ddep.Input('PlotHeightInput', 'value'),ddep.Input('PlotWidthInput', 'value')])

#Function which is called after the element is pressed
def SaveOptions(PlotFormat, PlotHeight, PlotWidth):


    return {
      'toImageButtonOptions': {
        'format': PlotFormat, 
        'filename': 'Radial Distribution Functions',
        'height': PlotHeight,
        'width': PlotWidth,
      }, 'modeBarButtonsToRemove' :['sendDataToCloud',
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
        'hoverCompareCartesian'],
        'displaylogo': False
    }


if __name__ == '__main__':
    app.run_server(debug=True)
