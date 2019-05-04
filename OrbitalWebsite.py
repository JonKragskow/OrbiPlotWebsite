import numpy as np
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import dash.dependencies as ddep
import dash_daq as daq
from subprocess import call
import json

#####################################################################################  Experimental  #####################################################################################

#####################################################################################
################################ Orbital Mathematics ################################
#####################################################################################


############################p Orbital Wavefunction functions#########################
def PrincipleAxisFuncp(n,r,f):
    out = n*np.exp(r/n)*0.003*np.sqrt(np.pi/3)/abs(f)
    return out

def fp(Orb,r):
    if Orb=='2p': out = 1/(2*np.sqrt(6))
 
    if Orb=='3p': out = (4 - 2*r/3)/(9*np.sqrt(6))

    if Orb=='4p': out = (20 - 5*r + (r/2)**2)/(32*np.sqrt(15))

    if Orb=='5p': out = (120 - 90*(2*r/5) + 18*(2*r/5)**2 - (2*r/5)**3)/(150*np.sqrt(30))

    if Orb=='6p': out = (840 - 840*r/3 + 252*(r/3)**2 - 28*(r/3)**3 + (r/3)**4)/(432*np.sqrt(210))
    return out

############################d Orbital Wavefunction functions#########################

def fd(Orb,r): 
    if Orb=='3d': out = 1/(9*np.sqrt(30))

    if Orb=='4d': out = (6-r/2)/(96*np.sqrt(5))

    if Orb=='5d': out = (42 - 28*r/5 + (2*r/5)**2)/(150*np.sqrt(70))

    if Orb=='6d': out = (336 - 168*(r/3) + 24*(r/3)**2 - (r/3)**3)/(864*np.sqrt(105))
    return out

############################f Orbital Wavefunction functions#########################

#def ff(Orb,r): 
#    if Orb=='4f': out = 
#
#    if Orb=='5f': out = 
#
#    if Orb=='6f': out = 
#    return out



def OrbDomain(Orb):

    DomainConstants=[]


#####2p surface bounds#####
    if Orb == '2p':
        p2LowerBound1 = 0.0305424
        p2UpperMinusLowerBound1 = 11.973154-0.03054255
       
        DomainConstants = [p2LowerBound1, p2UpperMinusLowerBound1]

#####3p surface bounds#####
    if Orb == '3p':

        p3LowerBound1 = 0.0521075737439532
        p3UpperMinusLowerBound1 =  5.5936882283563
        
        p3LowerBound2 = 6.40190923834268
        p3UpperMinusLowerBound2 =  14.3297971165078

        DomainConstants = [p3LowerBound1, p3UpperMinusLowerBound1, p3LowerBound2, p3UpperMinusLowerBound2]
    
#####4p surface bounds#####

    if Orb=='4p':
    
        p4LowerBound1 = 0.0791783
        p4UpperMinusLowerBound1 = 4.9947867

        p4LowerBound2 = 6.072543
        p4UpperMinusLowerBound2 = 6.791527

        p4LowerBound3 = 16.566121
        p4UpperMinusLowerBound3 = 13.745981

        DomainConstants = [p4LowerBound1, p4UpperMinusLowerBound1, p4LowerBound2, p4UpperMinusLowerBound2, p4LowerBound3, p4UpperMinusLowerBound3]
    
#####5p surface bounds#####

    if Orb=='5p':
    
        p5LowerBound1 = 0.1111277 
        p5UpperMinusLowerBound1 = 4.6587693
        
        p5LowerBound2 = 6.1052591
        p5UpperMinusLowerBound2 = 5.3146969
        
        p5LowerBound3 = 16.046327
        p5UpperMinusLowerBound3 = 6.007799789
        
        p5LowerBound4 = 37.13 #broken
        p5UpperMinusLowerBound4 = 26.9003 #broken

        DomainConstants = [p5LowerBound1, p5UpperMinusLowerBound1, p5LowerBound2, p5UpperMinusLowerBound2, p5LowerBound3, p5UpperMinusLowerBound3, p5LowerBound4, p5UpperMinusLowerBound4]
    
#####6p surface bounds#####
    
    if Orb=='6p':

        p6LowerBound1 = 0.011494042
        p6UpperMinusLowerBound1 = 5.227
        
        p6LowerBound2 = 5.35
        p6UpperMinusLowerBound2 = 7.5063
    
        p6LowerBound3 = 10.5611
        p6UpperMinusLowerBound3 = 13.0505
        
        p6LowerBound4 = 16.6383
        p6UpperMinusLowerBound4 = 24.7636
        
        p6LowerBound5 = 41.5770
        p6UpperMinusLowerBound5 = 43.0776
    
        DomainConstants = [p6LowerBound1, p6UpperMinusLowerBound1, p6LowerBound2, p6UpperMinusLowerBound2, p6LowerBound3, p6UpperMinusLowerBound3, p6LowerBound4, p6UpperMinusLowerBound4, p6LowerBound5, p6UpperMinusLowerBound5]

#####3d surface bounds#####
#    if Orb == '3d':
#
#        d3LowerBound1 = 
#        d3UpperMinusLowerBound1 =  
#       
#        DomainConstants = [d3LowerBound1, d3UpperMinusLowerBound1]
#    
######4d surface bounds#####
#
#    if Orb=='4d':
#    
#        d4LowerBound1 = 
#        d4UpperMinusLowerBound1 = 
#
#        d4LowerBound2 = 
#        d4UpperMinusLowerBound2 = 
#
#        DomainConstants = [d4LowerBound1, d4UpperMinusLowerBound1, d4LowerBound2, d4UpperMinusLowerBound2]
#    
######5d surface bounds#####
#
#    if Orb=='5d':
#    
#        d5LowerBound1 =
#        d5UpperMinusLowerBound1 = 
#        
#        d5LowerBound2 = 
#        d5UpperMinusLowerBound2 = 
#        
#        d5LowerBound3 = 
#        d5UpperMinusLowerBound3 = 
#        
#        DomainConstants = [d5LowerBound1, d5UpperMinusLowerBound1, d5LowerBound2, d5UpperMinusLowerBound2, d5LowerBound3, d5UpperMinusLowerBound3]
#    
######6d surface bounds#####
#    
#    if Orb=='6d':
#
#        d6LowerBound1 = 
#        d6UpperMinusLowerBound1 = 
#        
#        d6LowerBound2 = 
#        d6UpperMinusLowerBound2 = 
#    
#        d6LowerBound3 = 
#        d6UpperMinusLowerBound3 = 
#        
#        d6LowerBound4 = 
#        d6UpperMinusLowerBound4 = 
#    
#        DomainConstants = [d6LowerBound1, d6UpperMinusLowerBound1, d6LowerBound2, d6UpperMinusLowerBound2, d6LowerBound3, d6UpperMinusLowerBound3, d6LowerBound4, d6UpperMinusLowerBound4]
#
######4f surface bounds#####
#
#    if Orb=='4f':
#    
#        f4LowerBound1 = 
#        f4UpperMinusLowerBound1 = 
#
#        DomainConstants = [f4LowerBound1, f4UpperMinusLowerBound1]
#    
######5f surface bounds#####
#
#    if Orb=='5f':
#    
#        f5LowerBound1 = 
#        f5UpperMinusLowerBound1 =
#        
#        f5LowerBound2 = 
#        f5UpperMinusLowerBound2 = 
#        
#        f5LowerBound3 =
#        f5UpperMinusLowerBound3 = 
#        
#        DomainConstants = [f5LowerBound1, f5UpperMinusLowerBound1, f5LowerBound2, f5UpperMinusLowerBound2]
#    
######6f surface bounds#####
#    
#    if Orb=='6f':
#
#        f6LowerBound1 = 
#        f6UpperMinusLowerBound1 = 
#        
#        f6LowerBound2 = 
#        f6UpperMinusLowerBound2 = 
#    
#        f6LowerBound3 = 
#        f6UpperMinusLowerBound3 = 
#    
#        DomainConstants = [f6LowerBound1, f6UpperMinusLowerBound1, f6LowerBound2, f6UpperMinusLowerBound2, f6LowerBound3, f6UpperMinusLowerBound3]
#
    return DomainConstants

############################Orbital Wavefunction Calculation#########################


def OrbCalc(Orb):

    letter = Orb.strip('Orb')
    n = float(letter.strip('p'))

    #Arrays for parametric equations
    v = np.linspace(0, 1, 70)
    u = np.linspace(0, 2.*np.pi, 70)

    x = np.zeros([np.size(u),np.size(v)])
    y = np.zeros([np.size(u),np.size(v)])

    zp = np.zeros([np.size(u),np.size(v)])
    zm = np.zeros([np.size(u),np.size(v)])

    data = []

#Plot different radial portions of orbital

    t=0
    flag = True

    while t <= (2*n-4):

		#n     = 2 3 4 5 6
        #2n    = 4 6 8 10 12
        #2n-4  = 0 2 4 6 8
		#lobes = 2 4 6 8 10
		#t     = 0 2 4 6 8
		
        for vele,vval in enumerate(v):
            for uele,uval in enumerate(u):
                x[vele,uele] = np.sqrt((OrbDomain(Orb)[t+1]*vval + OrbDomain(Orb)[t])**2 - PrincipleAxisFuncp(n,OrbDomain(Orb)[t+1]*vval+OrbDomain(Orb)[t],fp(Orb,OrbDomain(Orb)[t+1]*vval+OrbDomain(Orb)[t]))**2)*np.cos(uval)
                y[vele,uele] = np.sqrt((OrbDomain(Orb)[t+1]*vval + OrbDomain(Orb)[t])**2 - PrincipleAxisFuncp(n,OrbDomain(Orb)[t+1]*vval+OrbDomain(Orb)[t],fp(Orb,OrbDomain(Orb)[t+1]*vval+OrbDomain(Orb)[t]))**2)*np.sin(uval)
                
                zp[vele,uele] = PrincipleAxisFuncp(n,OrbDomain(Orb)[t+1]*vval + OrbDomain(Orb)[t],fp(Orb,OrbDomain(Orb)[t+1]*vval+OrbDomain(Orb)[t]));
                zm[vele,uele] = -zp[vele,uele];
        
        t += 2
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
	O1s = 2.*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O1s**2.
	if mode==2:
		return O1s

def s_2(r,mode, zeff):
	rho = 2.*r/2. * zeff
	O2s = 1/(2*np.sqrt(2))*(2-rho)*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O2s**2.
	if mode==2:
		return O2s

def s_3(r,mode, zeff):
	rho = 2.*r/3. * zeff
	O3s = 1/(9*np.sqrt(3))*(6-6*rho+rho**2)*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O3s**2.
	if mode==2:
		return O3s

def s_4(r,mode, zeff):
	rho = 2.*r/4. * zeff
	O4s = (1./96)*(24-36*rho+12*rho**2-rho**3)*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O4s**2.
	if mode==2:
		return O4s

def s_5(r,mode, zeff):
	rho = 2.*r/5. * zeff
	O5s = (1/(300*np.sqrt(5)))*(120-240*rho+120*rho**2-20*rho**3+rho**4)*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O5s**2.
	if mode==2:
		return O5s

def s_6(r,mode, zeff):
	rho = 2.*r/6. * zeff
	O6s = (1/(2160*np.sqrt(6)))*(720-1800*rho+1200*rho**2-300*rho**3+30*rho**4-rho**5)*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O6s**2.
	if mode==2:
		return O6s

def p_2(r,mode, zeff):
	rho = 2.*r/2. * zeff
	O2p = 1/(2*np.sqrt(6))*rho*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O2p**2.
	if mode==2:
		return O2p

def p_3(r,mode, zeff):
	rho = 2.*r/3. * zeff
	O3p = 1/(9*np.sqrt(6))*rho*(4-rho)*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O3p**2.
	if mode==2:
		return O3p

def p_4(r,mode, zeff):
	rho = 2.*r/4. * zeff
	O4p = 1/(32*np.sqrt(15))*rho*(20-10*rho+rho**2)*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O4p**2.
	if mode==2:
		return O4p

def p_5(r,mode, zeff):
	rho = 2.*r/5. * zeff
	O5p = 1/(150*np.sqrt(30))*rho*(120-90*rho+18*rho**2-rho**3)*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O5p**2.
	if mode==2:
		return O5p

def p_6(r,mode, zeff):
	rho = 2.*r/6. * zeff
	O6p = 1/(432*np.sqrt(210))*rho*(840-840*rho+252*rho**2-28*rho**3+rho**4)*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O6p**2.
	if mode==2:
		return O6p

def d_3(r,mode, zeff):
	rho = 2.*r/3. * zeff
	O3d = 1/(9*np.sqrt(30))*rho**2*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O3d**2.
	if mode==2:
		return O3d

def d_4(r,mode, zeff):
	rho = 2.*r/4. * zeff
	O4d = 1/(96*np.sqrt(5))*(6-rho)*rho**2*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O4d**2.
	if mode==2:
		return O4d

def d_5(r,mode, zeff):
	rho = 2.*r/5. * zeff
	O5d = 1/(150*np.sqrt(70))*(42-14*rho+rho**2)*rho**2*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O5d**2.
	if mode==2:
		return O5d

def d_6(r,mode, zeff):
	rho = 2.*r/6. * zeff
	O6d = 1/(864*np.sqrt(105))*(336-168*rho+24*rho**2-rho**3)*rho**2*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O6d**2.
	if mode==2:
		return O6d

def f_4(r,mode, zeff):
	rho = 2.*r/4. * zeff
	O4f = 1/(96*np.sqrt(35))*rho**3*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O4f**2.
	if mode==2:
		return O4f

def f_5(r,mode, zeff):
	rho = 2.*r/5. * zeff
	O5f = 1/(300*np.sqrt(70))*(8-rho)*rho**3*np.exp(-rho/2) * zeff**1.5
	if mode==1:
		return r**2.*O5f**2.
	if mode==2:
		return O5f

def f_6(r,mode, zeff):
	rho = 2.*r/5. * zeff
	O6f = 1/(2592*np.sqrt(35))*(rho**2-18*rho+72)*rho**3*np.exp(-rho/2) * zeff**1.5
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
						values=['5p'],
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
						value='3DWF',
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
						values=[],
						),

					dcc.Checklist(id = 'ygridcheck',
						options=[
							{'label': ' y-Axis ', 'value': 'ygridval'}
						],
						values=[],
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
@app.callback(ddep.Output('RDF_Graph', 'figure'),[ddep.Input('OrbCheck', 'values'),ddep.Input('WFType', 'value'), ddep.Input('LineWidthSlider','value'), ddep.Input('TextSizeSlider','value'), ddep.Input('xgridcheck','values'), ddep.Input('ygridcheck','values'), ddep.Input('ZeffSlider','value'),ddep.Input('UpperxInput', 'value'),ddep.Input('LowerxInput', 'value')])

#Function which is called after the element is pressed
def UpdatePlot(Orbitals, WFType, Thickness, TextSize, xgridinput, ygridinput, zeff, upperxlim, lowerxlim):


	if xgridinput == ['xgridval']:
		xgrid = True
	else:
		xgrid = False

	if ygridinput == ['ygridval']:
		ygrid = True
	else:
		ygrid = False

	if WFType == 'RDF':
		WFFlag = 1
		WFName = 'Radial Distribution Function' #' 4\pi r R(r)'
	if WFType == 'RWF':
		WFFlag = 2
		WFName = 'Radial Wavefunction'#'$R(r)$' 
	if WFType == '3DWF':
		WFFlag = 3
		WFName = 'Radial Wavefunction'#'$R(r)$' 
	

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

	traces = []

	if len(Orbitals) > 0:
		for n in range(0,len(Orbitals)):
			id = Orbitals[n]
			traces.append(go.Scatter(x = np.linspace(lowerxlim,upperxlim,1000), y = functiondict[id](np.linspace(lowerxlim,upperxlim,1000), WFFlag, zeff), line = dict(width = Thickness), name = id, hoverinfo = 'none'))

	if WFFlag == 2 or WFFlag == 1:
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
	if WFFlag == 3 :
		return  {'data': OrbCalc(Orbitals[0]),'layout': layout}


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
