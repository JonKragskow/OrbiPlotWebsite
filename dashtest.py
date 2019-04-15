import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import numpy as np


#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__)
app.css.append_css({"external_url": "https://max.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css"})

colors = {
    'background': '#666555',
    'Header': '#265564',
    'text': '#7FDBFF'
}


#####################################################################################
################################ Orbital Mathematics ################################
#####################################################################################


############################p Orbital Wavefunction functions#########################
def PrincipleAxisFuncp(n,r,f):
    out = n*np.exp(r/n)*0.003*np.sqrt(np.pi/3)/abs(f)
    return out

def fp(Orb,r):
    if Orb=='2p': out = 1/(2*np.sqrt(6))
 
    if Orb=='3p':out = (4 - 2*r/3)/(9*np.sqrt(6))

    if Orb=='4p':out = (20 - 5*r + (r/2)**2)/(32*np.sqrt(15))

    if Orb=='5p':out = (120 - 90*(2*r/5) + 18*(2*r/5)**2 - (2*r/5)**3)/(150*np.sqrt(30))

    if Orb=='6p':out = (840 - 840*r/3 + 252*(r/3)**2 - 28*(r/3)**3 + (r/3)**4)/(432*np.sqrt(210))
    return out

############################d Orbital Wavefunction functions#########################
def zFunc(n,r,f):
    out = n*np.exp(r/n)*0.003*np.sqrt(np.pi/3)/abs(f)
    return out

def f3d(r):
    out = 1/(9*np.sqrt(30))
    return out

def f4d(r):
    out = (6-r/2)/(96*np.sqrt(5))
    return out

def f5d(r):
    out = (42 - 28*r/5 + (2*r/5)**2)/(150*np.sqrt(70))
    return out

def f6d(r):
    out = (336 - 168*(r/3) + 24*(r/3)**2 - (r/3)**3)/(864*np.sqrt(105))
    return out

############################f Orbital Wavefunction functions#########################
#def zFunc(n,r,f):
#    out = n*np.exp(r/n)*0.003*np.sqrt(np.pi/3)/abs(f)
#    return out
#
#
#def f4f(r):
#    out = 
#    return out
#
#def f5f(r):
#    out = 
#    return out
#
#def f6f(r):
#    out = 
#    return out

def OrbDomain(Orb):

    DomainConstants=[]


#####2p surface bounds#####
    if Orb=='2p':
        p2LowerBound1 = 0.0305424
        p2UpperMinusLowerBound1 = 11.973154-0.03054255
       
        DomainConstants = [p2LowerBound1,p2UpperMinusLowerBound1]

#####3p surface bounds#####
    if Orb=='3p':

        p3LowerBound1 = 0.0521075737439532
        p3UpperMinusLowerBound1 =  5.5936882283563
        
        p3LowerBound2 = 6.40190923834268
        p3UpperMinusLowerBound2 =  14.3297971165078

        DomainConstants = [p3LowerBound1,p3UpperMinusLowerBound1,p3LowerBound2,p3UpperMinusLowerBound2]
    
#####4p surface bounds#####

    if Orb=='4p':
    
        p4LowerBound1 = 0.0791783
        p4UpperMinusLowerBound1 = 4.9947867
        p4LowerBound2 = 6.072543
        p4UpperMinusLowerBound2 = 6.791527
        p4LowerBound3 = 16.566121
        p4UpperMinusLowerBound3 = 13.745981

        DomainConstants = [p4LowerBound1,p4UpperMinusLowerBound1,p4LowerBound2,p4UpperMinusLowerBound2,p4LowerBound3,p4UpperMinusLowerBound3]
    
#####5p surface bounds#####

    if Orb=='5p':
    
        p5LowerBound1 = 0.1111277 
        p5UpperMinusLowerBound1 = 4.6587693
        
        p5LowerBound2 = 6.1052591
        p5UpperMinusLowerBound2 = 5.3146969
        
        p5LowerBound3 = 16.046327
        p5UpperMinusLowerBound3 = 21.13525524-16.046327
        
        p5LowerBound4 =  35.12342545
        p5UpperMinusLowerBound4 = 26

        DomainConstants = [p5LowerBound1,p5UpperMinusLowerBound1,p5LowerBound2,p5UpperMinusLowerBound2,p5LowerBound3,p5UpperMinusLowerBound3,p5LowerBound4,p5UpperMinusLowerBound4]
    
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
    
        DomainConstants = [p6LowerBound1,p6UpperMinusLowerBound1,p6LowerBound2,p6UpperMinusLowerBound2,p6LowerBound3,p6UpperMinusLowerBound3,p6LowerBound4,p6UpperMinusLowerBound4,p6LowerBound5,p6UpperMinusLowerBound5]

    return DomainConstants

############################Orbital Wavefunction Calculation#########################


def OrbCalc(Orb):

    letter = Orb.strip('Orb')
    n = float(letter.strip('p'))

    #Arrays for parametric equations
    v = np.linspace(0, 1, 50)
    u = np.linspace(0, 2.*np.pi, 50)
    x = np.zeros([np.size(u),np.size(v)])
    y = np.zeros([np.size(u),np.size(v)])
    zp = np.zeros([np.size(u),np.size(v)])
    zm = np.zeros([np.size(u),np.size(v)])

    data = []

#Plot first radial portion of p orbital

    t=0

    while t <= (2*n-4):

#n     = 2 3 4 5 6
#lobes = 2 4 5 6 10
#t     = 0 2 4 6 8

        for vele,vval in enumerate(v):
            for uele,uval in enumerate(u):
                x[vele,uele] = np.sqrt((OrbDomain(Orb)[t+1]*vval + OrbDomain(Orb)[t])**2 - PrincipleAxisFuncp(n,OrbDomain(Orb)[t+1]*vval+OrbDomain(Orb)[t],fp(Orb,OrbDomain(Orb)[t+1]*vval+OrbDomain(Orb)[t]))**2)*np.cos(uval)
                y[vele,uele] = np.sqrt((OrbDomain(Orb)[t+1]*vval + OrbDomain(Orb)[t])**2 - PrincipleAxisFuncp(n,OrbDomain(Orb)[t+1]*vval+OrbDomain(Orb)[t],fp(Orb,OrbDomain(Orb)[t+1]*vval+OrbDomain(Orb)[t]))**2)*np.sin(uval)
                zp[vele,uele] = PrincipleAxisFuncp(n,OrbDomain(Orb)[t+1]*vval + OrbDomain(Orb)[t],fp(Orb,OrbDomain(Orb)[t+1]*vval+OrbDomain(Orb)[t]));
                zm[vele,uele] = -PrincipleAxisFuncp(n,OrbDomain(Orb)[t+1]*vval + OrbDomain(Orb)[t],fp(Orb,OrbDomain(Orb)[t+1]*vval+OrbDomain(Orb)[t]));
        
        t = t+2
        data0 = go.Surface(x=x, y=y, z=zp,colorscale= [[0, 'rgb(109,0,157)'], [1,'rgb(109,0,157)']],showscale=False)
        data1 = go.Surface(x=x, y=y, z=zm,colorscale= [[0, 'rgb(255,204,51)'], [1,'rgb(255,204,51)']],showscale=False)
    
        data.append(data0)
        data.append(data1)
        

    return data



def p_3(mode):
    r = np.arange(1,100)
    rho = 2.*r/3.
    O3p = 1/(9*np.sqrt(6))*rho*(4-rho)*np.exp(-rho/2)
    if mode == 1:
        data =  go.Scatter(x = r, y = r**2.*O3p**2., mode = 'lines', name = 'lines')
    if mode == 2:
        data =  go.Scatter(x = r, y = O3p, mode = 'lines', name = 'lines')
    return data
    


#####################################################################################
####################################Plot Options#####################################
#####################################################################################

layout = go.Layout(
    title='Parametric Plot',
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







#####################################################################################
#########################################html########################################
#####################################################################################


app.layout = html.Div(className = "container",children=[
        html.Div(className = 'row',children=[
        html.Div(className = "col-lg-6",children=[
            dcc.Graph(id='Orbital1', figure={'data': p_3(2),'layout': layout})
        ])

])#,
#        html.Div(className = 'row',children=[
#        html.Div(className = "col-lg-6",children=[
#            dcc.Graph(id='Orbital3', figure={'data': OrbCalc('6p'),'layout': layout})
#        ]),
#        html.Div(className = "col-lg-6",children=[
#            dcc.Graph(id='Orbital4', figure={'data': OrbCalc('6p'),'layout': layout})
#        ])    
#
#])#,
#        html.Div(className = 'row',children=[
#        html.Div(className = "col-lg-6",children=[
#            dcc.Graph(id='Orbital5', figure={'data': OrbCalc('3p'),'layout': layout})
#        ]),
#        html.Div(className = "col-lg-6",children=[
#            dcc.Graph(id='Orbital6', figure={'data': OrbCalc('3p'),'layout': layout})
#        ])    
#
#])
])


if __name__ == '__main__':
    app.run_server(debug=True, threaded=True)