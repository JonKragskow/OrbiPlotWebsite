import numpy as np
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import dash.dependencies as ddep
import dash_daq as daq
from subprocess import call
import json

#To do list
#undo button remove



###Functions for function


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
						values=['1s','2p','3d','4f'],
						labelStyle={'maxwidth' : '20px'}
					),

				]),

				html.Div(className = 'WF box', style = {'width' : '47%' , 'padding-left' : '2%' ,'float' : 'right', 'display' : 'inline-block'}, children = [

					html.H2(children = 'Orbital Options'),

					dcc.RadioItems(id = 'WFType',
						options=[
							{'label': 'Radial Distribution Function', 'value': 'RDF'},
							{'label': 'Radial Wave Function', 'value': 'RWF'}
						],
						value='RDF',
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
