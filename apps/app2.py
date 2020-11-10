import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash.dependencies as ddep
import dash_defer_js_import as dji
import plotly.graph_objs as go
import plotly.colors as pc
import numpy as np
import os

from app import app

navbar = dbc.NavbarSimple(
    id = "navbar",
    children=[
        dbc.NavItem(dbc.NavLink(id= "orb_tab", children = "Orbitals", href="/apps/app1")),
        dbc.NavItem(dbc.NavLink(id= "vib_tab", children = "Vibrations", href="/apps/app2",active=True)),
        dbc.NavItem(dbc.NavLink(id= "trans_tab", children = "Translations", href="/apps/app3")),
    ],
    brand="Waveplot",
    brand_href="/apps/app1",
    color="primary",
    dark=True,
)

vib_graph = dcc.Graph(
                    id='vib_graph', 
                    style = {
                        'responsive' : 'true',
                        'height' : '580px',
                        'automargin' : 'true'
                    }
                )

vib_options = [

    dbc.Row(
        [
            dbc.Col(
                children = [
                    dbc.InputGroup(
                        [
                            dbc.InputGroupAddon("k", addon_type="prepend"),
                            dbc.Input(id="fc_input", placeholder=100, value = 100, type = 'number', min = 0),
                            dbc.InputGroupAddon(r"N m⁻¹", addon_type="append"),
                        ]
                    )
                ],
                style={"height" : "100%"},
                className="h-100"
            ),
            dbc.Col(
                children = [
                    dbc.Row([
                        dbc.InputGroup(

                            [
                                dbc.DropdownMenu(
                                    children = [
                                        dbc.DropdownMenuItem("ω", id = "omega"),
                                        dbc.DropdownMenuItem("𝜈", id = "nubar", style={"text-decoration": "overline"})
                                    ],
                                    id="freq_type",
                                    label="ω",
                                    addon_type="prepend"
                                ),
                                dbc.Input(id="freq_input", placeholder=100, type = 'number'),
                                dbc.InputGroupAddon(r"cm⁻¹", addon_type="append"),
                            ]
                        )
                    ]),
                    ],
                style={"height" : "100%"},
                className="h-100"
            ),
        ],
        className="h-100",
    ),


]

vib_body = dbc.Container(
    dbc.Row(
        [
            dbc.Col(vib_graph, style={"height" : "100%"},className="h-100"),
            dbc.Col(
                html.Div(
                    style={
                        "height": "400px",
                        "width": "400px",
                        "position" : "relative"},
                    className='viewer_3Dmoljs',
                    **{
                        'data-pdb': '2POR',
                        'data-backgroundcolor':'0xffffff',
                        'data-style':'stick'
                    }
                    ,
                    
                )
            ),
            dbc.Col(vib_options, style={"height" : "100%"},className="h-100"),
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
        vib_body,
        footer,
    ]
)

@app.callback(
    [
        ddep.Output('vib_graph', 'figure'),
        ddep.Output('freq_type', 'label')
    ], 
    [
        ddep.Input("freq_input", "value"),
        ddep.Input("fc_input", "value"),
        ddep.Input("omega", "n_clicks"),
        ddep.Input("nubar", "n_clicks"),
    ]
)

def update_app(freq_input, fc_input, *args):
    """
    Updates the app, given the current state of the UI
    All inputs correspond (in the same order) to the list 
    of ddep.Input a few lines above^^^

    Input:
        vib_checklist (list, string)   :: names of orbitals

    Returns:
        figure (dict)             :: dictionary item for list of go.XX objects
        config (dict)             :: dictionary item for go.layout object, e.g. linewidth, colormap...

    """
    ctx = dash.callback_context

    if not ctx.triggered:
        button_id = "all"
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if button_id == "omega":
        freq = freq_input * 2*np.pi
        button_label = ["ω"]
    elif button_id == "nubar":
        freq = freq_input
        button_label = ["𝜈"]
    else:
        button_label = ["𝜈"]

    x = np.linspace(-0.5*10**-10,0.5*10**-10,1000)

    print(type(fc_input), flush=True)

    if fc_input == None:
        fc = 0
    else:
        fc = fc_input

    data = []
    data.append(
        go.Scatter(
            x = x*10**10,
            y = (0.5*fc*x**2)/(1.98630 * 10**-23),
            line = dict(width = 2),
            name = 'Plot',
            hoverinfo = 'none',
        )
    )

    layout = go.Layout(
                xaxis = {
                    'autorange' : True,
                    'showgrid'  : False,
                    'zeroline'  : False,
                    'showline'  : True,
                    'title' : {
                        'text' : "Displacement (" + u"\u212B" + ")"
                    },
                    'ticks' :'outside',
                    'showticklabels' : True
                },
                yaxis = {
                    'autorange'  : True,
                    'showgrid'   : False,
                    'zeroline'   : False,
                    'title' :{
                        'text' : "Energy (cm⁻¹)",
                    },
                    'title_standoff' : 100,
                    'showline' :True,
                    'ticks' :'outside',
                    'showticklabels' :True
                },
                margin=dict(l=90, r=30, t=30, b=60),
    )

    output = {
        "data" : data,
        "layout" : layout
    }

    return output, [button_label]
