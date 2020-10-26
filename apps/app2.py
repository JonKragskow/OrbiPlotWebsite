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

vib_tab = [

    html.Div(className = "item", 
        style = {
            'grid-column-start': '2',
            'grid-column-end': '3',
            'grid-row-start': '1',
            'grid-row-end': '2'
        },
        children=[
            dcc.Checklist(
                id = 'vib_checklist',
                style = {
                    'textAlign' : 'center', 
                },
                options=[
                    {"label": ' Plot ', "value": 'yes'}
                ],
                value=[],
            )        
        ]
    )
]




layout = html.Div(children=[

navbar,

html.Title(
    children='Waveplot'
),

html.Div(
        className = "container", 
        style = {
            'display' : 'grid',
            'grid-template-columns': r'50% 50%',
            'grid-template-rows' : r'100%',
            'margin': 'auto',
            'width' : '95%'
        },
        children=[
            html.Div(
                className = "item", 
                style = {
                    'grid-column-start': '1',
                    'grid-column-end': '2',
                    'grid-row-start': '1',
                    'grid-row-end': '2',
                    'justify-self': 'stretch',
                    'align-self': 'stretch'
                },
                children=[
                    dcc.Graph(
                        id='vib_plot_area', 
                        style = {
                            'responsive' : 'true',
                            'height' : '580px',
                            'automargin' : 'true'
                        }
                    )
                ]
            ),
            html.Div(children=vib_tab),
        ]
    ),
html.Footer(
    style = {
        'textAlign'       : 'center', 
        'font-size'       : 'smaller',
        'color'           : 'white',
        'background-color': '#3977AF'
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
),
])




@app.callback(
        [
        ddep.Output('vib_plot_area', 'figure')
        ], 
        [
        ddep.Input('vib_checklist', "value")
        ]
    )

def update_app(vib_checklist):
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

    return [vib_fig(vib_checklist)]

def vib_fig(vib_checklist):


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

    cols = def_cols + tol_cols + wong_cols

    data = []

    curr_ymax = 0.
    x = np.linspace(-20,20,1000)

    data.append(
        go.Scatter(
            x = x,
            y = x**2,
            line = dict(width = 2),
            name = 'Plot',
            hoverinfo = 'none',
            marker={"color":cols[0]}
        )
    )

    layout = go.Layout(
                xaxis = {
                    'autorange' : True,
                    'showgrid'  : False,
                    'zeroline'  : False,
                    'showline'  : True,
                    'title' : {
                        'text' : "Displacement",
                        #"text" : 'Distance (a₀)',
                        'font' :{'size' : 10} 
                    },
                    'ticks' :'outside',
                    'tickfont' :{'size' : 10},
                    'showticklabels' : True
                },
                yaxis = {
                    'autorange'  : True,
                    'showgrid'   : False,
                    'zeroline'   : False,
                    'fixedrange' : True,
                    'title' :{
                        'text' : 'Energy',
                        'font' :{
                            'size' : 10
                        }
                    },
                    'title_standoff' : 100,
                    'showline' :True,
                    'ticks' :'outside',
                    'tickfont' :{'size' : 10},
                    'showticklabels' :True
                },
                margin=dict(l=90, r=30, t=30, b=60),
    )

    output = {
        "data" : data,
        "layout" : layout
    }

    return output
