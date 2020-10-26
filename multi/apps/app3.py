import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash.dependencies as ddep
import dash_defer_js_import as dji
import plotly.graph_objs as go
import numpy as np
import os

from app import app

navbar = dbc.NavbarSimple(
    id = "navbar",
    children=[
        dbc.NavItem(dbc.NavLink(id= "orb_tab", children = "Orbitals", href="/apps/app1")),
        dbc.NavItem(dbc.NavLink(id= "vib_tab", children = "Vibrations", href="/apps/app2")),
        dbc.NavItem(dbc.NavLink(id= "trans_tab", children = "Translations", href="/apps/app3",active=True)),
    ],
    brand="Waveplot",
    brand_href="#",
    color="primary",
    dark=True,
)

layout = html.Div(children=navbar)