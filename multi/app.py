#! /usr/bin/env python3

import dash

app = dash.Dash(__name__, 
    external_stylesheets=["https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css"],
    external_scripts=['https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML'],
    suppress_callback_exceptions=True)

server = app.server

