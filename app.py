#! /usr/bin/env python3

import dash

app = dash.Dash(__name__, 
    external_stylesheets=["https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css"],
    suppress_callback_exceptions=False,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1"}
    ])



app.scripts.append_script({
"external_url": 'https://3Dmol.org/build/3Dmol-min.js'
})

server = app.server

