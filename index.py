import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from app import app
from apps import app1, app2, app3

url_bar_and_content_div = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])

app.layout = url_bar_and_content_div

app.validation_layout = html.Div([
    url_bar_and_content_div,
    app1.layout,
    app2.layout,
    app3.layout,
])

@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    print(pathname, flush=True)
    if pathname == '/apps/app1':
        return app1.layout
    elif pathname == '/apps/app2':
        return app2.layout
    elif pathname == '/apps/app3':
        return app3.layout
    else:
        return app1.layout

if __name__ == '__main__':
    app.run_server(debug=True)