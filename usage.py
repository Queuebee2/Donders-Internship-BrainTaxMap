
import dash
from dash.dependencies import Input, Output
import dash_html_components as html
import dash_core_components as dcc
import dash_table
import pandas as pd
from dash_network import Network
import dash_cytoscape as cyto
import plotly.graph_objects as go

from dateutil import parser

from braintaxmap.Neo4jManager import Neo4jManager, recordsToJSON
from braintaxmap.tools import fuzz, FUZZING
from braintaxmap.querymachine import QueryMachine
FUZZING=True

GRAPH_START_TERM = 'barrel cortex'
import json
import os
import logging
ERROR_LOG_FILE = 'braintaxmap' + os.sep +'logs' + os.sep+ 'dash-network.log'

# create logger
logger = logging.getLogger('dash-network-usage.main')
logger.setLevel(logging.INFO)
# create file handler which logs even debug messages
fh = logging.FileHandler(ERROR_LOG_FILE)
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)



# for testing... 
db = Neo4jManager()
qm = QueryMachine()

app = dash.Dash(__name__)

app.scripts.config.serve_locally = True
app.css.config.serve_locally = True



selected_colors = ['#006', '#060', '#600', '#A80', '#A08']



def net_data(selected, range_value, verbose=False): 
    if selected:
        if verbose: print('plotting ', selected[:25], 'slider:', range_value, end = '')
    selected_letter = selected and selected[0]
    
    """we need to rebuild something like this
    to represent it with the d3j network plotly graph 
    GRAPH:
     {'nodes': [
         {'id': 'work'}, 
         {'id': 'in'}, 
         {'id': 'progress'}
         ], 
     'links': [
        {'source': 'work', 'target': 'in'}, 
        {'source': 'in', 'target': 'progress'}, 
        {'source': 'progress', 'target': 'progress'}
        ]}
    
    """

    # TODO https://dash.plotly.com/sharing-data-between-callbacks

    amounts = 510
    if verbose: print(amounts)
    # records = db.find_anywhere('soma')
    records = db.websiteTester()

    nodes = list()
    links = list()

    if verbose: print('parsing nodes and links....')

    def create_id(node):
        if 'PMID' in node.keys():
            return node['PMID'], 1
        else:
            return node['name'], 0


    for index, record in enumerate(records):
        if index > range_value:
            break

        startnode = record['n']
        relations = record['r']
        targetnode = record['m']

        startnode, color1 = create_id(startnode)
        targetnode, color2 = create_id(targetnode)

        if {'id':startnode} not in nodes:nodes.append({'id':startnode}) 
        if {'id':targetnode} not in nodes:nodes.append({'id':targetnode})
        if {'source': startnode, 'target': targetnode} not in links: links.append({'source': startnode, 'target': targetnode})


    # for i, rec in enumerate(records):
        
    #     try:
    #         n1 = rec['start']['name']
    #         n2 = rec['target']['name']
    #     except Exception:
    #         continue

    #     # TODO MAKE THIS WAY WAY WAY WAY MORE EFFICIENT PLEASE
    #     if {'id':n1} not in nodes:nodes.append({'id':n1}) 
    #     if {'id':n2} not in nodes:nodes.append({'id':n2})
    #     if {'source': n1, 'target': n2} not in links: links.append({'source': n1, 'target': n2})

    #     if round( (i / amounts) * 10) % round(amounts/10) == 0:
    #         print('.',end='')

    if verbose: print('done!')

    return {
        'nodes': nodes,
        'links': links
    }

def cyto_data(selected):
    records = db.websiteTester()

    elements = []

    for start, relation, target in records:
        n1 = start['name']
        n2 = target['name']
        elements.append({'data': {'id': n1, 'label': n1}})
        elements.append({'data': {'id': n2, 'label': n2}})
        elements.append({'data': {'source': n1, 'target': n2}})


    return elements
    

#############################################################################
#############################################################################
#############################################################################

default_stylesheet = [
    {
        'selector': 'node',
        'style': {
            'background-color': '#BFD7B5',
            'label': 'data(label)'
        }
    }
]

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}

ALLOWED_TYPES = (
    "text", "number", "password", "email", "search",
    "tel", "url", "range", "hidden",
)

app.layout = html.Div([
    html.H2('Click a node to expand it, or the background to return'),

#############################################################################
    # network div
#############################################################################
    html.Div([
    # slider to decide amount of nodes available, see update_data callback
    dcc.Slider(
        id='input_amountOfNodesSlider',
        min=1,
        max=250,
        marks={10:'10',50:'50',100:'100'},
        step=1,
        value=15,
    ),
    # network (within loading)
    dcc.Loading(
        id="loading-network-dash",
        type="default",
        children=Network(
            id='net',
            data=net_data(f'{GRAPH_START_TERM}', 100)
            )
        ),
    # output information on network
    html.Div(id='output', children=[])
    ],
    style={'width': '50%',}
    ),

#############################################################################
#############################################################################
#############################################################################
    # playing around with CYTOSCAPE
    html.Div(
        id='cytoscape-container',
        children=[
        dcc.Dropdown(
            id='dropdown-update-layout-cyto',
            value='circle',
            clearable=False,
            options=[
                {'label': name.capitalize(), 'value': name}
                for name in ['grid', 'random', 'circle', 'cose', 'concentric']
        ]
        ),
        dcc.Loading(
            id="loading-network-cyto",
            type="default",
            children=cyto.Cytoscape(
                id='cytoscape-network',
                layout={'name': 'circle'},
                stylesheet=default_stylesheet,
                style={'width': '75%', 'height': '400px'},
                elements=cyto_data(f'{GRAPH_START_TERM}')
            )
        ),
        # the thing you click on
        dcc.Loading(
            id="loading-cyto-selected",
            type="default",
            children=html.Pre(id='cyto-selected-json', style=styles['pre'])
        ),
        ]),
    html.Div(
    [   # random inputfields
        dcc.Input(
            id="input_{}".format(_),
            type=_,
            placeholder="input type {}".format(_),
        )
        for _ in ALLOWED_TYPES
    ]
    + [html.H2('testing output stuff:')] +
     # output of random inputfields
     [html.Br(), html.Div(id="out-all-types",children=[])]
),

#############################################################################
    # playing around more div
#############################################################################
    html.Div([
        html.Div(id='hidden-article-storage', style={'display':'none'}),
        dcc.Input(id='query-input', type='text', value='barrel cortex',
         placeholder='type a query here', debounce=True),
        dcc.Dropdown(
            id='dropdown-select-result',  multi=True), 
        dcc.Slider(
            id='max-results-slider',
            min=1,
            max=1000000,
            marks={1:'1',500000:'500000',1000000:'1000000'},
            step=1,
            value=15,
        ),
        html.Div(id='article-amount-slider-text', children=[]),
        dcc.Loading(
            id="loading-query-results",
            type="default",
            children=[
                html.Div(id='query-results',children=[]),
                dash_table.DataTable( id='output-table'),
                html.Div(id='dropdown-selected-articles', children=[], style={'width': '75%'})
                ]
        ),
        
        ],style={'width': '50%',})
])

"""
    
    'dropdown-select-result' -> 'value' -> callback -> update something ->
    requires staets?
    'dropdown-select-result'
"""



@app.callback(Output('cyto-selected-json', 'children'),
              Input('cytoscape-network', 'tapNodeData'))
def displayTapNodeData(data):
    return json.dumps(data, indent=2)


@app.callback(Output('cytoscape-network', 'layout'),
              Input('dropdown-update-layout-cyto', 'value'))
def update_layout(layout, verbose=False):
    if verbose: print('chosen layout:', layout)
    return {
        'name': layout,
        'animate': True
    }


@app.callback(
Output('article-amount-slider-text', 'children'),
Input('max-results-slider','value'))
def show_amount_selected_articles(amount):
    return dcc.Markdown(f'max amount of articles to retrieve: **{amount}**')

@app.callback(Output('dropdown-selected-articles', 'children'),
                Output('hidden-article-storage', 'style'),
              Input('dropdown-select-result', 'value'),
              Input('hidden-article-storage', 'children'))
def present_selected_articles(selected, jsonified_article_data):
    if not selected:
        selected = [1]
    else:
        for s in selected:
            print(s)

    article_data = json.loads(jsonified_article_data)
    for article in article_data:
        if str(article['PMID']) in selected:
            print('this one is selected:', article['TI'])

    return [html.Div(id=f'article_{i}', children=[dcc.Markdown(f"# {i}. {r['TI']}   \n  - this is an article test -   \n {r['AB']}")]) for i, r in enumerate(article_data) if r['PMID'] in selected ], {'display':'none'}

###################################################################################
###################################################################################
###################################################################################
# using the last inputfield, do a query and show results
@app.callback([Output("query-results", "children"),
                Output("dropdown-select-result", "options"),
                Output("output-table", "columns"),
                Output("output-table", "data"),
                Output('hidden-article-storage', 'children')],
        [Input('query-input', "value"),
        Input('max-results-slider', 'value')])
def query_from_input(query, max_results):
    
    idlist = qm.search_pubmed_ids(query)

    articles = qm.get_pubmed_by_pmids(idlist, max_results)

    meshTerms = dict()
    words = dict()
    noMesh=[]
    noMeshAndBefore1970=[]
    noAb=[]
    first50 = []
    for i, record in enumerate(articles):
        if i >= max_results:
            break
        if i < 50:
            first50.append(record)
        if 'MH' not in record.keys():
            noMesh.append(record)
            publishyear = int(record['DP'].split(' ')[0])
            logger.info(f" {record['PMID']} 'publishdate':{publishyear}")
            if  publishyear < 1970:
                noMeshAndBefore1970.append(record)
        else:
            for term in record['MH']:
                try:
                    meshTerms[term] +=1
                except:
                    meshTerms[term] =1
        
        if 'AB' not in record.keys():
            noAb.append(record)
            continue
        else:
            for word in record['AB'].split(' '):
                try:
                    words[word] +=1
                except:
                    words[word] = 1
    
    words_counts = sorted(list(words.items()), key=lambda e: e[1],reverse=True)
    headers = [{'name':x, "id": x} for x in ['word', 'count']]
    df = pd.DataFrame(words_counts[:100], columns=['word', 'count'])
    md = dcc.Markdown(
        f"""
- you looked for _{query}_ and limited the amount of results to {max_results}
- there are {len(idlist)} results  
- {len(noMesh)} articles had no MeSH terms
- {len(noMeshAndBefore1970)} articles without MeSH terms were published before 1970  
- there are {len(words)} unique words  
- there are {len(meshTerms)} unique MeSH terms  
# word counts
"""
)
    options = [{'label':r['TI'], 'value':str(r['PMID'])} for (i,r) in enumerate(first50)]

    jsonified_articles = json.dumps(first50)

        
    return md, options, [{"name":i, "id":i} for i in df.columns], df.to_dict('records'), jsonified_articles

@app.callback(
    Output("out-all-types", component_property="children"),
    [Input("input_{}".format(_), component_property="value") for _ in ALLOWED_TYPES]

) 
def cb_render(*vals, verbose=False):
    if verbose: print('the values in cb_render are:', vals)
    return " | ".join((str(val) for val in vals if val))


@app.callback(Output('net', 'data'),
              [Input('net', 'selectedId'),
              Input("input_amountOfNodesSlider","value")])
def update_data(selected_id, range_value):
    return net_data(selected_id, range_value)

@app.callback(Output('output', 'children'),
              [Input('net', 'selectedId'), Input('net', 'data')])
def list_connections(selected_id, data):
    return dcc.Markdown('# You selected node __"`{}`"__ on a graph with `{}` nodes and `{}` links'.format(
        selected_id, len(data['nodes']), len(data['links'])))
        


if __name__ == '__main__':
    app.run_server(debug=True)
