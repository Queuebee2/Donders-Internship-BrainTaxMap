
import logging
import os
import json
import dash
import plotly.express as px 
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
from dash_html_components.Button import Button
import dash_table
import pandas as pd
from dash_network import Network
import dash_cytoscape as cyto
import plotly.graph_objects as go
import squarify

from dateutil import parser

from braintaxmap.Neo4jManager import Neo4jManager, recordsToJSON
from braintaxmap.tools import fuzz, FUZZING
from braintaxmap.querymachine import QueryMachine
FUZZING = True

GRAPH_START_TERM = 'barrel cortex'
MAX_TESTING_AMOUNT = 500

ERROR_LOG_FILE = 'braintaxmap' + os.sep + 'logs' + os.sep + 'dash-network.log'


# Logging
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
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

# STOPWORDS
# For now hardcoded. Could also do:
# import nltk
# nltk.download() # before doing import stopwords
# from nltk.corpus import stopwords
# and try printing the words using
# stopwords.words('english')
STOPWORDS_EN = ['i', 'me', 'my', 'myself', 'we', 'our', 'ours', 'ourselves', 'you', "you're", "you've", "you'll", "you'd", 'your', 'yours', 'yourself', 'yourselves', 'he', 'him', 'his', 'himself', 'she', "she's", 'her', 'hers', 'herself', 'it', "it's", 'its', 'itself', 'they', 'them', 'their', 'theirs', 'themselves', 'what', 'which', 'who', 'whom', 'this', 'that', "that'll", 'these', 'those', 'am', 'is', 'are', 'was', 'were', 'be', 'been', 'being', 'have', 'has', 'had', 'having', 'do', 'does', 'did', 'doing', 'a', 'an', 'the', 'and', 'but', 'if', 'or', 'because', 'as', 'until', 'while', 'of', 'at', 'by', 'for', 'with', 'about', 'against', 'between', 'into', 'through', 'during', 'before', 'after', 'above', 'below', 'to', 'from', 'up',
                'down', 'in', 'out', 'on', 'off', 'over', 'under', 'again', 'further', 'then', 'once', 'here', 'there', 'when', 'where', 'why', 'how', 'all', 'any', 'both', 'each', 'few', 'more', 'most', 'other', 'some', 'such', 'no', 'nor', 'not', 'only', 'own', 'same', 'so', 'than', 'too', 'very', 's', 't', 'can', 'will', 'just', 'don', "don't", 'should', "should've", 'now', 'd', 'll', 'm', 'o', 're', 've', 'y', 'ain', 'aren', "aren't", 'couldn', "couldn't", 'didn', "didn't", 'doesn', "doesn't", 'hadn', "hadn't", 'hasn', "hasn't", 'haven', "haven't", 'isn', "isn't", 'ma', 'mightn', "mightn't", 'mustn', "mustn't", 'needn', "needn't", 'shan', "shan't", 'shouldn', "shouldn't", 'wasn', "wasn't", 'weren', "weren't", 'won', "won't", 'wouldn', "wouldn't"]

# append other stuff later. STOPWORDS will hold all filtering words for now
STOPWORDS = STOPWORDS_EN

"""Todo
        - this is a possible feature to the dash-cytoscape or biopython or both libraries if we can make something that automatically puts records to be ready for the thingy
        - restructure project again
        - figure out import issues 
        - the todos (find them)
        - look at networkX and try Neo4j -> (magic) -> networkX -> somevisualiser (dash-cyto/dash-network)
        - this would be insanely useful aswell (collapse nodes into one node) https://github.com/iVis-at-Bilkent/cytoscape.js-expand-collapse
        - build again from https://github.com/ned2/slapdash
        - CAN IT ACTUALLY SUPPORT NETWORKX? AS OF 2019 IT COULDNT YET
    
    Loose ideas
        # ADD PUBMED TO GRAPH - BUTTON
        - live pubmed search, user can add entries to the graph. Clicking them would start expanding the nodes, maybe even something that searches for a shortest path between two loose nodes by
        looking at articles but that would be very computation heavy, so we have to provide the pre-made graph db via an api or somtehing.
        that means even if the first search is live on pubmed,  [(maybe show whether the entry exists in the larger graph db)] , we check the db for it and start expanding from there (cheaper??? because
        nodes already have edges to other nodes. How would the user define custom settings? should we do EVERYTHING on the database and let user select which ones are relevant, or let user
        build a custom one )

    HARD ISSUES
        can't have multiple callbacks edit the same output
        https://community.plotly.com/t/duplicate-callback-outputs-issue/43017
"""

# for testing...
db = Neo4jManager()
qm = QueryMachine()

app = dash.Dash(__name__)

app.scripts.config.serve_locally = True
app.css.config.serve_locally = True


selected_colors = ['#006', '#060', '#600', '#A80', '#A08']


def net_data(selected, range_value, verbose=False):
    """ from a query (not impl) and a max amount, query the database (testing) and return the graph"""
    if selected:
        if verbose:
            print('plotting ', selected[:25], 'slider:', range_value, end='')
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
    if verbose:
        print(amounts)
    # records = db.find_anywhere('soma')
    records = db.websiteTester(range_value)

    nodes = list()
    links = list()

    if verbose:
        print('parsing nodes and links....')

    def create_id(node):
        if 'PMID' in node.keys():
            return node['PMID'], 1
        else:
            return node['name'], 0

    for index, (startnode, relation, targetnode) in enumerate(records):
        if index > range_value:
            break
        
        startnode, color1 = create_id(startnode)
        targetnode, color2 = create_id(targetnode)

        if {'id': startnode} not in nodes:
            nodes.append({'id': startnode})
        if {'id': targetnode} not in nodes:
            nodes.append({'id': targetnode})
        if {'source': startnode, 'target': targetnode} not in links:
            links.append({'source': startnode, 'target': targetnode})

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

    if verbose:
        print('done!')

    return {
        'nodes': nodes,
        'links': links
    }


def cyto_data(selected, amount):
    """ from a query (not impl) and a max amount, query the database (testing) and return the graph"""
    # limit max amount to testing amount
    amount = MAX_TESTING_AMOUNT if MAX_TESTING_AMOUNT < amount else amount

    records = db.websiteTester(amount)

    elements = []

    # TODO ; functions to check what kind of node is coming in

    for start, relation, target in records:
        # look again
        # elems = [ {'data': {'id': , 'label' .... ?}}, {'data': {'source': n1[ID], 'target': n2[ID!!@!@!@!@]}}]
        n1 = start['PMID'] if 'PMID' in start.keys() else start['name']
        n2 = target['PMID'] if 'PMID' in target.keys() else target['name']

        elements.append({'data': {'id': n1, 'label': n1}})
        elements.append({'data': {'id': n2, 'label': n2}})
        elements.append({'data': {'source': n1, 'target': n2}})

    print(f'found {len(elements)} nodes for <websiteTester>')
    return elements

def treemap():
    
    # x = 0.
    # y = 0.
    # width = 100.
    # height = 100.
    
    # values = [500, 433, 78, 25, 25, 7]

    # normed = squarify.normalize_sizes(values, width, height)
    # rects = squarify.squarify(normed, x, y, width, height)

    # # Choose colors from http://colorbrewer2.org/ under "Export"
    # color_brewer = ['rgb(166,206,227)','rgb(31,120,180)','rgb(178,223,138)',
    #                 'rgb(51,160,44)','rgb(251,154,153)','rgb(227,26,28)']
    # shapes = []
    # annotations = []
    # counter = 0

    # for r in rects:
    #     shapes.append( 
    #         dict(
    #             type = 'rect', 
    #             x0 = r['x'], 
    #             y0 = r['y'], 
    #             x1 = r['x']+r['dx'], 
    #             y1 = r['y']+r['dy'],
    #             line = dict( width = 2 ),
    #             fillcolor = color_brewer[counter]
    #         ) 
    #     )
    #     annotations.append(
    #         dict(
    #             x = r['x']+(r['dx']/2),
    #             y = r['y']+(r['dy']/2),
    #             text = values[counter],
    #             showarrow = False
    #         )
    #     )
    #     counter = counter + 1
    #     if counter >= len(color_brewer):
    #         counter = 0

    # figure = {
    # 'data': [go.Scatter(
    #     x = [ r['x']+(r['dx']/2) for r in rects ], 
    #     y = [ r['y']+(r['dy']/2) for r in rects ],
    #     text = [ str(v) for v in values ], 
    #     mode = 'text',
    #     )
    # ],
    # 'layout': go.Layout(
    #     height=700, 
    #     width=700,
    #     xaxis={'showgrid':False, 'zeroline':False, 'showticklabels': False},
    #     yaxis={'showgrid':False, 'zeroline':False, 'showticklabels': False},
    #     shapes=shapes,
    #     annotations=annotations,
    #     hovermode='closest',
    #     )
    # }
    # return figure

    df = pd.DataFrame(
        {'searchword':{x:'barrel cortex' for x in range(10)},
         'word':{x:f'word-{x}' for x in range(10)},
         'count':{x:x*x for x in range(10)}
        }
)
    figure = px.treemap(
       df, path=['searchword', 'word'], values='count'
    )
    return figure



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
            marks={10: '10', 50: '50', 100: '100'},
            step=1,
            value=15,
        ),
        # network (within loading)
        dcc.Loading(
            id="loading-network-dash",
            type="default",
            children=Network(
                # Allowed arguments: data, dataVersion, height, id, linkWidth,
                # maxLinkWidth, maxRadius, nodeRadius, selectedId, width
                id='net',
                data=net_data(f'{GRAPH_START_TERM}', 100)

            )
        ),
        # output information on network
        html.Div(id='dash-selected-node-output', children=[])
    ],

        style={'width': '75%', 'height': '700px', 'border': '3px solid grey'}
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
                ],
                style={'width': '75%'}),
            dcc.Input(id="searchterm-cyto-input", type='text', value=f'{GRAPH_START_TERM}',
                      placeholder='type a query here', debounce=True),
            html.Button('test add elem', id='btn-add-test', n_clicks=0),
            html.Button('test remove elem', id='btn-remove-test', n_clicks=0),
            dcc.Loading(
                id="loading-network-cyto",
                type="default",
                children=cyto.Cytoscape(
                    id='cytoscape-network',
                    layout={'name': 'circle'},
                    stylesheet=default_stylesheet,
                    style={'width': '75%', 'height': '700px',
                           'border': '3px solid grey'},
                    elements=[]
                )
            ),
            # Information on the node you selected
            dcc.Loading(
                id="loading-cyto-selected",
                type="default",
                children=html.Pre(id='cyto-selected-json', style=styles['pre'])
            ),
        ], style={'border': '3px solid grey'}),
    html.Div(
        [   # random inputfields
            dcc.Input(
                id="input_{}".format(_),
                type=_,
                placeholder="input type {}".format(_),
            )
            for _ in ALLOWED_TYPES
        ]
        + [html.P('type in the testboxes to let values appear down here')] +
        # output of random inputfields
        [html.Br(), html.Div(id="out-all-types", children=[])],
        style={'border': '2px solid grey'}),

    #############################################################################
    # playing around more div
    #############################################################################
    html.Div([
        # treemap
        dcc.Loading(id='loading-treemap-results', type='circle',
        children=[dcc.Graph(id='treemap-container', figure=treemap()),
                html.Button('print treemap', id='btn-treemap', n_clicks=0),
                html.Div(id='treemap-selected-node-output', children=[])
                ],
                style={'width': '75%'}
        ),
        html.Div(id='hidden-article-storage', style={'display': 'none'}),

        # todo load from article
        dcc.Markdown(
            ['# Perform a search', 'Here you can query the pubmed database']),
        dcc.Input(id='query-input', type='text', value='barrel cortex',
                  placeholder='type a query here', debounce=True),
        dcc.Checklist(id='query-checklist',
                      options=[
                          {'label': 'Filter Stopwords',
                              'value': 'filter-stopwords'},
                          {'label': 'placeholder', 'value': 'placeholder'}
                      ]),
        dcc.Dropdown(
            id='dropdown-select-result',  multi=True),
        dcc.Slider(
            id='max-results-slider',
            min=1,
            max=1000000,
            marks={1: '1', 500000: '500000', 1000000: '1000000'},
            step=1,
            value=150,
        ),
        html.Div(id='article-amount-slider-text', children=[]),
        dcc.Loading(
            id="loading-query-results",
            type="default",
            children=[
                html.Div(id='query-results', children=[]),
                dash_table.DataTable(id='output-table',
                                     # todo add fancy options
                                     # todo change height (its taking way too much space. limit rows to ~10?)
                                     style_table={"selector": ".dash-spreadsheet-inner",
                                                  "rule": 'max-height: "calc(100vh - 226px)"'}
                                     ),
                html.Div(id='dropdown-selected-articles',
                         children=[], style={'width': '75%'})
            ]
        ),

    ], style={'width': '50%'})
], style={'width': '50%'})


"""
notes
The dcc.Graph component has four attributes that can change through user-interaction: 
 hoverData, clickData, selectedData, relayoutData. 
 These properties update when you hover over points, 
 click on points, or select regions of points in a graph.
"""
@app.callback(Output('treemap-selected-node-output', 'children'),
        Input('treemap-container', 'hoverData'))
def treemap_selected(treemap_figure):
    # https://community.plotly.com/t/treemap-click-events/31337
    # Input( .. plotly_treemapclick?
    print(treemap_figure)
    print(dir(treemap_figure))
    return dcc.Markdown(str(treemap_figure))


"""
    
    'dropdown-select-result' -> 'value' -> callback -> update something ->
    requires staets?
    'dropdown-select-result'
"""



@app.callback(
    Output('cytoscape-network', 'elements'),
    Input('btn-add-test', 'n_clicks'),
    Input('btn-remove-test', 'n_clicks'),
    Input('searchterm-cyto-input', 'value'),
    Input('query-input', "value"),
    Input('max-results-slider', 'value'),
    State('cytoscape-network', 'elements'))
def update_cytograph(btn_add, btn_remove, query1, query2, amount_val, elements):
    """
    TODO 
        - cant have multiple callbacks output to same element
          https://community.plotly.com/t/duplicate-callback-outputs-issue/43017

        - cant have same input and output (how will we add someting to the chart??)
          https://community.plotly.com/t/same-input-and-output-in-a-callback/20002

        - should be renamed to update_cyto something

        updates the cytoscape graph according to the settings chosen
    """
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]

    if 'btn-add-test' in changed_id:
        secretvalue = str(len(elements)**3 if len(elements) < 12 else 42)
        try:
            if len(elements) == 0:
                elements += [{'data': {'id': f'1',
                                       'label': f'1', 'secretvalue': secretvalue}}]
            else:
                prevnode = int(elements[-1]['data']['id'])
                elements += [{'data': {'source': f'{prevnode}',
                                       'target': f'{prevnode+1}'}}]
                elements += [{'data': {'id': f'{prevnode+1}',
                                       'label': f'{prevnode+1}', 'secretvalue': secretvalue}}]
        except ValueError:
            'stop testing this now please'

    elif 'btn-remove-test' in changed_id:
        elements = elements[:-2]
    elif 'max-results-slider' in changed_id:
        elements = elements
    elif ('query-input' in changed_id) or 'searchterm-cyto-input' in changed_id:
        query = query1 or query2
        if query == None:
            return [{'id': 'do a search', 'label': 'do a query first'}]
        else:
            elements = cyto_data(query, amount_val)
    else:
        pass

    return elements


@app.callback(Output('cyto-selected-json', 'children'),
              Input('cytoscape-network', 'tapNodeData'))
def displayTapNodeData(data):
    # TODO ; functions to check what kind of node is coming in,
    # then markdown it and display it niceley.
    # Article nodes can show all the information on the article
    # other nodes might be short for now
    return json.dumps(data, indent=2)


@app.callback(Output('cytoscape-network', 'layout'),
              Input('dropdown-update-layout-cyto', 'value'))
def update_layout(layout, verbose=False):
    if verbose:
        print('chosen layout:', layout)
    return {
        'name': layout,
        'animate': True
    }


@app.callback(
    Output('article-amount-slider-text', 'children'),
    Input('max-results-slider', 'value'))
def show_amount_selected_articles(amount):
    return dcc.Markdown(f'max amount of articles to retrieve: **{amount}**')


@app.callback(Output('dropdown-selected-articles', 'children'),
              Input('dropdown-select-result', 'value'),
              Input('hidden-article-storage', 'children'))
def present_selected_articles(selected, jsonified_article_data):
    """ when articles are found, the dropdown is populated with their titles
        all chosen titles are to be displayed
        articles come from storage in the hidden-article-storage div(50 articles max now)
    """
    if not selected:
        selected = [1]
    else:
        for s in selected:
            print(s)

    article_data = json.loads(jsonified_article_data)
    for article in article_data:
        if str(article['PMID']) in selected:
            print('this one is selected:', article['TI'])

    markup_articles = []
    for i, r in enumerate(article_data):

        # this article was selected to be displayed
        if r['PMID'] in selected:
            meshterms = r['MH'] if 'MH' in r.keys() else 'No mesh Terms'
            markup_articles.append(
                html.Div(
                    id=f'article_display{i}',
                    children=[
                        dcc.Markdown([
                            # TITLE
                            f"# {i}. {r['TI']} PMID:{r['PMID']}",
                            f"##### this is an article test",
                            f"##### Abstract",
                            f"{r['AB']}",
                            f"##### MeSH",
                            f"{meshterms}",
                            f"",
                            f"",
                            f"",
                            f"",
                        ])
                    ]
                )
            )

    return markup_articles


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
               Input('max-results-slider', 'value'),
               Input('query-checklist', 'value')])
def query_from_input(query, max_results, checklist):
    # TODO just move all under here to a seperate function please (or 5 ...)

    filter_stopwords = False  # todo -> move out, can be done after searching.
    if checklist:  # apparently parameters can also be None
        for check_val in checklist:
            if check_val == 'filter-stopwords':
                filter_stopwords = True

    idlist = qm.search_pubmed_ids(query)

    articles = qm.get_pubmed_by_pmids(idlist, max_results)

    meshTerms = dict()
    words = dict()
    noMesh = []
    noMeshAndBefore1970 = []
    noAb = []
    first50 = []
    for i, record in enumerate(articles):
        if i >= max_results:
            break
        if i < 100:
            first50.append(record)

        if 'MH' not in record.keys():
            noMesh.append(record)
            publishyear = int(record['DP'].split(' ')[0])
            if publishyear < 1970:
                logger.info(f" {record['PMID']} 'publishdate':{publishyear}")
                noMeshAndBefore1970.append(record)
        else:
            for term in record['MH']:
                try:
                    meshTerms[term] += 1
                except:
                    meshTerms[term] = 1

        if 'AB' not in record.keys():
            noAb.append(record)
            continue
        else:
            for w in record['AB'].split(' '):
                # cant just lower everything, you wont be able to see genes and stuff!
                word = w.lower()
                if filter_stopwords and word in STOPWORDS:
                    continue
                try:
                    words[word] += 1
                except:
                    words[word] = 1

    words_counts = sorted(list(words.items()),
                          key=lambda e: e[1], reverse=True)
    headers = [{'name': x, "id": x} for x in ['word', 'count']]
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
    options = [{'label': r['TI'], 'value':str(
        r['PMID'])} for (i, r) in enumerate(first50)]

    # at the moment we only dump 50 (top) articles to show it's quite expensive to put them all in..
    # should have a temp storage in db but that will break with multiple users (this is cached in client browser)
    jsonified_articles = json.dumps(first50)

    return md, options, [{"name": i, "id": i} for i in df.columns], df.to_dict('records'), jsonified_articles


@app.callback(
    Output("out-all-types", component_property="children"),
    [Input("input_{}".format(_), component_property="value")
     for _ in ALLOWED_TYPES]

)
def cb_render(*vals, verbose=False):
    if verbose:
        print('the values in cb_render are:', vals)
    return " | ".join((str(val) for val in vals if val))


@app.callback(Output('net', 'data'),
              [Input('net', 'selectedId'),
               Input("input_amountOfNodesSlider", "value")])
def update_data(selected_id, range_value):
    return net_data(selected_id, range_value)


@app.callback(Output('dash-selected-node-output', 'children'),
              [Input('net', 'selectedId'), Input('net', 'data')])
def list_connections(selected_id, data):
    return dcc.Markdown('# You selected node __"`{}`"__ on a graph with `{}` nodes and `{}` links'.format(
        selected_id, len(data['nodes']), len(data['links'])))


if __name__ == '__main__':
    app.run_server(debug=True)
