#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 17:11:36 2024

@author: jitong2023
@description: dash app
"""

import math
import dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from dash import html
from dash import dcc
import dash_bootstrap_components as dbc
import dash_daq as daq
import dash_bio
import textwrap

from layout_helper import run_standalone_app
from plot_contactMap import predict_structure, plot_contactMap, plot_contactMapAlign, parseRNAfold_prob


def header_colors():
    return {
        'bg_color': '#3d85c6',
        'font_color': 'white'
    }
app_name = 'RNA ContactMap'


initial_sequences = {
    'PDB_01019': {
        'sequence': 'AUGGGCCCGGGCCCAAUGGGCCCGGGCCCA',
        'structure': '.((((((())))))).((((((()))))))',
        'options': {
            'applyForce': True,
            'name': 'PDB_01019'
        }
    }
}


def description():
    return 'RNA secondary structure analysis.'


def layout():
    return html.Div(
        id='forna-body',
        className='app-body',
        children=[
            html.Div(
                id='forna-control-tabs',
                className='control-tabs',
                children=[
                    dcc.Tabs(id='forna-tabs', value='what-is', children=[
                        dcc.Tab(
                            label='About',
                            value='what-is',
                            children=html.Div(className='control-tab', children=[
                                html.H4(className='what-is', children='RNA ContactMap'),
                                dcc.Markdown('''
                                The RNA ContactMap Visualizer offers a novel method for visualizing RNA structures. \n
                                In the "Add Sequence" tab, you can enter a sequence by specifying the nucleotide sequence. \n
                                In the "Information" tab, you can: \n
                                Select which template sequence will be displayed by RNA contact maps, which link bases with colors 
                                indicating base pairing probability. \n
                                Select which query sequence to be aligned with the template and specify the alignment position 
                                by the slider. \n
                                Select sequences to be displayed by Forna, as well as obtain information about the sequences 
                                that you have already created.
                                ''')
                            ])
                        ),

                        dcc.Tab(
                            label='Add Sequence',
                            value='add-sequence',
                            children=html.Div(className='control-tab', children=[
                                html.Div(
                                    title='Enter the first nucleotide sequence.',
                                    className='app-controls-block',
                                    children=[
                                        html.Div(className='fullwidth-app-controls-name', children='Sequence'),
                                        html.Div(
                                            className='app-controls-desc',
                                            children='Specify the template nucleotide sequence.'
                                        ),
                                        dcc.Textarea(
                                            id='new-sequence',
                                            placeholder=initial_sequences['PDB_01019']['sequence'],
                                            rows=2,
                                            style = {'marginLeft':'10px', 'width': '90%'}
                                        ),
                                        html.Br(),
                                        html.Br(),
                                        
                                        html.Div(className='fullwidth-app-controls-name', children='ID'),
                                        html.Div(
                                            className='app-controls-desc',
                                            children='Specify a unique ID for this sequence.'
                                        ),
                                        dcc.Input(
                                            id='new-id', placeholder='PDB_01019',
                                            style = {'marginLeft':'10px', 'width': '90%'}
                                        ),
                                        html.Br(),
                                        html.Br(),
                                    
                                        html.Div(id='new-error-message'),
                                        html.Button(id='new-submit', children='Submit sequence', n_clicks=0, style = {'marginLeft':'10px'})
                                    ]
                                )
                                
                            ])
                        ),
                        dcc.Tab(
                            label='Information',
                            value='show-sequences',
                            children=html.Div(className='control-tab', children=[

                                html.Div(
                                    className='app-controls-block',
                                    children=[
                                        html.Div(
                                            className='fullwidth-app-controls-name',
                                            children='Sequences to display in ContactMap'
                                        ),
                                        html.Br(),
                                        html.Div(
                                            className='app-controls-desc',
                                            children='Choose the template sequence by ID.'
                                        ),
                                        dcc.Dropdown(id='template-id'),
                                        
                                        html.Div(
                                            className='app-controls-desc',
                                            children='Choose the query sequence by ID.'
                                        ),
                                        dcc.Dropdown(id='query-id'),
                                        html.Br(),
                                        html.Div(className='fullwidth-app-controls-name', children='Align position'),
                                        html.Div(
                                            className='app-controls-desc',
                                            children='Specify alignment position for query on template.'
                                        ),
                                        dcc.Slider(0, 1, step=0.005, marks={ 0: "5'", 1:"3'"}, value=0, id='align-pos'),
                                        html.Div(id='align-error-message', style = {'marginLeft':'10px'})
                                        
                                    ]
                                ),
                                html.Hr(),
                                html.Div(
                                    className='app-controls-block',
                                    children=[
                                        html.Div(
                                            className='fullwidth-app-controls-name',
                                            children='Sequences to display in Forna'
                                        ),
                                        html.Div(
                                            className='app-controls-desc',
                                            children='Choose the sequences to display by ID.'
                                        ),
                                        dcc.Dropdown(
                                            id='forna-sequences-display',
                                            clearable=True,
                                            multi=True
                                        )
                                    ]
                                ),
                                html.Hr(),
                                html.Div(
                                    className='app-controls-block',
                                    children=[
                                        html.Div(
                                            className='app-controls-block',
                                            children=[
                                                html.Div(
                                                    className='fullwidth-app-controls-name',
                                                    children='Sequence information by ID'
                                                ),
                                                html.Div(
                                                    className='app-controls-desc',
                                                    children='Search for a sequence by ID to get more information.'
                                                ),
                                                dcc.Dropdown(
                                                    id='forna-sequences-info-search'
                                                ),
                                                html.Br(),
                                                html.Div(id='forna-sequence-info')
                                            ]
                                        )
                                    ]
                                )
                                

                            ])
                        )  
                    ])
                ]),
            html.Div(id='forna-display-tabs', className='display-tabs', children = [
                dcc.Tabs(children=[
                    dcc.Tab(label='ContactMap Template', children=[
                        html.Div(id='templatemap-container', style={'display': 'inline-block', 'width': '60%', 'marginTop':'40px', 'padding': '10px'}, 
                            #children=[dcc.Graph(id = 'template-map')]),
                            children=[html.Img(id = 'template-map', src = '')])
                    ]),
                    dcc.Tab(label='ContactMap Alignment', children=[
                        html.Div(id='alignmap-container', style={'display': 'inline-block', 'width':'60%', 'marginTop':'40px', 'padding': '10px'}, 
                            #children=[dcc.Graph(id = 'align-map')]), 
                            children=[html.Img(id = 'align-map', src = '')])
                    ]),
                    dcc.Tab(label='Forna Visualizer', children=[
                        html.Div(id='forna-container', children=[
                            dash_bio.FornaContainer(
                                id='forna',
                                height=500,
                                width=800
                            )
                        ])
                    ])
                ])
            ]),
            dcc.Store(id='contactmap-sequences', data={})

        ]
    )


def callbacks(_app):

    @_app.callback(
        [Output('contactmap-sequences', 'data'),
         Output('new-error-message', 'children')],
        [Input('new-submit', 'n_clicks')],
        [State('new-sequence', 'value'),
         State('new-id', 'value'),
         State('contactmap-sequences', 'data')]
    )
    def add_sequence(nclicks, new_sequence, new_seqid, database):

        if nclicks == 0:
            raise PreventUpdate

        if new_sequence is None:
            raise PreventUpdate

        new_sequence = new_sequence.upper().replace('T', 'U')
        if not all(base in ['A', 'U', 'G', 'C'] for base in new_sequence):
            new_error_msg = html.P(
                'Sequence {} contains undefined characters. '.format(new_seqid) +
                'Please make sure the input sequence only contains ACGTU.',
                style={'color': 'red'}
            )
            return database, new_error_msg

                
        if new_seqid not in database.keys():
            [seq, structure, scores, probRecord] = predict_structure(new_sequence, new_seqid)
            new_error_msg = html.P(
                'Successfully added {}!'.format(new_seqid),
                style={'color': 'green'}
            )
            database[new_seqid] = [seq, structure, '', probRecord]
        else:
            new_error_msg = html.P(
                'You already have a sequence with this ID. ' +
                'Please choose a different ID, or check the next tab ' +
                'to see which IDs have already been taken.',
                style={'color': 'red'}
            )
        
        return database, new_error_msg


    @_app.callback(
        [Output('forna-sequences-display', 'options'),
         Output('forna-sequences-info-search', 'options'),
         Output('template-id', 'options'),
         Output('query-id', 'options'),
         ],
        [Input('contactmap-sequences', 'data')]
    )
    def update_sequences(data):

        if data is None:
            raise PreventUpdate

        new_options = [
            {'label': sequence_id,
             'value': sequence_id}
            for sequence_id in data.keys()
        ]

        return new_options, new_options, new_options, new_options

    @_app.callback(
        Output('forna-sequence-info', 'children'),
        [Input('forna-sequences-info-search', 'value')],
        [State('contactmap-sequences', 'data')]
    )
    def update_sequence_info(sequence_id, data):
        if data is None or sequence_id is None:
            raise PreventUpdate

        return html.Div(
            [
                'Sequence:\n',
                textwrap.fill(data[sequence_id][0], width=30),
                html.Br(),
                'Structure:\n',
                textwrap.fill(data[sequence_id][1], width=30)
            ]
        )

    @_app.callback(
        Output('forna', 'sequences'),
        [Input('forna-sequences-display', 'value')],
        [State('contactmap-sequences', 'data')]
    )
    def update_shown_sequences(selected_sequence_ids, stored_sequences):

        if selected_sequence_ids is None or stored_sequences is None:
            raise PreventUpdate

        sequences = []

        for sequence_id in selected_sequence_ids:
            stored_sequence_info = {
                'sequence': stored_sequences[sequence_id][0],
                'structure': stored_sequences[sequence_id][1],
                'options': {'name': sequence_id}
            }
            sequences.append(stored_sequence_info)

        return sequences

    @_app.callback(
        #Output('template-map', 'figure'),
        Output('template-map', 'src'),
        [Input('template-id', 'value')],
        [State('contactmap-sequences', 'data')]
    )
    def update_template_contactMap(template_id, database):
        if template_id is None:
            raise PreventUpdate
        seqinfo = database[template_id]
        fig = plot_contactMap(seqinfo)
        return fig

    @_app.callback(
        #Output('align-map', 'figure'),
        [Output('align-map', 'src'),
         Output('align-error-message', 'children')],
        [Input('template-id', 'value'),
         Input('query-id', 'value'),
         Input('align-pos', 'value')],
        [State('contactmap-sequences', 'data')]
    )
    def update_align_contactMap(template_id, query_id, align_pos, database):
        if template_id is None or query_id is None:
            raise PreventUpdate
        seqinfo = database[template_id]
        query_seqinfo = database[query_id]
        align_pos = math.ceil((len(seqinfo[0]) + len(query_seqinfo[0])) * align_pos) - len(query_seqinfo[0])
        align_pos, fig = plot_contactMapAlign(seqinfo, query_seqinfo, align_pos)

        if align_pos<=0:
            align_pos -= 1

        error_msg = html.P(
            'Query 5\' aligned to Template position {}'.format(align_pos),
            style={'color': 'green'}
        )
        return fig, error_msg

app = run_standalone_app(layout, callbacks, header_colors, app_name)
server = app.server

if __name__ == '__main__':
    app.run_server(debug=True, port=8050)

