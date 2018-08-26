import plotly.graph_objs as go
from plotly.offline import offline
from plotly import tools

cpk_colors = dict(Ar='cyan', B='salmon', Ba='darkgreen', Be='darkgreen', Br='darkred', C='black', Ca='darkgreen',
                  Cl='green', Cs='violet', F='green', Fe='darkorange', Fr='violet', H='white', He='cyan',
                  I='darkviolet', K='violet', Kr='cyan', Li='violet', Mg='darkgreen', N='blue', Na='violet', Ne='cyan',
                  O='red', P='orange', Ra='darkgreen', Rb='violet', S='yellow', Sr='darkgreen', Ti='gray', Xe='cyan')
cpk_color_rest = 'pink'


def gen_trace(adjacency_list: dict, elements: list, x_coordinates: list,
         y_coordinates: list, z_coordinates: list) -> None:
    """Creates atom and bond traces for the plot"""
    def gen_atom_trace():
        """Creates an atom trace for the plot"""
        ids = [f'{el} ({id_})' for id_, el in enumerate(elements)]
        colors = [cpk_colors.get(element, cpk_color_rest) for element in elements]
        markers = dict(color=colors, line=dict(color='lightgray', width=2), size=7, symbol='circle', opacity=0.8)
        trace = go.Scatter3d(x=x_coordinates, y=y_coordinates, z=z_coordinates, mode='markers', marker=markers,
                             text=ids, textposition='top center')
        return trace

    def gen_bond_trace():
        """"Creates a bond trace for the plot"""
        trace = go.Scatter3d(x=[], y=[], z=[], hoverinfo='none', mode='lines',
                             marker=dict(color='grey', size=7, opacity=1))
        adjascent_atoms = ((atom, neighbour) for atom, neighbours in adjacency_list.items()
                           for neighbour in neighbours)
        for i, j in adjascent_atoms:
            trace['x'] += (x_coordinates[i], x_coordinates[j], None)
            trace['y'] += (y_coordinates[i], y_coordinates[j], None)
            trace['z'] += (z_coordinates[i], z_coordinates[j], None)
        return trace

    return gen_atom_trace(), gen_bond_trace()


def plot(trace1, trace2, annotations1 = None, annotations2 = None, plot_name = 'plot'):
    """Creates two 3D scatter plots"""
    if not annotations1:
        annotations1 = []
    if not annotations2:
        annotations2 = []
    fig = tools.make_subplots(rows=1, cols=2, specs=[[{'is_3d': True}, {'is_3d': True}]], print_grid=False,
                              horizontal_spacing=0.01, )
    fig.append_trace(trace1[0], 1, 1)
    fig.append_trace(trace1[1], 1, 1)
    fig.append_trace(trace2[0], 1, 2)
    fig.append_trace(trace2[1], 1, 2)
    axis_params = dict(showgrid=False, showticklabels=False, zeroline=False, titlefont=dict(color='white'))
    fig['layout'].update(showlegend=False, margin=dict(r=0, l=0, b=0, t=0))
    fig['layout']['scene1'].update(dict(xaxis=axis_params, yaxis=axis_params, zaxis=axis_params,
                                        annotations=annotations1))
    fig['layout']['scene2'].update(dict(xaxis=axis_params, yaxis=axis_params, zaxis=axis_params,
                                        annotations=annotations2))
    offline.plot(fig, show_link=False, filename=plot_name + '.html')