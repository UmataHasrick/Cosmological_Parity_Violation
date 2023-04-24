import matplotlib.pyplot as plt

plot_vertices_colorscheme = ['black','red','yellow','blue']


def Vertices_plot(ax, color_index, vertices_position):
    """To plot a single 

    Args:
        ax (_type_): _description_
        color_index (_type_): _description_
        vertices_position (_type_): _description_

    Returns:
        _type_: _description_
    """    
    color = plot_vertices_colorscheme[color_index]

    xs = vertices_position[:,0]
    ys = vertices_position[:,1]
    zs = vertices_position[:,2]
    
    ax.scatter(xs, ys, zs, s=20, c=color, depthshade=True)
 
    return ax

def Vertices_single_plot(ax, color_index, vertice_position):
    color = plot_vertices_colorscheme[color_index]

    xs = vertice_position[0]
    ys = vertice_position[1]
    zs = vertice_position[2]
    
    ax.scatter(xs, ys, zs, s=20, c=color, depthshade=True)
 
    return ax

def Line_single_plot(ax, color_index, line_start, line_end):
    xs = [line_start[0], line_end[0]]
    ys = [line_start[1], line_end[1]]
    zs = [line_start[2], line_end[2]]
    ax.plot(xs, ys, zs, c = plot_vertices_colorscheme[color_index], ls = '-')
    return ax

def Line_multiple_plot(ax, multiple_tetrahedron):
    length = multiple_tetrahedron.shape[0]

    for i in range(length):
        for j in range(1,4):
            ax = Line_single_plot(ax, j, multiple_tetrahedron[i][0], multiple_tetrahedron[i][j])

    return ax

def Tetrahedron_plot(ax, multiple_tetrahedron):
    for i in range(4):
        ax = Vertices_plot(ax, i, multiple_tetrahedron[::,i])
    
    for j in range(1,4):
        ax = Line_multiple_plot(ax, multiple_tetrahedron)

    return ax