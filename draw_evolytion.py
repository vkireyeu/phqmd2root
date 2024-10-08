#!/usr/bin/env python3

import sys
try:
    import os
    import argparse
    import ROOT
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    # EXPERIMENTAL
    import pandas as pd
    import plotly.express as px
    import plotly.graph_objects as go

except ModuleNotFoundError as err:
    sys.exit(err)

parser = argparse.ArgumentParser(
                    prog        = 'draw_evolution',
                    description = 'PHSD data processing program')
parser.add_argument('-i', '--input', metavar = 'INPUT', help = 'Input ROOT file', default = 'out.root')
parser.add_argument('-d', '--dir',   metavar = 'DIR',   help = 'Output directory for plots', default = 'plots')



LC_R  = '#E02B3E' # RED
LC_O  = '#E0A336' # ORANGE
LC_P  = '#801FE0' # PURPLE
LC_B  = '#1777E3' # BLUE
LC_G  = '#09E33A' # GREEN
LC_E1 = '#2da4a8'
LC_E2 = '#657a57'


# This is the main function, all other subroutines are called from here
def main():
    global INP_FILE     # Input ROOT file

    # Command line arguments parsing
    args = parser.parse_args()
    input_file      = args.input           # Input file
    outdir          = args.dir

    print(f'Input file: {input_file}')
    check_file_exist(input_file)                  # Check if the input file exists
    ROOT.gSystem.Load('libPHQMDevent')
    os.system(f'mkdir -pv {outdir}')

    INP_FILE = ROOT.TFile.Open(input_file,'READ') # Open and check (if not 'zombie') the input file
    check_root_file(INP_FILE)

    TSACA  = 0
    DTSACA = 0
    inputPHSD  = INP_FILE['inputPHSD']
    for entry in inputPHSD:
        run = entry.run[0]
        TSACA  = run.GetTSACA()
        DTSACA = run.GetDTSACA()
    
    tree  = INP_FILE['Events']
    print(f'Entries: {tree.GetEntries()}')
    for index_ev,entry in enumerate(tree):
        if not index_ev ==  0:
            continue
        event = entry.event[0] # Event in the tree
        mstbsteps = event.GetMstBSteps()
        mstfsteps = event.GetMstFSteps()
        if not len(mstbsteps) == len(mstfsteps):
            break
        common_slices = len(mstbsteps)
        plot_data   = [[]]
        cl_data   = [[]]
        frames = []
        for index_slice in range(common_slices):
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.set_xlim([-75, 75])
            ax.set_ylim([-75, 75])
            ax.set_zlim([-130, 130])
            ax.set_xlabel('x, fm')
            ax.set_ylabel('y, fm')
            ax.set_zlabel('z, fm')
            baryons  = mstbsteps[index_slice]
            clusters = mstfsteps[index_slice]
            print(f'Event {index_ev}, time step {index_slice+1}') 
            plot_data.append([9999, 9999, 9999, 20, LC_G, index_slice+1, 'Cluster'])
            cl_data.append([9999, 9999, 9999, 20, LC_G, index_slice+1, 'Cluster'])
            go_x = [[],[],[],[]]
            go_y = [[],[],[],[]]
            go_z = [[],[],[],[]]
            go_s = [[],[],[],[]]
            go_c = [[],[],[],[]]
            trace = []
            for index_cl,cluster in enumerate(clusters):
                ax.plot(cluster.GetR().X(), cluster.GetR().Y(), cluster.GetR().Z(), marker='o', markersize=2+0.1*cluster.GetA(), color=LC_G)
                plot_data.append([cluster.GetR().X(), cluster.GetR().Y(), cluster.GetR().Z(), 4+0.2*cluster.GetA(), LC_G, index_slice+1, 'Cluster'])
                go_x[0].append(cluster.GetR().X())
                go_y[0].append(cluster.GetR().Y())
                go_z[0].append(cluster.GetR().Z())
                go_c[0].append(LC_G)
                go_s[0].append(7+0.1*cluster.GetA())
            plot_data.append([9999, 9999, 9999, 1, LC_O, index_slice+1, 'Hyperons'])
            for index_bar,baryon in enumerate(baryons):
                if not baryon.IsBound():
                    col = 0
                    pname =''
                    nlist = -1
                    if baryon.GetPdg() == 2212:
                        col=LC_R
                        pname='p'
                        nlist = 1
                    elif baryon.GetPdg() == 2112:
                        col=LC_B
                        pname='n'
                        nlist = 2
                    elif baryon.GetPdg() == 3122 or baryon.GetPdg() == 3212:
                        col=LC_O
                        pname='Hyperons'
                        nlist = 3
                    ax.plot(baryon.GetR().X(), baryon.GetR().Y(), baryon.GetR().Z(), marker='o', markersize=1, color=col)
                    plot_data.append([baryon.GetR().X(), baryon.GetR().Y(), baryon.GetR().Z(), 1, col, index_slice+1, pname])
                    go_x[nlist].append(baryon.GetR().X())
                    go_y[nlist].append(baryon.GetR().Y())
                    go_z[nlist].append(baryon.GetR().Z())
                    go_c[nlist].append(col)
                    go_s[nlist].append(5)
            for part in range(4):
                if part == 0:
                    col=LC_G
                    pname='Cluster'
                    go_x[part].append(999)
                    go_y[part].append(999)
                    go_z[part].append(999)
                    go_c[part].append(LC_G)
                    go_s[part].append(10)
                if part == 1:
                    col=LC_R
                    pname='p'
                elif part == 2:
                    col=LC_B
                    pname='n'
                elif part == 3:
                    col=LC_O
                    pname=r'$\Lambda$'
                    go_x[part].append(999)
                    go_y[part].append(999)
                    go_z[part].append(999)
                    go_c[part].append(LC_O)
                    go_s[part].append(5)
                trace.append(go.Scatter3d(x=go_x[part], y=go_y[part], z=go_z[part], mode='markers', 
                                                marker=dict(size=go_s[part], color=col, opacity=0.8, line=dict(width=0)),
                                                name=pname))
            frames.append(go.Frame(data=[trace[1],trace[2],trace[3],trace[0]],
                                   name=f'time = {TSACA + index_slice*DTSACA} fm/c'))

            plt.title(f'Au+Au, $E_{{kin}}$ = 1.5 AGeV, b = {event.GetB():5.2f} fm, time = {TSACA + index_slice*DTSACA} fm/c')
            entry_p = mpatches.Patch(color=LC_R, label='p')
            entry_n = mpatches.Patch(color=LC_B, label='n')
            entry_l = mpatches.Patch(color=LC_O, label=r'$\Lambda$')
            entry_c = mpatches.Patch(color=LC_G, label='Cluster')
            ax.legend(handles=[entry_p,entry_n,entry_l,entry_c], loc='center left', bbox_to_anchor=(-0.3,0.5))
            plt.savefig(f'{outdir}/im{index_slice:03d}.png', dpi=300)
            plt.close()


        figo = go.Figure(data=[go.Scatter3d(x=[0, 999], y=[0, 999], z=[0, 999],),
                               go.Scatter3d(x=[0, 999], y=[0, 999], z=[0, 999],),
                               go.Scatter3d(x=[0, 999], y=[0, 999], z=[0, 999],),
                               go.Scatter3d(x=[0, 999], y=[0, 999], z=[0, 999],)],
                         frames=frames)
        figo.update_layout(scene = dict(xaxis = dict(title='x, fm', autorange=False, range=[-75,75]),
                                        yaxis = dict(title='y, fm', autorange=False, range=[-75,75]),
                                        zaxis = dict(title='z, fm', autorange=False, range=[-130,130]))
                           )

        def frame_args(duration):
            return {
                    "frame": {"duration": duration},
                    "mode": "immediate",
                    "fromcurrent": True,
                    "transition": {"duration": duration, "easing": "linear"},
                    }

        sliders = [
            {"pad": {"b": 30, "t": 60},
             "len": 0.9,
             "x": 0.1,
             "y": 0,
             "steps": [
                         {"args": [[f.name], frame_args(0)],
                          "label": f'{TSACA + k*DTSACA}',
                          "method": "animate",
                          } for k, f in enumerate(figo.frames)
                      ]
             }
        ]

        figo.update_layout(
            updatemenus = [{"buttons":[
                        {
                            "args": [None, frame_args(10)],
                            "label": "Play", 
                            "method": "animate",
                        },
                        {
                            "args": [[None], frame_args(0)],
                            "label": "Pause", 
                            "method": "animate",
                      }],
                    "direction": "left",
                    "pad": {"r": 10, "t": 70},
                    "type": "buttons",
                    "x": 0.1,
                    "y": 0,
                }
             ],
             sliders=sliders
        )
        figo.layout.updatemenus[0].buttons[0].args[1]['frame']['duration'] = 80
        figo.layout.updatemenus[0].buttons[0].args[1]['transition']['duration'] = 10
        figo.update_layout(scene_aspectmode='cube', autosize=False, width=800, height=800, sliders=sliders,
                          legend=dict(font=dict(size=24),itemwidth=30))

        for k in range(len(figo.frames)):
            figo.frames[k].layout.update(
                autosize=False, width=720, height=720,
                annotations=[
                    dict(
                        x=0, y=1.15,
                        xref='paper', yref='paper',
                        text=f'$\\text{{Au+Au, }} E_{{kin}} = 1.5 \\text{{ AGeV, b = {event.GetB():5.2f} fm, {TSACA + k*DTSACA} fm/c}}$',
                        showarrow=False,
                        font=dict(
                            size=24
                        )
                    )
                ]
            )

        # ~ figo.write_html(f'{outdir}/figo.html', include_plotlyjs = 'cdn', include_mathjax = 'cdn', auto_play = False)
        figo.write_html(f'{outdir}/figo.html', include_plotlyjs = 'cdn', include_mathjax = 'cdn',
                        animation_opts = dict(frame=dict(duration=80), transition=dict(duration=10)))
        # EXPERIMENTAL
        df = pd.DataFrame(plot_data, columns = ['x', 'y', 'z', 'size', 'color', 'step', 'particle'])
        cldf = pd.DataFrame(cl_data, columns = ['x', 'y', 'z', 'size', 'color', 'step', 'particle'])
        fig = px.scatter_3d(df, x='x', y='y', z='z', color='particle', size='size', animation_frame='step',
                            range_x=[-75, 75], range_y=[-75, 75], range_z=[-130,130], size_max=50,
                            color_discrete_map={'Cluster': LC_G, 'p': LC_R, 'n': LC_B, 'Hyperons': LC_O})

        fig.layout.scene.aspectratio = {'x':1, 'y':1, 'z':1}
        fig.update_layout(autosize=False, width=800, height=800)
        fig.write_html(f'{outdir}/fig_express.html')

    os.system(f'ffmpeg -y -framerate 5 -i {outdir}/im%03d.png  -c:v libx264 -crf 0 {outdir}/../mst_evolution.mp4')
    os.system(f'magick -delay 15 -loop 0 {outdir}/im*.png -resize 960 {outdir}/../mst_evolution.gif')



##  This subroutine checks if the file exists
#  \param fname The file name to check
def check_file_exist(fname):
    if(not os.path.isfile(fname)):
        print(f'--- {fname} does not exist')
        exit(-1)



##  This subroutine checks if the ROOT file can be opened and is not 'zombie'
#  \param file The ROOT file to check
def check_root_file(file):
    if (not file or file.IsZombie()):
        print(f'--- Can not open {file}')
        exit(-1)


# This is a truly main function
if __name__ == '__main__':
    main()
    print('All done')
