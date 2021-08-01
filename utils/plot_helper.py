import plotly.graph_objs as go
from plotly.offline import plot


def simple_plot_2d(dfs, names, plot_path):
    trace_list = []
    ranges = [[-40, 100], [-40, 20], [-40, 20]]
    for i in range(len(dfs)):
        trace = go.Scatter(
            x=dfs[i].iloc[:, 0],
            y=dfs[i].iloc[:, 1],
            text=dfs[i].index,
            mode='markers',
            marker=dict(
                size=3,
            ),
            name=names[i]
        )
        trace_list.append(trace)
    layout = go.Layout(
        scene=dict(
            xaxis=dict(
                title=dict(
                    text='PC1',
                    font=dict(family="sans-serif", size=16, color="black")
                ),
                range=[-40, 100],
            ),
            yaxis=dict(
                title=dict(
                    text='PC2',
                    font=dict(family="sans-serif", size=16, color="black")
                ),
                range=[-40, 20],
            ),
        ),
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=0
        ),
        legend=go.layout.Legend(
            traceorder="normal",
            font=dict(
                family="sans-serif",
                size=20,
                color="black"
            ),
            bordercolor="Grey",
            borderwidth=1
        ),
        template="plotly_white"
    )

    fig = go.Figure(data=trace_list, layout=layout)
    plot(fig, filename=plot_path)


def simple_plot(dfs, names, plot_path):
    trace_list = []
    ranges = [[-40, 100], [-40, 20], [-40, 20]]
    for i in range(len(dfs)):
        trace = go.Scatter3d(
            x=dfs[i].iloc[:, 0],
            y=dfs[i].iloc[:, 1],
            z=dfs[i].iloc[:, 2],
            text=dfs[i].index,
            mode='markers',
            marker=dict(
                size=5,
            ),
            name=names[i]
        )
        trace_list.append(trace)

    layout = go.Layout(
        scene=dict(
            xaxis=dict(
                title=dict(
                    text='PC1',
                    font=dict(family="sans-serif", size=16, color="black")
                ),
                range=ranges[0],
            ),
            yaxis=dict(
                title=dict(
                    text='PC2',
                    font=dict(family="sans-serif", size=16, color="black")
                ),
                range=ranges[1],
            ),
            zaxis=dict(
                title=dict(
                    text='PC3',
                    font=dict(family="sans-serif", size=16, color="black")
                ),
                range=ranges[2],
            ),
        ),
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=0
        ),
        legend=go.layout.Legend(
            traceorder="normal",
            font=dict(
                family="sans-serif",
                size=20,
                color="black"
            ),
            bordercolor="Grey",
            borderwidth=1
        ),
        template="plotly_white"
    )

    fig = go.Figure(data=trace_list, layout=layout)
    plot(fig, filename=plot_path)


def plot_activity_data(df_active, df_na, activity_column, colorbar_title, cid_column, ranges, output_path):
    active_trace = go.Scatter3d(
        x=df_active['PC-1'],
        y=df_active['PC-2'],
        z=df_active['PC-3'],
        mode='markers',
        text=df_active["SMILES"],
        # texttemplate="%{text:.2f}",
        textfont=dict(size=9),
        marker=dict(
            size=6,
            color=df_active[activity_column],
            colorscale='Bluered',
            colorbar_title=colorbar_title,
            opacity=1
        ),
        hovertext=df_active[cid_column],
        customdata=df_active[activity_column],
        # hovertemplate="CID: %{hovertext}<br>Score: %{customdata}<br>(%{x:.2f}, %{y:.2f}, %{z:.2f})",
        hovertemplate="CID: %{hovertext}<br>SMILES: %{text}<br>Score: %{customdata}",
        name="With activity data",
    )
    na_trace = go.Scatter3d(
        x=df_na['PC-1'],
        y=df_na['PC-2'],
        z=df_na['PC-3'],
        mode='markers',
        text=df_na["SMILES"],
        marker=dict(
            size=5,
            color='grey',
            opacity=0.05
        ),
        hovertext=df_na[cid_column],
        # hovertemplate="CID: %{hovertext}<br>(%{x:.2f}, %{y:.2f}, %{z:.2f})",
        hovertemplate="CID: %{hovertext}<br>SMILES: %{text}",
        name="No activity data",
    )
    fig = go.Figure(data=[active_trace, na_trace])
    fig.update_layout(showlegend=False,
                      template='plotly_white',
                      scene=dict(
                          xaxis=dict(
                              title=dict(
                                  text='PC1',
                                  font=dict(size=16)
                              ),
                              range=ranges[0], ),
                          yaxis=dict(
                              title=dict(
                                  text='PC2',
                                  font=dict(size=16)
                              ),
                              range=ranges[1], ),
                          zaxis=dict(
                              title=dict(
                                  text='PC3',
                                  font=dict(size=16)
                              ),
                              range=ranges[2]),
                      ),
                      )
    plot(fig, filename=output_path)


def plot_activity_data_pc2_pc3(df_active, activity_column, colorbar_title, cid_column, ranges, output_path):
    active_trace = go.Scatter(
        x=df_active['PC-2'],
        y=df_active['PC-3'],
        mode='markers+text',
        textposition="bottom center",
        text=df_active["SMILES"],
        # texttemplate="%{text:.2f}",
        textfont=dict(size=12),
        marker=dict(
            size=12,
            color=df_active[activity_column],
            colorscale='Bluered',
            colorbar_title=colorbar_title,
            opacity=1
        ),
        hovertext=df_active[cid_column],
        customdata=df_active[activity_column],
        # hovertemplate="CID: %{hovertext}<br>Score: %{customdata}<br>(%{x:.2f}, %{y:.2f}, %{z:.2f})",
        hovertemplate="CID: %{hovertext}<br>SMILES: %{text}<br>Score: %{customdata}",
        name="With activity data",
    )
    fig = go.Figure(data=[active_trace])
    fig.update_layout(showlegend=False,
                      template='plotly_white',
                      scene=dict(
                          xaxis=dict(
                              title=dict(
                                  text='PC2',
                                  font=dict(size=16)
                              ),
                              range=ranges[0], ),
                          yaxis=dict(
                              title=dict(
                                  text='PC3',
                                  font=dict(size=16)
                              ),
                              range=ranges[1], ),
                      ),
                      )
    plot(fig, filename=output_path)
