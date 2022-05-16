import pygmt
fig = pygmt.Figure()

with fig.subplot(nrows=2, ncols=3, figsize=("15c", "6c"), frame="lrtb"):
    for i in range(2):  # row number starting from 0
        for j in range(3):  # column number starting from 0
            index = i * 3 + j  # index number starting from 0
            with fig.set_panel(panel=index):  # sets the current panel
                fig.text(
                    position="MC",
                    text=f"index: {index}; row: {i}, col: {j}",
                    region=[0, 1, 0, 1],
                )
fig.show()