import re
import pandas as pd
import panel as pn
from src.chemistry import draw_compound, get_smiles

pn.extension('tabulator')

placeholder = 'Enter Chembridge Compound ID(s) separated by commas or new lines'
compound_input = pn.widgets.input.TextAreaInput(name='Compound ID(s)',
                                                placeholder=placeholder,
                                                height=400)
compound_image = pn.pane.PNG(object=None, height=400, width=400)
compound_slider = pn.widgets.IntSlider(name='Select Compound',
                                       start=0,
                                       end=1,
                                       disabled=True)
# dataframe widget
df_dict = {'Compound ID': [], 'SMILES': []}
df = pd.DataFrame(df_dict, dtype='object')
compound_table = pn.widgets.Tabulator(df, height=400, width=800)


def update_compounds(event):
    compounds = re.findall(r'\d+', compound_input.value)
    if compounds:
        compound_slider.disabled = False
        compound_slider.start = 0
        compound_slider.end = len(compounds) - 1
        compound_slider.value = 0
        compound = compounds[0]
        image = draw_compound(compound)
        compound_image.object = image

        df_dict = {'Compound ID': compounds,
                   'SMILES': [get_smiles(compound) for compound in compounds]}
        df = pd.DataFrame(df_dict)
        compound_table.value = df
    else:
        compound_slider.disabled = True
        compound_slider.value = 0
        compound_image.object = None


def update_slider(event):
    compounds = re.findall(r'\d+', compound_input.value)
    if compounds:
        compound = compounds[event.new]
        image = draw_compound(compound)
        compound_image.object = image
        # compound_table.param.set_property(selected=[event.new])
        compound_table.selection = [event.new]
    else:
        compound_image.object = None


compound_input.param.watch(update_compounds, 'value')
compound_slider.param.watch(update_slider, 'value')


app = pn.Column(
    pn.Row(pn.Column(compound_input, compound_slider, styles=dict(background='WhiteSmoke')),
           pn.Column(compound_image)),
    pn.Row(compound_table)
    )

app.servable()
