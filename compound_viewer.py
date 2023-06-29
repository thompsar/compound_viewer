import re
import panel as pn
import param
from src.chemistry import ChemLibrary

pn.extension('tabulator', template='fast')


class ChemGUI(param.Parameterized):
    # note, below loads both libraries into memory simultaneously
    # probably not the best as the libraries grow.
    libraries = {'cns': ChemLibrary('cns'), 'diverset': ChemLibrary('diverset')}
    default_string = 'Select library then enter Chembridge Compound ID(s)'
    visible_columns = ['Compound', 'SMILES', 'Molecular Weight', 'LogP']
    # Widgets
    lib_select = param.Selector(objects=libraries, default=libraries['cns'])
    compound_input = param.String()
    selected_compound = param.String()
    load_button = param.Action(lambda x: x.param.trigger('load_button'), label='Load Compounds')
    save_button = param.Action(lambda x: x.param.trigger('save_button'), label='Save Spreadsheet')

    def __init__(self, **params):
        super().__init__(**params)
        library_df = self.libraries['cns'].library[self.visible_columns]
        self.compound_table = pn.widgets.Tabulator(library_df,
                                                   height=800,
                                                   width=800,
                                                   disabled=True,
                                                   visible=False,
                                                   show_index=False)
        self.compound_table.on_click(self.on_click)
        self.compound_image = pn.pane.PNG(object=None, width=300, name='Compound Image')

    @param.depends('lib_select', watch=True)
    def update_compound_df(self):
        self.compound_table.value = self.lib_select.library[self.visible_columns]
        self.compound_input = ''
        self.selected_compound = ''
        self.compound_table.visible = False
        self.compound_image.object = None

    @param.depends('load_button', watch=True)
    def load_compounds(self):
        compounds = re.findall(r'\d+', self.compound_input)
        if compounds:
            new_df = self.lib_select.get_compounds(compounds)[self.visible_columns].reset_index()
            self.compound_table.value = new_df
            self.compound_table.visible = True

    @param.depends('save_button', watch=True)
    def save_spreadsheet(self):
        # NOTE: see link (1) at bottom of this file if you want to implment saving
        # the compound image to the excel file
        self.compound_table.value.to_excel('compound_list.xlsx', index=False)

    def on_click(self, event):
        compound = self.compound_table.value.loc[event.row, 'Compound']
        image = self.lib_select.draw_compound(compound)
        self.selected_compound = '### Compound ID: '+compound
        self.compound_image.object = image


chemGUI = ChemGUI()

pn.Column(
    pn.Param(chemGUI.param,
             parameters=['lib_select', 'compound_input'],
             widgets={'lib_select': pn.widgets.RadioButtonGroup,
                      'compound_input': {'type': pn.widgets.TextAreaInput,
                                         'placeholder': chemGUI.default_string,
                                         'height': 250,
                                         'name': 'Compound ID(s)'}},
             show_name=False),
    pn.Param(chemGUI.param,
             parameters=['load_button', 'save_button'],
             default_layout=pn.Row,
             margin=(-2, 5),
             show_name=False),
    pn.pane.Markdown(chemGUI.param.selected_compound, margin=(0, 10)),
    chemGUI.compound_image
).servable(target='sidebar')

pn.Row(chemGUI.compound_table).servable(target='main', title='ChemBridge Compound Viewer')

# References:
# https://stackoverflow.com/questions/39563065/how-to-insert-an-image-in-an-excel-sheet-using-openpyxl
