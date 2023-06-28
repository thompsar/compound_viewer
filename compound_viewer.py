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
    compound_df = param.DataFrame(libraries['cns'].library[visible_columns], precedence=-1)
    load_button = param.Action(lambda x: x.param.trigger('load_button'), label='Load Compounds')
    save_button = param.Action(lambda x: x.param.trigger('save_button'), label='Save Spreadsheet')

    @param.depends('lib_select', watch=True)
    def update_compound_df(self):
        self.compound_df = self.lib_select.library[self.visible_columns]
        self.compound_input = ''
        self.param.compound_df.precedence = -1

    @param.depends('load_button', watch=True)
    def load_compounds(self):
        compounds = re.findall(r'\d+', self.compound_input)
        if compounds:
            self.compound_df = self.lib_select.get_compounds(compounds)[self.visible_columns]
            self.param.compound_df.precedence = 0

    @param.depends('save_button', watch=True)
    def save_spreadsheet(self):
        self.compound_df.to_excel('compound_list.xlsx', index=False)


chemGUI = ChemGUI()

pn.Param(chemGUI.param,
         parameters=['lib_select', 'compound_input'],
         widgets={'lib_select': pn.widgets.RadioButtonGroup,
                  'compound_input': {'type': pn.widgets.TextAreaInput,
                                     'placeholder': chemGUI.default_string,
                                     'height': 400,
                                     'name': 'Compound ID(s)'}},
         show_name=False).servable(target='sidebar')

pn.Param(chemGUI.param,
         parameters=['load_button', 'save_button'],
         default_layout=pn.Row,
         margin=(-2, 5),
         show_name=False).servable(target='sidebar')

pn.Param(chemGUI.param,
         parameters=['compound_df'],
         widgets={'compound_df': {'type': pn.widgets.Tabulator,
                                  'disabled': True}},
         show_name=False).servable(target='main', title='ChemBridge Compound Viewer')
