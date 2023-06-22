import re
import panel as pn
from src.chemistry import ChemLibrary

pn.extension('tabulator')

class ChemGUI():
    
    def __init__(self):
        # Variables
        self.placeholder = 'Select library then enter Chembridge Compound ID(s) separated by commas or new lines'
        self.libraries = ['cns', 'diverset']
        self.visible_columns = ['Compound', 'SMILES', 'Molecular Weight', 'LogP']
        # Widgets
        self.library_select = pn.widgets.RadioButtonGroup(name='Library',
                                                    options=self.libraries,
                                                    value=self.libraries[0])
        self.compound_library = ChemLibrary(self.library_select.value)

        self.compound_input = pn.widgets.input.TextAreaInput(name='Compound ID(s)',
                                                        placeholder=self.placeholder,
                                                        height=400)

        self.compound_load = pn.widgets.Button(name='Load Compounds', button_type='primary')
        self.library_subset = self.compound_library.get_compounds([])
        self.compound_table = pn.widgets.Tabulator(self.library_subset[self.visible_columns], height=800, width=800)

        self.library_select.param.watch(self.update_library, 'value')
        self.compound_load.on_click(self.update_compounds)

    def update_library(self, event):
            
            self.compound_library = ChemLibrary(self.library_select.value)
            self.compound_input.value = ''
            self.library_subset = self.compound_library.get_compounds([])
            self.compound_table.value = self.library_subset[self.visible_columns]

    def update_compounds(self, event):
        compounds = re.findall(r'\d+', self.compound_input.value)
        if compounds:
            self.library_subset = self.compound_library.get_compounds(compounds)
            self.compound_table.value = self.library_subset[self.visible_columns]
        

    

chemGUI = ChemGUI()

pn.Row(
     pn.Column(chemGUI.library_select,
               chemGUI.compound_input,
               chemGUI.compound_load,
               styles=dict(background='WhiteSmoke')),
     pn.Column(chemGUI.compound_table)
     ).servable()
