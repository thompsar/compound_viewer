import warnings
import requests
import pandas as pd
import numpy as np
from bs4 import BeautifulSoup
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors


class ChemLibrary:
    """
    Class for loading info about ChemBridge compound library
    """

    def __init__(self, library_name):
        self.library_name = library_name.lower()
        self.load_library()

    def load_library(self):
        """
        Loads chemical library file
        """
        self.library_file = self.library_name + '.csv'
        try:
            self.library = pd.read_csv(self.library_file, dtype='object')
        except FileNotFoundError:
            raise FileNotFoundError('Library not found, supported libraries are cns and diverset')

    def get_compounds(self, compound_ids):
        """Gets compound(s) info from ChemBridge servers if not locally avaialble in library file

        Parameters
        ----------
        compound_id : list of str, int or single str, int
            Chembridge compound ID(s)

        Returns
        -------
        if list:
        pd.DataFrame: chem info including compound name, SMILES string, other details and molfile
        if str or int:
        pd.Series: chem info including compound name, SMILES string, other details and molfile
        """
        # convert compound_id to string
        if not isinstance(compound_ids, list):
            compound_ids = [compound_ids]
        compound_ids = [str(compound_id) for compound_id in compound_ids]
        update_library = False
        for compound_id in compound_ids:
            chem_info = self.library[self.library['Compound'] == compound_id]

            if chem_info.empty:
                warnings.warn('Compound ' + compound_id + ' not found in library')
                continue
            if chem_info.isna().values.any():
                update_library = True
                print('Loading '+compound_id+' structure from ChemBridge')
                # get basic compound info
                self.fetch_info(compound_id, chem_info.index)
                # get compound structure
                self.fetch_structure(compound_id, chem_info.index)
                # generate smiles string and update MW
                self.generate_smiles(chem_info.index)
                # update chem_info
                chem_info = self.library[self.library['Compound'] == compound_id]
        # save updated library
        if update_library:
            self.library.to_csv(self.library_file, index=False)
        # if len(compound_ids) == 1:
        #     return chem_info.squeeze()
        return self.library[self.library['Compound'].isin(compound_ids)]

    def fetch_info(self, compound_id, index):
        """
        Gets compound info from ChemBridge servers
        Use only if compound info is not available locally or if you want to update the info
        """
        info_url = 'https://www.hit2lead.com/screening-compounds/'
        # get basic compound info
        page = requests.get(info_url + compound_id, timeout=5)
        soup = BeautifulSoup(page.content, 'html.parser')
        table = soup.find('div', {'class': 'matter'})
        self.library.loc[index, 'Compound Name'] = table.find('p').text.strip()
        dt_tags = soup.find_all('dt')
        for dt_tag in dt_tags:
            category = dt_tag.text.strip()
            dd_tag = dt_tag.find_next_sibling('dd')
            value = dd_tag.text.strip()
            self.library.loc[index, category] = value

    def fetch_structure(self, compound_id, index):
        """
        Gets compound structure from ChemBridge servers
        Use only if compound structure is not available locally or if you want to update the info
        """
        structure_url = 'https://www.hit2lead.com/search.asp?db=SC&SearchPage=structure&id='
        # get compound structure
        page = requests.get(structure_url + compound_id, timeout=5)
        soup = BeautifulSoup(page.content, 'html.parser')
        molfile = soup.find('input', {'id': 'molfile'}).get('value')
        self.library.loc[index, 'Structure'] = molfile

    def generate_smiles(self, index):
        """
        Generates a SMILES string for a compound
        Also updates MW as Chembridge often rounds MW to nearest integer

        For printing SMILES run get_compound(compound_ID).SMILES
        """

        molfile = self.library.loc[index, 'Structure'].values[0]
        mol = Chem.MolFromMolBlock(molfile)
        self.library.loc[index, 'SMILES'] = Chem.MolToSmiles(mol)
        self.library.loc[index, 'Molecular Weight'] = np.round(Descriptors.MolWt(mol),2)

    def draw_compound(self, compound_id, transparent=False, legend=None):
        """Draws Chembridge compound

        Parameters
        ----------
        compound_id : str, int
            Chembridge compound ID
        transparent : bool, optional
            draw grid with transparent background, by default False
        legend : str, optional
            include legend (such as compound_id in str), font size broken

        Returns
        -------
        image : PIL image of compounds in grid
        """

        dopts = Draw.rdMolDraw2D.MolDrawOptions()
        dopts.legendFontSize = 80

        chem_info = self.get_compounds(compound_id)
        molfile = chem_info['Structure']

        mol = Chem.MolFromMolBlock(molfile)
        img = Draw.MolToImage(mol, size=(600, 600), legend=legend, drawOptions=dopts)
        if transparent:
            img = make_transparent(img)
        return img

    def compound_grid(self, compound_ids, mols_per_row=6, transparent=False):
        """Draws grid of Chembridge compounds

        Parameters
        ----------
        compound_ids : list of str, int
            Chembridge compound IDs
        molsPerRow : int, optional
            number of molecules drawn per row, by default 6
        transparent : bool, optional
            draw grid with transparent background, by default False

        Returns
        -------
        image : PIL image of compounds in grid
        """
        # options for drawing molecules located here:
        # http://rdkit.org/docs/cppapi/structRDKit_1_1MolDrawOptions.html
        dopts = Draw.rdMolDraw2D.MolDrawOptions()
        dopts.legendFontSize = 80
        dopts.drawMolsSameScale = False

        mols = {}
        for compound_id in compound_ids:
            chem_info = self.get_compounds(compound_id)
            molfile = chem_info['Structure']
            # force str for compound_id in dict key to eliminate dupes
            # and prevent error in Draw.MolsToGridImage legends
            mols[str(compound_id)] = Chem.MolFromMolBlock(molfile)

        img = Draw.MolsToGridImage(list(mols.values()),
                                   molsPerRow=mols_per_row,
                                   subImgSize=(600, 600),
                                   legends=list(mols.keys()),
                                   returnPNG=False,
                                   drawOptions=dopts)

        if transparent:
            img = make_transparent(img)
        return img


def make_transparent(img):
    """
    Makes white background transparent
    """
    img = img.convert("RGBA")
    datas = img.getdata()

    # make all white pixels transparent
    # https://stackoverflow.com/questions/765736/how-to-use-pil-to-make-all-white-pixels-transparent
    newdata = []
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            newdata.append((255, 255, 255, 0))
        else:
            newdata.append(item)

    img.putdata(newdata)

    return img
