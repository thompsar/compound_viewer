import requests
import pandas as pd
from bs4 import BeautifulSoup
from rdkit import Chem
from rdkit.Chem import Draw


def get_structure(compound_id):
    """Gets strutucture from ChemBridge servers if not locally avaialble in chem_structures folder

    Parameters
    ----------
    compound_id : str, int
        Chembridge compound ID

    Returns
    -------
    molfile : Accelyrs molecule coordinates file
    """
    # convert compound_id to string
    compound_id = str(compound_id)
    info_url = 'https://www.hit2lead.com/screening-compounds/'
    structure_url = 'https://www.hit2lead.com/search.asp?db=SC&SearchPage=structure&id='

    try:
        chem_info = pd.read_csv('chem_structures/'+compound_id+'_structure.csv',
                                header=0, index_col=0).squeeze('columns')
    except FileNotFoundError:
        print('Loading structure from ChemBridge')
        chem_info = pd.Series()
        chem_info['Compound ID'] = compound_id
        # get basic compound info
        page = requests.get(info_url + compound_id, timeout=5)
        soup = BeautifulSoup(page.content, 'html.parser')
        table = soup.find('div', {'class': 'matter'})
        chem_info['Compound Name'] = table.find('p').text.strip()
        dt_tags = soup.find_all('dt')
        for dt_tag in dt_tags:
            category = dt_tag.text.strip()
            dd_tag = dt_tag.find_next_sibling('dd')
            value = dd_tag.text.strip()
            chem_info[category] = value

        # get compound structure
        page = requests.get(structure_url + compound_id, timeout=5)
        soup = BeautifulSoup(page.content, 'html.parser')
        molfile = soup.find('input', {'id': 'molfile'}).get('value')

        chem_info['Structure'] = molfile
        chem_info.to_csv('chem_structures/'+compound_id+'_structure.csv')
    return chem_info


def get_smiles(compound_id):
    """
    Gets SMILES string for Chembridge compound
    """
    chem_info = get_structure(compound_id)
    molfile = chem_info['Structure']
    mol = Chem.MolFromMolBlock(molfile)
    return Chem.MolToSmiles(mol)


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


def draw_compound(compound_id, transparent=False, legend=None):
    """
    Draws Chembridge compound
    """
    dopts = Draw.rdMolDraw2D.MolDrawOptions()
    dopts.legendFontSize = 80

    chem_info = get_structure(compound_id)
    molfile = chem_info['Structure']

    mol = Chem.MolFromMolBlock(molfile)
    img = Draw.MolToImage(mol, size=(600, 600), legend=legend, drawOptions=dopts)
    if transparent:
        img = make_transparent(img)
    return img


def compound_grid(compound_ids, mols_per_row=6, transparent=False):
    """Draws grid of Chembridge compounds

    Parameters
    ----------
    compound_ids : list of str, int
        Chembridge compound IDs
    molsPerRow : int, optional
        number of molecules drawn per row, by default 6
    transparent : bool, optional
        draw grid with transparent backgroun, by default False

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
        chem_info = get_structure(compound_id)
        molfile = chem_info['Structure']
        # force str for compound_id in dict key to eliminate dupes
        # and prevent error in Draw.MolsToGridImage legends
        mols[str(compound_id)] = Chem.MolFromMolBlock(molfile)

    img = Draw.MolsToGridImage(list(mols.values()), molsPerRow=mols_per_row, subImgSize=(600, 600),
                               legends=list(mols.keys()), returnPNG=False, drawOptions=dopts)

    if transparent:
        img = make_transparent(img)
    return img
