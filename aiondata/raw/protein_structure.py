from ..datasets import ExcelDataset, CsvDataset, CachedDataset
from Bio import PDB
import pypdb
from pypdb.clients.search.operators import text_operators
from pypdb.clients.search.search_client import (
    QueryGroup,
    LogicalOperator,
    ReturnType,
    perform_search_with_graph,
)
import requests


class FoldswitchProteinsTableS1A(ExcelDataset):
    """(A) List of pairs (PDBIDs), lengths and the sequence of the fold-switching region.
    (For those pairs not having the second fold solved in PDB, only the first PDB is reported).

    From Paper: AlphaFold2 fails to predict protein fold switching
    https://doi.org/10.1002/pro.4353
    """

    SOURCE = "https://raw.githubusercontent.com/tomshani/aiondata/tom-branch/data/pro4353-sup-0002-tables1%20/Table_S1A_final.xlsx"
    COLLECTION = "foldswitch_proteins"


class FoldswitchProteinsTableS1B(ExcelDataset):
    """(B) RMSD, TM-scores for the whole protein and only fold-switching fragment,
    as well as sequence identities between the fold-switching pairs.
    wTM-score/wRMSD indicate TM-scores/RMSDs considering whole protein chains.
    fsTM-score/fsRMSD indicate TM-scores/RMSDs considering fold-switching regions only.

    From Paper: AlphaFold2 fails to predict protein fold switching
    https://doi.org/10.1002/pro.4353
    """

    SOURCE = "https://raw.githubusercontent.com/tomshani/aiondata/tom-branch/data/pro4353-sup-0002-tables1%20/Table_S1B_final.xlsx"
    COLLECTION = "foldswitch_proteins"


class FoldswitchProteinsTableS1C(ExcelDataset):
    """(C) List of fold-switching protein pairs (PDBID and chain) used for the analysis,
    first column corresponds to Fold1 and second to Fold2, followed by TM-scores of the predictions.
    Tables attached separately.

    From Paper: AlphaFold2 fails to predict protein fold switching
    https://doi.org/10.1002/pro.4353
    """

    SOURCE = "https://raw.githubusercontent.com/tomshani/aiondata/tom-branch/data/pro4353-sup-0002-tables1%20/Table_S1C_final.xlsx"
    COLLECTION = "foldswitch_proteins"


class CodNas91(CsvDataset):
    """
    Paper: Impact of protein conformational diversity on AlphaFold predictions
    https://doi.org/10.1093/bioinformatics/btac202
    We selected 91 proteins (Supplementary Table S1) with different degrees of conformational diversity expressed as the range of pairwise global Cα-RMSD between their conformers in the PDB (Fig. 1).
    All the pairs of conformers for each protein are apo–holo pairs selected from the CoDNaS database (Monzon et al., 2016) and bibliography. Manual curation for each protein confirmed that structural deformations were associated with a given biological process based on experimental evidence.
    This step is essential to ensure that conformational diversity is not associated with artifacts, misalignments, missing regions, or the presence of flexible ends. When more than two conformers were known, we selected the apo–holo pair showing the maximum Cα-RMSD (maxRMSD).
    Other considerations were absence of disorder, PDB resolution, absence of mutations and sequence differences. We previously observed that when conformational diversity is derived from experimentally based conformers, different ranges of RMSD are obtained between them depending on the structure determination method (Monzon et al., 2017a).
    Here we considered a continuum of protein flexibility measured as the RMSD between apo and holo forms as shown in Figure 1.
    """

    SOURCE = "https://raw.githubusercontent.com/tomshani/aiondata/tom-branch/data/Supplementary_Table_1_91_apo_holo_pairs.csv"


class PDBHandler(CachedDataset):
    """
    A class for handling PDB files.

    Attributes:
    - COLLECTION: The collection name for the PDB files.

    Methods:
    - get_pdb: Retrieves PDB files from the PDB database.
    - get_pdb_info: Retrieves information about a specific PDB file.
    - get_ligand_info: Retrieves information about ligands in a specific PDB file.
    - searchpdb: Performs a search in the PDB database based on specified criteria.
    """

    COLLECTION = "PDB_files"

    def __init__(self):
        """
        Initializes the PDBHandler object.

        Parameters:
        - None

        Returns:
        - None
        """
        self.pdb_list = PDB.PDBList()
        self.save_dir = self.get_cache_path()

    def get_pdb(self, pdb_ids, file_format="pdb"):
        """
        Retrieves PDB files from the PDB database.

        Parameters:
        - pdb_ids: A string or a list of PDB IDs.
        - file_format: The format of the retrieved PDB files (default: 'pdb').

        Returns:
        - None
        """
        if isinstance(pdb_ids, str):
            self.pdb_list.retrieve_pdb_file(
                pdb_id, pdir=self.save_dir, file_format=file_format
            )
        else:
            for pdb_id in pdb_ids:
                self.pdb_list.retrieve_pdb_file(
                    pdb_id, pdir=self.save_dir, file_format=file_format
                )

    def get_pdb_info(self, pdb_id):
        """
        Retrieves information about a specific PDB file.

        Parameters:
        - pdb_id: The ID of the PDB file.

        Returns:
        - The information about the PDB file.
        """
        res = pypdb.get_all_info(pdb_id)
        return res

    def search_pdb(
        self,
        title=None,
        fromdb=None,
        organism=None,
        Uniprot_accession=None,
        nonpolymer=None,
        ComparisonType=None,
    ):
        """
        Perform a search in the Protein Data Bank (PDB) based on the specified criteria.

        Args:
            title (str, optional): Title name to search for. Defaults to None.
            fromdb (str, optional): Database name to search in. Defaults to None.
            organism (str, optional): Organism taxonomy ID to search for. Defaults to None.
            Uniprot_accession (str, optional): Uniprot accession number to search for. Defaults to None.
            nonpolymer (int, optional): Number of non-polymer entities to compare. Defaults to None.
            ComparisonType (str, optional): Comparison type for nonpolymer comparison. Must be 'Greater' or 'Less'. Defaults to None.

        Returns:
            list: List of search results from the Protein Data Bank.

        Raises:
            ValueError: If ComparisonType is not 'Greater' or 'Less' when nonpolymer is provided.

        Examples:
        search_pdb(Uniprot_accession="P04637",title="Solution",organism="9606",fromdb="UniProt")
        search_pdb(nonpolymer=1,ComparisonType="Less")

        """

        # title name search
        if title:
            title = text_operators.ContainsPhraseOperator(
                value=title, attribute="struct.title"
            )

        # Uniprot accession number search
        if Uniprot_accession:
            Uniprot_accession = text_operators.InOperator(
                values=[Uniprot_accession],
                attribute="rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
            )

        # organism search
        # option: "9606"
        if organism:
            organism = text_operators.InOperator(
                values=[organism],
                attribute="rcsb_entity_source_organism.taxonomy_lineage.id",
            )

        # is in certain db
        # "uniprot"
        if fromdb:
            fromdb = text_operators.ExactMatchOperator(
                value=fromdb,
                attribute="rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
            )

        if nonpolymer is not None:
            if ComparisonType not in ["Greater", "Less"]:
                raise ValueError("ComparisonType must be 'Greater' or 'Less'")
            if ComparisonType == "Greater":
                ComparisonType = text_operators.ComparisonType.GREATER
            else:
                ComparisonType = text_operators.ComparisonType.LESS
            nonpolymer = text_operators.ComparisonOperator(
                value=nonpolymer,
                attribute="rcsb_entry_info.nonpolymer_entity_count",
                comparison_type=ComparisonType,
            )

        queries = [fromdb, title, organism, Uniprot_accession, nonpolymer]
        queries = [query for query in queries if query is not None]

        if len(queries) > 1:
            search_operator = QueryGroup(
                queries=queries, logical_operator=LogicalOperator.AND
            )
        else:
            search_operator = queries[0]

        results = perform_search_with_graph(
            query_object=search_operator, return_type=ReturnType.ENTRY, verbosity=True
        )

        return results
    
    

    def fetch_PDB_uniprot_accession(pdb_id):
        """
        Fetches the UniProt accession number for a given PDB ID.
        Problem, this PDB ID does not use chain information. If you need to use chain information, do not remove the chain and run, this may give you the wrong UniProt accession.
        
        Parameters:
            pdb_id (str): The PDB ID of the protein.

        Returns:
            str: The UniProt accession number if available, otherwise None.
        """
        url = "https://data.rcsb.org/graphql"
        
        query = """
        {
        entries(entry_ids: ["%s"]) {
            polymer_entities {
            rcsb_id
            rcsb_polymer_entity_container_identifiers {
                reference_sequence_identifiers {
                database_accession
                database_name
                }
            }
            }
        }
        }
        """ % pdb_id
        
        try:
            response = requests.post(url, json={'query': query})
            response.raise_for_status()  # Raises an HTTPError for bad responses (4xx or 5xx)
            data = response.json()
            
            # Extracting the UniProt accession number
            try:
                db_accession = data['data']['entries'][0]['polymer_entities'][0]['rcsb_polymer_entity_container_identifiers']['reference_sequence_identifiers'][0]['database_accession']
                db_name = data['data']['entries'][0]['polymer_entities'][0]['rcsb_polymer_entity_container_identifiers']['reference_sequence_identifiers'][0]['database_name']
                return db_accession, db_name
            except (KeyError, IndexError) as e:
                print(f"Error extracting data: {e}")
                return None
        
        except requests.RequestException as e:
            print(f"Request failed: {e}")
            return None



    def fetch_uniprot_sequence(uniprot_id):
        """
        Fetches the sequence for a given UniProt accession number.

        Parameters:
            uniprot_id (str): The UniProt accession number.
        
        Returns:
            str: The protein sequence if available, otherwise None.
        
        """
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        
        response = requests.get(url)
        
        if response.status_code == 200:
            # Split the response text into lines
            lines = response.text.split('\n')
            
            # The first line is the header, we'll skip it
            sequence = ''.join(lines[1:])
            
            return sequence.strip()
        else:
            print(f"Failed to fetch sequence. Status code: {response.status_code}")
            print(response.text)
            return None

