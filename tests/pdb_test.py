import pytest
from unittest.mock import patch
import pypdb
from aiondata import PDBHandler


@pytest.fixture
def pdb_handler():
    return PDBHandler()


def test_get_pdb(pdb_handler):
    """Test that the PDB files are retrieved."""
    pdb_ids = ["8IRB", "100D"]
    with patch.object(
        pdb_handler.pdb_list, "retrieve_pdb_file"
    ) as mock_retrieve_pdb_file:
        pdb_handler.get_pdb(pdb_ids)
        mock_retrieve_pdb_file.assert_any_call(
            pdb_ids[0], pdir=pdb_handler.save_dir, file_format="pdb"
        )
        mock_retrieve_pdb_file.assert_any_call(
            pdb_ids[1], pdir=pdb_handler.save_dir, file_format="pdb"
        )


def test_get_pdb_info(pdb_handler):
    """Test that the information about a PDB file is retrieved."""
    pdb_id = "4OS0"
    with patch.object(pypdb, "get_all_info") as mock_get_all_info:
        pdb_handler.get_pdb_info(pdb_id)
        mock_get_all_info.assert_called_with(pdb_id)


def test_search_pdb(pdb_handler):
    """Test that the PDB files are searched for."""
    title = "Solution"
    organism = "9606"
    Uniprot_accession = "P04637"
    fromdb = "UniProt"
    nonpolymer = 1
    experiment = "SOLUTION NMR"
    ComparisonType = "Less"
    with patch(
        "aiondata.raw.protein_structure.perform_search_with_graph"
    ) as mock_perform_search_with_graph:
        pdb_handler.search_pdb(
            title=title,
            organism=organism,
            Uniprot_accession=Uniprot_accession,
            fromdb=fromdb,
            experiment=experiment,
            nonpolymer=nonpolymer,
            ComparisonType=ComparisonType,
        )
        mock_perform_search_with_graph.assert_called()


def test_fetch_PDB_uniprot_accession():
    """Test that UniProt accessions are fetched correctly using a PDB ID"""
    pdbid = "3sn6"

    # Sample mock response to simulate the data returned from the API
    mock_response = {
        "data": {
            "entries": [
                {
                    "polymer_entities": [
                        {
                            "rcsb_id": "3SN6_1",
                            "rcsb_polymer_entity_container_identifiers": {
                                "reference_sequence_identifiers": [
                                    {
                                        "database_accession": "P04896",
                                        "database_name": "UniProt",
                                    }
                                ]
                            },
                        },
                        {
                            "rcsb_id": "3SN6_2",
                            "rcsb_polymer_entity_container_identifiers": {
                                "reference_sequence_identifiers": [
                                    {
                                        "database_accession": "P54311",
                                        "database_name": "UniProt",
                                    }
                                ]
                            },
                        },
                    ]
                }
            ]
        }
    }

    with patch("requests.post") as mock_post:
        mock_post.return_value.json.return_value = mock_response
        mock_post.return_value.raise_for_status = (
            lambda: None
        )  # To bypass raise_for_status

        # Call the function under test
        result = pdb_handler.fetch_PDB_uniprot_accession(pdbid)

        # Assert that the requests.post was called with the correct URL
        mock_post.assert_called_with(
            "https://data.rcsb.org/graphql", json={"query": ANY}
        )

        # Verify the result is as expected
        expected_result = {"3SN6_1": ["P04896"], "3SN6_2": ["P54311"]}
        assert result == expected_result


def test_fetch_uniprot_sequence(pdb_handler):
    """Test that a Sequence are fetched using a Uniprot accession"""
    uniprot = "P00504"
    with patch(
        "aiondata.raw.protein_structure.PDBHandler.fetch_uniprot_sequence"
    ) as mock_uniprot:
        pdb_handler.fetch_uniprot_sequence(uniprot)
        mock_uniprot.assert_called()


@pytest.mark.integration
def test_fetch_uniprot_sequence_with_internet_integration(pdb_handler):
    """Test that a Sequence are fetched using a Uniprot accession, by actually calling the UniProt API"""
    uniprot = "P00504"
    result = pdb_handler.fetch_uniprot_sequence(uniprot)
    assert "ADFREDGDSRKVNLGVGAYRTDEGQPWVLPVVRKVEQLIAGDGSLNHEYLPILGLPEFRA" in result


@pytest.mark.integration
def test_fetch_PDB_uniprot_accession_with_internet_integration(pdb_handler):
    """Test that Uniprot accession are fetched using a PDB ID, by actually calling the PDB API"""
    pdbid = "IAAT"
    result = pdb_handler.fetch_PDB_uniprot_accession(pdbid)
    # TODO: Write a correct assert here
