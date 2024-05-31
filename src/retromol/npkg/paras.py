# -*- coding: utf-8 -*-

"""Functions for interaction witn PARAS web service."""

import json
import logging
import typing as ty
from enum import Enum

import requests


class ParasInputType(Enum):
    """Input type for PARAS web service.
    
    :cvar FASTA: Fasta format.
    :cvar GENBANK: Genbank format.
    """
    FASTA = "fasta"
    GENBANK = "genbank"


class ParasSubstrateChoice(Enum):
    """Substrate choice for PARAS web service.
    
    :cvar ALL_SUBSTRATES: All substrates.
    :cvar COMMON_SUBSTRATES: Common substrates.
    """
    ALL_SUBSTRATES = "allSubstrates"
    COMMON_SUBSTRATES = "commonSubstrates"


def predict_specificity(
    src: str,
    selected_input_type: ParasInputType = ParasInputType.FASTA,
    save_active_site_signatures: bool = False,
    save_extended_signatures: bool = False,
    save_adenylation_domain_sequences: bool = False,
    selected_substrate_choice: ParasSubstrateChoice = ParasSubstrateChoice.ALL_SUBSTRATES,
    num_predictions_to_report: int = 1,
    use_structure_guided_profile_alignment: bool = False,
    first_separator: str = "|",
    second_separator: str = "_",
    third_separator: str = "-",
) -> ty.List[ty.Tuple[str, ty.List[ty.Tuple[str, float]]]]:
    """Predict specificity of a given sequence.

    :param src: Sequence to predict specificity for.
    :type src: str
    :param selected_input_type: Input type for the sequence, defaults to ParasInputType.FASTA
    :type selected_input_type: ParasInputType, optional
    :param save_active_site_signatures: Save active site signatures, defaults to False
    :type save_active_site_signatures: bool, optional
    :param save_extended_signatures: Save extended signatures, defaults to False
    :type save_extended_signatures: bool, optional
    :param save_adenylation_domain_sequences: Save adenylation domain sequences, defaults to False
    :type save_adenylation_domain_sequences: bool, optional
    :param selected_substrate_choice: Substrate choice, defaults to ParasSubstrateChoice.ALL_SUBSTRATES
    :type selected_substrate_choice: ParasSubstrateChoice, optional
    :param num_predictions_to_report: Number of predictions to report, defaults to 1
    :type num_predictions_to_report: int, optional
    :param use_structure_guided_profile_alignment: Use structure guided profile alignment, defaults to False
    :type use_structure_guided_profile_alignment: bool, optional
    :param first_separator: First separator, defaults to "|"
    :type first_separator: str, optional
    :param second_separator: Second separator, defaults to "_"
    :type second_separator: str, optional
    :param third_separator: Third separator, defaults to "-"
    :type third_separator: str, optional
    :return: predictions
    :rtype: ty.Dict[str, ty.List[ty.Tuple[str, float]]]
    """
    logger = logging.getLogger(__name__)

    url = "https://paras.bioinformatics.nl/api/submit_paras"
    data = {
        "data": {
            "src": src,
            "selectedInputType": selected_input_type.value,
            "saveActiveSiteSignatures": save_active_site_signatures,
            "saveExtendedSignatures": save_extended_signatures,
            "saveAdenylationDomainSequences": save_adenylation_domain_sequences,
            "selectedSubstrateChoice": selected_substrate_choice.value,
            "numPredictionsToReport": num_predictions_to_report,
            "useStructureGuidedProfileAlignment": use_structure_guided_profile_alignment,
            "firstSeparator": first_separator,
            "secondSeparator": second_separator,
            "thirdSeparator": third_separator,
        }
    }
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, data=json.dumps(data), headers=headers).json()
    predictions = []

    if payload := response.get("payload", None):
        if results := payload.get("results", None):
            
            for domain in results:
                domain_id = domain["domain_id"]
                domain_predictions = domain["data"]["predictions"]
                predictions.append((domain_id, domain_predictions))

    return predictions
