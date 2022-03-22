#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full Apache License, Version 2.0 can be found in `LICENSE` at the project root.

import sys
from typing import List
import logging
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from nanodesign.data.dna_structure import DnaStructure as Structure
from nanodesign.data.strand import DnaStrand as Strand

def init_logging():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler()
    formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    return logger
    

logger = init_logging()


def _max_melting_temperatures(dna_structure: Structure, c_mg:float, c_na:float) -> List[float]:
    """ max_melting_temperature: highest domain metling temperature of a staple
        using santaLucia nearest neighbor model 
    """
    melting_temperatures = []
    for domains in (staple.domain_list for staple in dna_structure.strands if not staple.is_scaffold):
        domain_temperatures = []
        for domain in domains:
            sequence = domain.sequence
            if "N" in sequence:
                logger.debug("Skipped domain %s with sequence %s. Likely a passivization", domain.id, domain.sequence)
            else:
                domain_temperatures.append(mt.Tm_NN(Seq(sequence), Na=c_na, Mg=c_mg))
           
        melting_temperatures.append(domain_temperatures)
    return [max(domain_ts) for domain_ts in melting_temperatures]


def compute_alpha_value(dna_structure: Structure, threshold:float, c_mg:float, c_na:float) -> float:
    """ calculates alpha value for a given critical temperature."""
    # NOTE: force recompute of domain data
    for strand in dna_structure.strands:
        strand.domain_list = []
    dna_structure._aux_data_computed = False
    dna_structure.compute_aux_data()

    domain_max_t = _max_melting_temperatures(dna_structure, c_mg=c_mg, c_na=c_na)

    domains_critical = sum(1 for melting_temp in domain_max_t if melting_temp >= threshold)
    return domains_critical / len(domain_max_t)


def permutate_optimal_scaffold_start(dna_structure: Structure, threshold:float, c_mg:float, c_na:float):
    scaffold = get_scaffold(dna_structure)
    logging.getLogger("nanodesign").setLevel(logging.WARNING)
    
    alpha_values = []
    for shift in range(len(scaffold.tour)):
        alpha = compute_alpha_value(dna_structure=dna_structure, threshold=threshold, c_mg=c_mg, c_na=c_na)
        alpha_values.append(alpha)
        start_base = scaffold.tour[shift]
        start = (start_base.h, start_base.p)
        logger.debug("Shift: %s-(%s). Alpha-%s value: %s", shift, start, threshold, alpha)
        shift_scaffold_sequence(scaffold=scaffold, step=1)

    max_alpha_value = max(alpha_values)
    max_index = alpha_values.index(max_alpha_value)

    start_base = scaffold.tour[max_index]
    start = (start_base.h, start_base.p)
    return (start, max_alpha_value)

def get_scaffold(dna_structure: Structure) -> Strand:
    scaffolds = [strand for strand in dna_structure.strands if strand.is_scaffold]
    if len(scaffolds) > 1:
        logger.critical("Found more than one scaffold. Abort.")
        sys.exit(1)
    return scaffolds[0]


def shift_scaffold_sequence(scaffold: Strand, step=1) -> None:
    """Shifts the scaffold sequence by step number of bases and recalculates sequences and auxillary data""" 
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}

    sequence = "".join(base.seq for base in scaffold.tour).upper()

    if step == 0:
        shifted_sequence = sequence
    else:
        shifted_sequence = sequence[-step:] + sequence[:-step] 

    for base in scaffold.tour:
        new_char, shifted_sequence = shifted_sequence[0], shifted_sequence[1:]
        base.seq = new_char
        if base.across is not None:
            base.across.seq = complement[new_char]
