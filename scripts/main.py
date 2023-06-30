#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full Apache License, Version 2.0 can be found in `LICENSE` at the project root.

import logging
import click

from pathlib import Path

from nanodesign.converters.converter import Converter
from nanodesign.converters.dna_sequence_data import dna_sequence_data
from nanodesign.converters.cadnano.reader import CadnanoReader
from nanodesign.converters.cadnano.writer import CadnanoWriter
from nanodesign.converters.cadnano.convert_design import CadnanoConvertDesign
from nanodesign.data.parameters import DnaParameters

from .insert_uv_welding import welding
from .stapler import Stapler
from .alpha_value import compute_alpha_value, permutate_optimal_scaffold_start


def init_logging():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler()
    formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    return logger


logger = init_logging()


def print_version(ctx, _, value):
    """click print version."""
    if not value or ctx.resilient_parsing:
        return
    click.echo("1.15.1")  #TODO: add get_version
    ctx.exit()


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-v",
    "--version",
    is_flag=True,
    help="Show __version__ and exit.",
    callback=print_version,
    expose_value=False,
    is_eager=True,
)
def cli():
    """nanodesign script interface.
    """
   

@cli.command()
@click.argument("design", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("sequence", type=str)
@click.option("-t", "--temperature",type=float, default=40.0, show_default=True, help="Alpha value threshold temperature in C.")
@click.option("-m", "--c-mg", "c_mg", type=float, default=20.0, show_default=True, help="Magnesium ion concentration in mM.")
@click.option("-n", "--c-ma", "c_na", type=float, default=5.0, show_default=True, help="Sodium ion concentration in mM.")
@click.option("-p", "--permutate", "search_start", is_flag=True, help="Find optimal scaffold start based on maximizing the alpha by sequence permutation.")
def alpha_value(design, sequence, temperature, c_mg, c_na, search_start):
    """\b
    Computes the alpha value at a given temperature 
    The alpha value is the percentage of staples that contain at least one domain with a melting temperature higher
        than the specified temperature

    Can find scaffold start postion with maximum alpha value.
    Note: only works for single scaffold designs
    \b
    DESIGN is the name of the design file [.json]
    SEQUENCE is the scaffold strand sequence file 
    """
    # parse command line arguments
    logger.info("Starting alpha value script")
    if Path(sequence).exists():
        seq_file_name = str(Path(sequence).absolute())
        seq_name = None
    elif sequence in dna_sequence_data:
        seq_file_name = None
        seq_name = sequence
    else:
        raise FileNotFoundError

    # Read cadnano file and create dna structure.
    converter = Converter()
    converter.modify = True  # NOTE: removes deletions and insertions.
    converter.read_cadnano_file(
        file_name=str(design),
        seq_file_name=seq_file_name,
        seq_name=seq_name,
    )
    dna_structure = converter.dna_structure

    alpha_value = compute_alpha_value(dna_structure, threshold=temperature, c_mg=c_mg, c_na=c_na)
    logger.info("Base alpha-%s value: %s", temperature, alpha_value)

    if search_start:
        optimal_start, optimal_alpha_value= permutate_optimal_scaffold_start(dna_structure, threshold=temperature, c_mg=c_mg, c_na=c_na)
        logger.info("Optimal scaffold starting position: %s,  alpha-%s value: %s", optimal_start, temperature, optimal_alpha_value)


@cli.command()
@click.argument("design", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("sequence", type=str)
@click.argument("fill_char", type=str)
def uv_welding(design, sequence, fill_char):
    """\b
    add additional fill characters to staple sequences.
    \b
    DESIGN is the name of the design file [.json]
    SEQUENCE is the scaffold strand sequence file 
    FILL_CHAR chararcter(s) to be used to signify inserted bases
    """

    # parse command line arguments
    logger.info("Starting UV Welding modifications")
    if Path(sequence).exists():
        seq_file_name = str(Path(sequence).absolute())
        seq_name = None
    elif sequence in dna_sequence_data:
        seq_file_name = None
        seq_name = sequence
    else:
        raise FileNotFoundError

    # Read cadnano file and create dna structure.
    converter = Converter()
    converter.modify = True # NOTE: removes deletions and insertions.
    converter.read_cadnano_file(
        file_name=str(design),
        seq_file_name=seq_file_name,
        seq_name=seq_name,
    )
    dna_structure = converter.dna_structure
    dna_structure.get_domains()

    output_lines = welding(dna_structure, fill_char=fill_char)

    with open("UVw_sequence.csv", mode="w") as seq_file:
        seq_file.writelines(output_lines)


@cli.command()
@click.argument("design", type=click.Path(exists=True, resolve_path=True, path_type=Path))
def stapler(design):
    """\b
    create new cadnano design with automated staple breaking.
    \b
    DESIGN is the name of the design file [.json]
    """

    # parse command line arguments
    logger.info("Starting stapler script")


    converter = Converter()
    cadnano_reader = CadnanoReader()
    converter.cadnano_design = cadnano_reader.read_json(design)

    dna_parameters = DnaParameters()
    converter.cadnano_convert_design = CadnanoConvertDesign(dna_parameters)
    converter.dna_structure = converter.cadnano_convert_design.create_structure(
        converter.cadnano_design
    )
    _ = converter.dna_structure.get_domains()


    stapler = Stapler(converter.dna_structure, converter.cadnano_design)
    stapler.template = [
        [7, 7, 14, 7, 7, 7],
        [7, 14, 7, 7, 7],
        [14, 7, 7, 7],
        [7, 7, 7, 14, 7, 7],
        [7, 7, 7, 14, 7],
        [7, 7, 7, 14],
        [7, 7, 14, 7, 7],
        [7, 14, 7, 7],
        [14, 7, 7],
        [7, 7, 14, 7],
        [7, 7, 14],
        [7, 14, 7],
    ]

    stapler.uncovered_domain_penalty = 10.0
    stapler.template_overhang_penalty_factor = 0.5
    stapler.step_probabilities = (0.1, 0.1, 0.35, 0.35, 0.1)
    stapler.standard_domain_length = 5

    stapler.initialize_system()

    # NOTE: For an overnight run, the paramaters (1e8, 1e4, 0.9999998047) should get to an error free structure.
    stapler.temperature = 10.0
    nr_steps = 10000000
    nr_steps_timescale = 10000
    stapler.generate(nr_steps, nr_steps_timescale, 0.999975)

    converter.dna_structure = converter.cadnano_convert_design.create_structure(
        converter.cadnano_design
    )
    converter.dna_structure.get_domains()
    logger.debug("Staple domains:")
    for strands in converter.dna_structure.strands[1:]:
        lengths = []
        for domain in strands.domain_list:
            lengths.append(len(domain.base_list))
        logger.debug(lengths)

    # Write a caDNAno JSON file.
    design_out = design.with_stem(f"{design.stem}_recode")
    logger.info("Write modified structure to file %s", design_out)
    cadnano_writer = CadnanoWriter(converter.dna_structure)
    cadnano_writer.write(design_out)