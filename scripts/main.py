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
    click.echo("1.0.1")  #TODO: add get_version
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
@click.option("-t", "--temperature", "temperature", type=float, default=10.0, show_default=True, help="Temperature in C.")
@click.option("-r", "--temperature-rate", "temperature_rate", type=float, default=0.999975, show_default=True, help="Temperature adjustment rate.")
@click.option("-n", "--steps", "steps", type=int, default=int(1e7), show_default=True, help="Number of steps.")
@click.option("-m", "--steps-timescale", "steps_timescale", type=int, default=int(1e4), show_default=True, help="Number of timescale steps.")

def stapler(design, steps, temperature, steps_timescale, temperature_rate):
    """\b
    create new cadnano design with automated staple breaking.
    \b
    DESIGN is the name of the design file [.json]
    """

    logger.info("Starting stapler script")
    converter = Converter()
    converter.read_cadnano_file(
                file_name=design,
                seq_file_name=None,  # TODO
                seq_name=None,
            )
    converter.dna_structure.compute_aux_data()

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

    # NOTE: general parameter set
    stapler.uncovered_domain_penalty = 10.0
    stapler.template_overhang_penalty_factor = 0.5
    stapler.step_probabilities = (0.1, 0.1, 0.35, 0.35, 0.1)
    stapler.standard_domain_length = 5

    stapler.initialize_system()
    # NOTE: For an overnight run, the parameters (1e8, 1e4, 0.9999998047) should get to an error free structure.
    stapler.temperature = temperature
    design = stapler.generate(nr_steps=steps, nr_steps_timescale=steps_timescale, temperature_adjust_rate=temperature_rate)
    # TODO: Check following edge condition. Reported by JP. Get sample file? (JMS 11/21/16)
    #    structure with only two staples with two segments each, one double crossover, crashes:
    #    **** ERROR: Reached a visited base.Traceback (most recent call last):
    #    File "testScript3.py", line 357, in <module>
    #    main()
    #    File "testScript3.py", line 346, in main
    #    converter.dna_structure = converter.cadnano_convert_design.create_structure(converter.cadnano_design)
    #    File "/converters/cadnano/convert_design.py", line 170, in create_structure
    #    self._set_strands_colors(strands)
    #    File "/converters/cadnano/convert_design.py", line 1124, in _set_strands_colors
    #    for strand in strands:
    #    TypeError: 'NoneType' object is not iterable
    
    # regenerate structure from modified design
    dna_structure = converter.cadnano_convert_design.create_structure(design)
    dna_structure.compute_aux_data()

    logger.debug("Staple domains:")
    for strand in dna_structure.strands[1:]:
        logger.debug("Strand %s", strand.id, [len(domain.base_list) for domain in strand.domain_list])

    # Write a caDNAno JSON file.
    design_out = design.with_stem(f"{design.stem}_recode")
    logger.info("Write modified structure to file %s", design_out)
    cadnano_writer = CadnanoWriter(dna_structure)
    cadnano_writer.write(design_out)