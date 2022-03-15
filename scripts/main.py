#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full Apache License, Version 2.0 can be found in `LICENSE` at the project root.

import logging
import click

from pathlib import Path

from nanodesign.converters.converter import Converter
from nanodesign.converters.dna_sequence_data import dna_sequence_data

from .insert_uv_welding import welding


logger = logging.getLogger(__name__)


def print_version(ctx, _, value):
    """click print version."""
    if not value or ctx.resilient_parsing:
        return
    click.echo("1.0")  #TODO: add get_version
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