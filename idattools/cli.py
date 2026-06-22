#!/usr/bin/env python
# *- coding: utf-8 -*-

import click
from pathlib import Path

import idattools
from idattools.idat import IDATreader, IDATmixer

import pandas as pd


def main():
    CLI()


@click.version_option(
    idattools.__version__ + "\n\n" +
    idattools.__license_notice__ + "\n\nCopyright (C) 2024  " +
    idattools.__author__ + ".\n\nFor more info please visit:\n" +
    idattools.__homepage__
)
@click.group()
def CLI():
    pass


@CLI.command(name="view", short_help="View IDAT details (with [small] data summary)")
@click.argument('idat_file', type=click.Path(exists=True))
@click.option('-n', type=click.IntRange(min=1), default=10, help="Number of lines to print.", show_default=True)
@click.option('--no-header', is_flag=True, default=False,
              help="Print only the probe data table, omitting the metadata header.")
def CLI_view(idat_file, n, no_header):
    idat_r = IDATreader(Path(idat_file))

    pd.set_option('display.min_rows', n)
    pd.set_option('display.max_rows', n)
    pd.set_option('display.width', 240)
    pd.set_option('display.max_columns', 500)

    if no_header:
        print(str(idat_r.data.per_probe_matrix))
    else:
        print(str(idat_r.data))


@CLI.command(name="mix", short_help="Mix two IDAT files at a given ratio")
@click.argument('idat_file_reference', type=click.Path(exists=True))
@click.argument('idat_file_mixed_in', type=click.Path(exists=True))
@click.argument('idat_file_output', type=click.Path(exists=False))
@click.option('-r', '--mix-ratio', type=click.FloatRange(min=0, max=1), default=0.5,
              help="Fraction of mixed-in file values to be mixed into reference file. "
                   "E.g. 0.25 results in 75% of reference and 25% of mixed-in file.",
              show_default=True)
@click.option('--geometric-mean', is_flag=True, default=False,
              help="Mix intensities in log space (geometric mean: I_ref^(1-r) * I_mix^r) "
                   "instead of the default linear average. See README for rationale.")
def CLI_mix(idat_file_reference, idat_file_mixed_in, idat_file_output, mix_ratio, geometric_mean):
    idat_ref = IDATreader(Path(idat_file_reference))
    idat_mix = IDATreader(Path(idat_file_mixed_in))

    idattools.log.debug(
        "Mixing: " + idat_ref.data.get_sentrix_id() +
        " [" + str(round((1 - mix_ratio) * 100, 2)) + "%]" +
        " and " +
        idat_mix.data.get_sentrix_id() +
        " [" + str(round(mix_ratio * 100, 2)) + "%]" +
        (" [geometric mean]" if geometric_mean else " [linear]")
    )

    m = IDATmixer(idat_ref.data)
    m.mix(idat_mix.data, mix_ratio, Path(idat_file_output), geometric_mean=geometric_mean)
