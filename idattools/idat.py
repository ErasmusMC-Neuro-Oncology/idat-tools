#!/usr/bin/env python

# https://github.com/FoxoTech/methylprep/blob/master/methylprep/files/idat.py#L228
# https://github.com/HenrikBengtsson/illuminaio/blob/develop/R/readIDAT_nonenc.R
# https://github.com/bioinformed/glu-genetics/blob/dcbbbf67a308d35e157b20a9c76373530510379a/glu/lib/illumina.py#L44-L61

import idattools
from .utils import *

from pathlib import Path
import os
import re
import random
import warnings

import numpy as np
import pandas as pd



section_names = {
    102: 'PROBE_IDS',
    103: 'PROBE_STD_DEVS',
    104: 'PROBE_MEAN_INTENSITIES',
    107: 'PROBE_N_BEADS',
    200: 'PROBE_MID_BLOCK', # also contains the probe_id's for some reason?
    300: 'ARRAY_RUN_INFO',
    400: 'ARRAY_RED_GREEN', # really concerned about this one, always [0]
    401: 'ARRAY_MANIFEST',
    402: 'ARRAY_BARCODE', # e.g. '203927450093'
    403: 'ARRAY_CHIP_TYPE', # e.g. 'BeadChip 8x5'
    404: 'ARRAY_CHIP_LABEL', # e.g. 'R01C01'
    405: 'ARRAY_OLD_STYLE_MANIFEST',
    406: 'ARRAY_SAMPLE_ID',
    407: 'ARRAY_DESCRIPTION',
    408: 'ARRAY_PLATE',
    409: 'ARRAY_WELL',
    410: 'ARRAY_UNKNOWN_1', # [1][0][0][0] <- could be int for 1, could be byte set, tuple of bytes is safest and easier to track new values down
    510: 'ARRAY_UNKNOWN_2',
    1000: 'ARRAY_N_PROBES'
}



class IDATdata(object):

    def __init__(self):
        self.file_magic = None
        self.idat_version = None
        self.section_index_order = None # the order of sections in index is typically different from implementation in the body
        self.section_physical_order = None        
        self.array_n_probes = None
        self.per_probe_matrix = None
        self.array_red_green = None
        self.array_manifest = None
        self.array_barcode = None
        self.array_chip_type = None
        self.array_chip_label = None
        self.array_old_style_manifest = None
        self.array_unknown_1 = None
        self.array_sample_id = None
        self.array_description = None
        self.array_plate = None
        self.array_well = None
        self.array_unknown_2 = None
        self.array_run_info: None


    def __str__(self):
        out = ""

        out += "# array_n_probes:       " + str(self.array_n_probes) + "\n"
        out += "# total intensity:      " + str(sum(self.per_probe_matrix['probe_mean_intensities'])) + "\n"
        out += "# manifest:             '" + str(self.array_manifest) + "'\n"
        out += "# manifest (old style): '" + str(self.array_old_style_manifest) + "'\n"
        out += "# unknown #1:           [" + "][".join([str(_) for _ in self.array_unknown_1]) + "]\n"
        out += "# sample id:            '" + str(self.array_sample_id) + "'\n"
        out += "# description:          '" + str(self.array_description) + "'\n"
        out += "# plate:                '" + str(self.array_plate) + "'\n"
        out += "# well:                 '" + str(self.array_well) + "'\n"
        out += "# unknown #2:           '" + str(self.array_unknown_2) + "'\n"
        out += "# run info:\n"
        for i in range(len(self.array_run_info)):
            out += "# "+str(i+1)+". [" + "] [".join([str(_) for _ in self.array_run_info[i]]) + "]\n"

        out += "\n"
        out += self.file_magic 
        out += " v"
        out += str(self.idat_version)
        out += ": "
        out += self.get_sentrix_id()
        out += " (R/G: "
        out += str(self.array_red_green)
        out += ", "
        out += self.array_chip_type
        out += ")"
        out += "\n"
        out += str(self.per_probe_matrix)

        return out


    @beartype
    def set_file_magic(self, file_magic: str) -> str:
        if file_magic != "IDAT":
            raise Exception("Invalid file format")
        else:
            self.file_magic = file_magic
        
        return self.file_magic


    @beartype
    def set_idat_version(self, idat_version: int) -> int:
        if idat_version != 3:
            raise Exception("Only tested with idat version 3")
        else:
            self.idat_version = idat_version
        
        return self.idat_version


    @beartype
    def set_array_n_probes(self, array_n_probes: int) -> int:
        if array_n_probes <= 0:
            raise Exception("Invalid number of probes: " + str(array_n_probes))
        else:
            self.array_n_probes = array_n_probes
        
        return self.array_n_probes


    @beartype
    def set_section_index_order(self, section_index_order: list[str]) -> list[str]:
        for _ in section_index_order:
            if _ not in section_names.values():
                raise Exception("Unknown section: "+str(_))
        
        self.section_index_order = section_index_order
        return self.section_index_order


    @beartype
    def set_section_physical_order(self, section_physical_order: list[str]) -> list[str]:
        for _ in section_physical_order:
            if _ not in section_names.values():
                raise Exception("Unknown section: "+str(_))
        
        self.section_physical_order = section_physical_order
        return self.section_physical_order


    @beartype
    def set_per_probe_matrix(self, per_probe_matrix: pd.DataFrame) -> pd.DataFrame:
        if per_probe_matrix.shape[0] != self.array_n_probes:
            raise Exception("Matrix (nrow: "+str(per_probe_matrix.shape[0])+") does no fit size of the array (n="+str(self.array_n_probes)+")")
    
        if not per_probe_matrix['probe_ids'].equals(per_probe_matrix['probe_mid_block']):
            raise Exception("Discrepance between probe_ids and probe_mid_block")

        if per_probe_matrix['probe_ids'].dtype != dtype("uint32"):
            raise Exception("Wrong format for probe_ids")

        if per_probe_matrix['probe_std_devs'].dtype != dtype("uint16"):
            raise Exception("Wrong format for probe_std_devs")

        if per_probe_matrix['probe_mean_intensities'].dtype != dtype("uint16"):
            raise Exception("Wrong format for probe_mean_intensities")

        if per_probe_matrix['probe_n_beads'].dtype != dtype("uint8"):
            raise Exception("Wrong format for probe_n_beads")

        if per_probe_matrix['probe_mid_block'].dtype != dtype("uint32"):
            raise Exception("Wrong format for probe_mid_block")


        self.per_probe_matrix = per_probe_matrix
        return self.per_probe_matrix


    @beartype
    def set_array_red_green(self, array_red_green: int) -> int:
        # checks here
        
        self.array_red_green = array_red_green
        return self.array_red_green


    @beartype
    def set_array_manifest(self, array_manifest: str) -> str:
        # checks here
        
        self.array_manifest = array_manifest
        return self.array_manifest


    @beartype
    def set_array_barcode(self, array_barcode: str) -> str:
        if not re.match(r"^[0-9]+$", array_barcode):
            raise Exception("Incorrect barcode: " + array_barcode)
        
        self.array_barcode = array_barcode
        return self.array_barcode


    @beartype
    def set_array_chip_type(self, array_chip_type: str) -> str:
        if array_chip_type not in ["BeadChip 8x5", "BeadChip 12x8"]:
            raise Exception("[%s] This tool is not tested with other chip types than BeadChip 8x5", array_chip_type)
        
        self.array_chip_type = array_chip_type
        return self.array_chip_type


    @beartype
    def set_array_chip_label(self, array_chip_label: str) -> str:
        if re.match(r"^R[0-9]+C[0-9]+_[0-9]$", array_chip_label):
            warnings.warn(f"Obscure but accepted array chip label: {array_chip_label}", UserWarning) # Seen happening in GSM8784882 and GSM8784883 ~ probably hacked and not original format ...
        elif not re.match(r"^R[0-9]+C[0-9]+$", array_chip_label):
            raise Exception(f"Odd label: {array_chip_label}")
        
        self.array_chip_label = array_chip_label
        return self.array_chip_label


    @beartype
    def set_array_old_style_manifest(self, array_old_style_manifest: str) -> str:
        # checks here
        
        self.array_old_style_manifest = array_old_style_manifest
        return self.array_old_style_manifest


    @beartype
    def set_array_unknown_1(self, array_unknown_1: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
        # checks here
        
        self.array_unknown_1 = array_unknown_1
        return self.array_unknown_1


    @beartype
    def set_array_sample_id(self, array_sample_id: str) -> str:
        # checks here
        
        self.array_sample_id = array_sample_id
        return self.array_sample_id


    @beartype
    def set_array_description(self, array_description: str) -> str:
        # checks here
        
        self.array_description = array_description
        return self.array_description


    @beartype
    def set_array_plate(self, array_plate: str) -> str:
        # checks here
        
        self.array_plate = array_plate
        return self.array_plate


    @beartype
    def set_array_well(self, array_well: str) -> str:
        # checks here
        
        self.array_well = array_well
        return self.array_well


    @beartype
    def set_array_unknown_2(self, array_unknown_2: str) -> str:
        # checks here
        
        self.array_unknown_2 = array_unknown_2
        return self.array_unknown_2


    @beartype
    def set_array_run_info(self, array_run_info: list[tuple[str, str, str, str, str]]) -> list[tuple[str, str, str, str, str]]:
        # checks here
        
        self.array_run_info = array_run_info
        return self.array_run_info

    def get_sentrix_id(self):
        return self.array_barcode + "_" + self.array_chip_label



class IDATreader:

    @beartype
    def __init__(self, idat_filename: Path):
        self.data = IDATdata()
        
        self.idat_filename = idat_filename
        self.parse()
    
    @beartype
    def parse_file_magic(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['FILE_MAGIC'])
        
        return self.data.set_file_magic(read_char(fh_in, 4))
    
    @beartype
    def parse_idat_version(self, fh_in: BufferedReader, section_seek_index: dict) -> int:
        fh_in.seek(section_seek_index['IDAT_VERSION'])
        
        return self.data.set_idat_version(read_long(fh_in))
        
    @beartype
    def parse_section_index(self, fh_in: BufferedReader, section_seek_index: dict) -> dict:
        section_index_order = []
        section_physical_order = {}

        fh_in.seek(section_seek_index['SECTION_INDEX_N'])
        n_sections = read_int(fh_in)

        for i in range(n_sections):
            section_type_int = read_short(fh_in)
            section_file_offset = read_long(fh_in)

            if section_type_int not in section_names:
                raise Exception("Unimplemented section type: "+str(section_type_int))
            else:
                section_type = section_names[section_type_int]
            
            section_index_order.append(section_type)
            section_physical_order[section_file_offset] = section_type
            section_seek_index[section_type] = section_file_offset

        self.data.set_section_index_order(section_index_order)
        self.data.set_section_physical_order([section_physical_order[_] for _ in sorted(section_physical_order.keys())])

        return section_seek_index
    
    @beartype
    def parse_array_n_probes(self, fh_in: BufferedReader, section_seek_index: dict) -> int:
        fh_in.seek(section_seek_index['ARRAY_N_PROBES'])
        
        return self.data.set_array_n_probes(read_int(fh_in))


    @beartype
    def parse_probe_ids(self, fh_in: BufferedReader, section_seek_index: dict) -> np.ndarray:
        fh_in.seek(section_seek_index['PROBE_IDS'])
        
        if self.data.array_n_probes is None:
            self.parse_array_n_probes(fh_in, section_seek_index)
        
        probe_ids = read_numpy_vector(fh_in, np.dtype('<u4'), self.data.array_n_probes) # layout-check: (4207470- 210) / 1051815 = 4
        
        if np.any(probe_ids <= 0):
            raise Exception("Wrong probe id's found")
        
        if np.any((probe_ids[1:] - probe_ids[:-1]) <= 0):
            raise Exception("probe id's are not unique or not incremental")
        
        return probe_ids

    @beartype
    def parse_probe_std_devs(self, fh_in: BufferedReader, section_seek_index: dict) -> np.ndarray:
        fh_in.seek(section_seek_index['PROBE_STD_DEVS'])
        
        if self.data.array_n_probes is None:
            self.parse_array_n_probes(fh_in, section_seek_index)
        
        probe_std_devs = read_numpy_vector(fh_in, np.dtype('<u2'), self.data.array_n_probes) # layout-check: (6311100 - 4207470) / 1051815 = 2
        
        if np.any(probe_std_devs < 0):
            raise Exception("Wrong std dev found (0 or negative)")
        
        return probe_std_devs

    @beartype
    def parse_probe_mean_intensities(self, fh_in: BufferedReader, section_seek_index: dict) -> np.ndarray:
        fh_in.seek(section_seek_index['PROBE_MEAN_INTENSITIES'])
        
        if self.data.array_n_probes is None:
            self.parse_array_n_probes(fh_in, section_seek_index)
        
        probe_mean_intensities = read_numpy_vector(fh_in, np.dtype('<u2'), self.data.array_n_probes) # layout-check: (8414730 - 6311100) / 1051815 = 2
        
        if np.any(probe_mean_intensities < 0):
            raise Exception("Wrong median probe intensity found (negative)")
        
        return probe_mean_intensities

    @beartype
    def parse_probe_n_beads(self, fh_in: BufferedReader, section_seek_index: dict) -> np.ndarray:
        fh_in.seek(section_seek_index['PROBE_N_BEADS'])
        
        if self.data.array_n_probes is None:
            self.parse_array_n_probes(fh_in, section_seek_index)
        
        probe_n_beads = read_numpy_vector(fh_in, np.dtype('<u1'), self.data.array_n_probes) # layout-check: (9466545 - 8414730) / 1051815 = 1
        
        if np.any(probe_n_beads < 0):
            raise Exception("Wrong number of beads per probe found (0 or negative)")
        
        return probe_n_beads

    @beartype
    def parse_probe_mid_block(self, fh_in: BufferedReader, section_seek_index: dict) -> np.ndarray:
        fh_in.seek(section_seek_index['PROBE_MID_BLOCK'])
        
        if self.data.array_n_probes is None:
            self.parse_array_n_probes(fh_in, section_seek_index)
        
        if self.data.array_n_probes != read_int(fh_in):
            raise Exception("Weird discrepancy between number of probes and size of mid block")
        
        probe_mid_block = read_numpy_vector(fh_in, np.dtype('<u4'), self.data.array_n_probes) # layout-check: (13673809 - (9466545 + 4)) / 1051815 = 4
        
        if np.any(probe_mid_block <= 0):
            raise Exception("Wrong probe id's found")
        
        if np.any((probe_mid_block[1:] - probe_mid_block[:-1]) <= 0):
            raise Exception("probe id's are not unique or not incremental")
        
        return probe_mid_block

    @beartype
    def parse_per_probe_matrix(self, fh_in: BufferedReader, section_seek_index: dict) -> pd.DataFrame:
        per_probe_matrix = pd.DataFrame({
            'probe_ids': self.parse_probe_ids(fh_in, section_seek_index),
            'probe_std_devs': self.parse_probe_std_devs(fh_in, section_seek_index),
            'probe_mean_intensities': self.parse_probe_mean_intensities(fh_in, section_seek_index),
            'probe_n_beads': self.parse_probe_n_beads(fh_in, section_seek_index),
            'probe_mid_block': self.parse_probe_mid_block(fh_in, section_seek_index)
            })

        if not per_probe_matrix['probe_ids'].equals(per_probe_matrix['probe_mid_block']):
            raise Exception("Discrepance between probe_ids and probe_mid_block")

        self.per_probe_matrix = per_probe_matrix

        return self.data.set_per_probe_matrix(per_probe_matrix)

    @beartype
    def parse_array_red_green(self, fh_in: BufferedReader, section_seek_index: dict) -> int:
        fh_in.seek(section_seek_index['ARRAY_RED_GREEN'])
        
        red_green = read_int(fh_in)
        
        #if red_green != 0:
        #    raise Exception("Only seen 0 so far, but probably good...")
        
        return self.data.set_array_red_green(red_green)

    @beartype
    def parse_array_manifest(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_MANIFEST'])
        
        return self.data.set_array_manifest(read_string(fh_in))
    
    @beartype
    def parse_array_barcode(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_BARCODE'])
        
        return self.data.set_array_barcode(read_string(fh_in))

    @beartype
    def parse_array_chip_type(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_CHIP_TYPE'])
        
        return self.data.set_array_chip_type(read_string(fh_in))

    @beartype
    def parse_array_chip_label(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_CHIP_LABEL'])
        
        try:
            return self.data.set_array_chip_label(read_string(fh_in))
        
        except Exception as e:
            raise Exception(f"File: {self.idat_filename} -- an error occurred: {e}")

    @beartype
    def parse_array_old_style_manifest(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_OLD_STYLE_MANIFEST'])
        
        return self.data.set_array_old_style_manifest(read_string(fh_in))

    @beartype
    def parse_array_unknown_1(self, fh_in: BufferedReader, section_seek_index: dict) -> tuple[int, int, int, int]:
        fh_in.seek(section_seek_index['ARRAY_UNKNOWN_1'])
        
        return self.data.set_array_unknown_1((read_byte(fh_in), read_byte(fh_in), read_byte(fh_in), read_byte(fh_in)))

    @beartype
    def parse_array_sample_id(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_SAMPLE_ID'])
        
        return self.data.set_array_sample_id(read_string(fh_in))

    @beartype
    def parse_array_description(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_DESCRIPTION'])
        
        return self.data.set_array_description(read_string(fh_in))

    @beartype
    def parse_array_plate(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_PLATE'])
        
        return self.data.set_array_plate(read_string(fh_in))

    @beartype
    def parse_array_well(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_WELL'])
        
        return self.data.set_array_well(read_string(fh_in))

    @beartype
    def parse_array_unknown_2(self, fh_in: BufferedReader, section_seek_index: dict) -> str:
        fh_in.seek(section_seek_index['ARRAY_UNKNOWN_2'])
        
        return self.data.set_array_unknown_2(read_string(fh_in))
    
    @beartype
    def parse_array_run_info(self, fh_in: BufferedReader, section_seek_index: dict) -> list[tuple[str, str, str, str, str]]:
        fh_in.seek(section_seek_index['ARRAY_RUN_INFO'])
        
        run_info = []
        
        for i in range(read_int(fh_in)): # blocks containing 5 consequtive strings
            run_info.append((read_string(fh_in), read_string(fh_in), read_string(fh_in), read_string(fh_in), read_string(fh_in)))
        
        return self.data.set_array_run_info(run_info)

    @beartype
    def parse(self) -> int:
        section_seek_index = { # static entries - todo: make class
            'FILE_MAGIC': 0,
            'IDAT_VERSION': 4,
            'SECTION_INDEX_N': 12
        }
        
        with open(self.idat_filename, "rb") as fh_in:
            self.parse_file_magic(fh_in, section_seek_index)
            self.parse_idat_version(fh_in, section_seek_index)
            self.parse_section_index(fh_in, section_seek_index)
            self.parse_array_n_probes(fh_in, section_seek_index)
            self.parse_per_probe_matrix(fh_in, section_seek_index)
            self.parse_array_red_green(fh_in, section_seek_index)
            self.parse_array_manifest(fh_in, section_seek_index)
            self.parse_array_barcode(fh_in, section_seek_index)
            self.parse_array_chip_type(fh_in, section_seek_index)
            self.parse_array_chip_label(fh_in, section_seek_index)
            self.parse_array_old_style_manifest(fh_in, section_seek_index)
            self.parse_array_unknown_1(fh_in, section_seek_index)
            self.parse_array_sample_id(fh_in, section_seek_index)
            self.parse_array_description(fh_in, section_seek_index)
            self.parse_array_plate(fh_in, section_seek_index)
            self.parse_array_well(fh_in, section_seek_index)
            self.parse_array_unknown_2(fh_in, section_seek_index)
            self.parse_array_run_info(fh_in, section_seek_index)

        return 0



class IDATwriter(IDATdata):

    @beartype
    def __init__(self, idat_data: IDATdata):
        if isinstance(idat_data, IDATdata):
            self.data = idat_data
        elif isinstance(idat_data, IDATreader):
            self.data = idat_data.data
        elif isinstance(idat_data, IDATwriter):
            self.data = idat_data.data
        else:
            raise Exception("Unclear input type")

    @beartype
    def write(self, idat_filename: Path):
        section_seek_index = { # static entries - todo: make class
            'FILE_MAGIC': 0,
            'IDAT_VERSION': 4,
            'SECTION_INDEX_N': 12
        }
        
        offset = 0

        with open(idat_filename, 'wb') as fh_out:
            offset += write_char(fh_out, self.data.file_magic)
            offset += write_long(fh_out, self.data.idat_version)
            offset += write_int(fh_out, len(self.data.section_index_order))

            section_sizes = {
                "ARRAY_N_PROBES": 4,
                "PROBE_IDS": (4 * self.data.array_n_probes),
                "PROBE_STD_DEVS": (2 * self.data.array_n_probes),
                "PROBE_MEAN_INTENSITIES": (2 * self.data.array_n_probes),
                "PROBE_N_BEADS": (1 * self.data.array_n_probes),
                "PROBE_MID_BLOCK": 4 + (4 * self.data.array_n_probes),
                "ARRAY_RUN_INFO": 4 + sum([ sum([binary_string_len(_) for _ in __ ]) for __ in self.data.array_run_info]),
                "ARRAY_RED_GREEN": 4,
                "ARRAY_MANIFEST": binary_string_len(self.data.array_manifest),
                "ARRAY_BARCODE": binary_string_len(self.data.array_barcode),
                "ARRAY_CHIP_TYPE": binary_string_len(self.data.array_chip_type),
                "ARRAY_CHIP_LABEL": binary_string_len(self.data.array_chip_label),
                "ARRAY_OLD_STYLE_MANIFEST": binary_string_len(self.data.array_old_style_manifest),
                "ARRAY_UNKNOWN_1": 1 + 1 + 1 + 1,
                "ARRAY_SAMPLE_ID": binary_string_len(self.data.array_sample_id),
                "ARRAY_DESCRIPTION": binary_string_len(self.data.array_sample_id),
                "ARRAY_PLATE": binary_string_len(self.data.array_plate),
                "ARRAY_WELL": binary_string_len(self.data.array_plate),
                "ARRAY_UNKNOWN_2": binary_string_len(self.data.array_plate)
            }

            offset_virtual = offset # should be 16
            offset_virtual += len(self.data.section_index_order) * (2 + 8)

            for section in self.data.section_index_order:
                section_code = [_ for _ in section_names.items() if _[1] == section][0][0]
                
                sections_before = self.data.section_physical_order[0:self.data.section_physical_order.index(section)]
                sections_before_sizes = [section_sizes[_] for _ in sections_before]

                offset_virtual_section = offset_virtual + sum(sections_before_sizes)
                
                offset += write_short(fh_out, section_code)
                offset += write_long(fh_out, offset_virtual_section)


            for section in self.data.section_physical_order: # keep original order of sections in file
                if section == "ARRAY_N_PROBES":
                    offset += write_int(fh_out, self.data.array_n_probes)
                elif section == "PROBE_IDS":
                    offset += write_numpy_vector(fh_out, self.data.per_probe_matrix["probe_ids"].to_numpy("<u4"))
                elif section == "PROBE_STD_DEVS":
                    offset += write_numpy_vector(fh_out, self.data.per_probe_matrix["probe_std_devs"].to_numpy("<u2"))
                elif section == "PROBE_MEAN_INTENSITIES":
                    offset += write_numpy_vector(fh_out, self.data.per_probe_matrix["probe_mean_intensities"].to_numpy("<u2"))
                elif section == "PROBE_N_BEADS":
                    offset += write_numpy_vector(fh_out, self.data.per_probe_matrix["probe_n_beads"].to_numpy("<u1"))
                elif section == "PROBE_MID_BLOCK":
                    offset += write_int(fh_out, self.data.array_n_probes)
                    offset += write_numpy_vector(fh_out, self.data.per_probe_matrix["probe_mid_block"].to_numpy("<u4"))
                elif section == "ARRAY_RED_GREEN":
                    offset += write_int(fh_out, self.data.array_red_green)
                elif section == "ARRAY_MANIFEST":
                    offset += write_string(fh_out, self.data.array_manifest)
                elif section == "ARRAY_BARCODE":
                    offset += write_string(fh_out, self.data.array_barcode)
                elif section == "ARRAY_CHIP_TYPE":
                    offset += write_string(fh_out, self.data.array_chip_type)
                elif section == "ARRAY_CHIP_LABEL":
                    offset += write_string(fh_out, self.data.array_chip_label)
                elif section == "ARRAY_OLD_STYLE_MANIFEST":
                    offset += write_string(fh_out, self.data.array_old_style_manifest)
                elif section == "ARRAY_UNKNOWN_1":
                    offset += write_char(fh_out, chr(self.data.array_unknown_1[0]))
                    offset += write_char(fh_out, chr(self.data.array_unknown_1[1]))
                    offset += write_char(fh_out, chr(self.data.array_unknown_1[2]))
                    offset += write_char(fh_out, chr(self.data.array_unknown_1[3]))
                elif section == "ARRAY_SAMPLE_ID":
                    offset += write_string(fh_out, self.data.array_sample_id)
                elif section == "ARRAY_DESCRIPTION":
                    offset += write_string(fh_out, self.data.array_description)
                elif section == "ARRAY_PLATE":
                    offset += write_string(fh_out, self.data.array_plate)
                elif section == "ARRAY_WELL":
                    offset += write_string(fh_out, self.data.array_well)
                elif section == "ARRAY_UNKNOWN_2":
                    offset += write_string(fh_out, self.data.array_unknown_2)
                elif section == "ARRAY_RUN_INFO":
                    n = len(self.data.array_run_info)
                    offset += write_int(fh_out, n)
                    for i in range(len(self.data.array_run_info)):
                        for j in range(5):
                            offset += write_string(fh_out, self.data.array_run_info[i][j])
                else:
                    raise Exception("Not implemented section: " + str(section))


class IDATmixer:
    @beartype
    def __init__(self, idat_reference: IDATdata):
        if isinstance(idat_reference, IDATdata):
            self.data_idat_ref = idat_reference
        elif isinstance(idat_reference, IDATreader):
            self.data_idat_ref = idat_reference.data
        else:
            raise Exception("Unclear input type (idat_reference)")

    @beartype
    def mix(self, idat_mixed_in: IDATdata, mixed_in_fraction: float, output_file: Path, geometric_mean: bool = False):
        if isinstance(idat_mixed_in, IDATdata):
            pass # ok
        elif isinstance(idat_mixed_in, IDATreader):
            idat_mixed_in = idat_mixed_in.data
        else:
            raise Exception("Unclear input type (idat_mixed_in)")

        idattools.log.debug("Initializing new IDATdata object")
        mixed_data = IDATdata()
        # check file_magic
        if self.data_idat_ref.file_magic != idat_mixed_in.file_magic:
            raise Exception("Different file magic's between reference and mixed-in sample - should never happen")
        else:
            mixed_data.set_file_magic(self.data_idat_ref.file_magic)

        # check idat_version
        if self.data_idat_ref.idat_version != idat_mixed_in.idat_version:
            raise Exception("Different idat versions between reference and mixed-in sample")
        else:
            mixed_data.set_idat_version(self.data_idat_ref.idat_version)

        if self.data_idat_ref.section_index_order !=  idat_mixed_in.section_index_order:
            idattools.log.warning("Order in which indexes are written differs, using the order from 'ref'")
        mixed_data.set_section_index_order(self.data_idat_ref.section_index_order)

        if self.data_idat_ref.section_physical_order !=  idat_mixed_in.section_physical_order:
            idattools.log.warning("Order in which sections are written differs, using the order from 'ref'")
        mixed_data.set_section_physical_order(self.data_idat_ref.section_physical_order)

        if self.data_idat_ref.array_red_green != idat_mixed_in.array_red_green:
            raise Exception("Differences in red_green value - odd, since this value seems to be always 0")
        else:
            mixed_data.set_array_red_green(self.data_idat_ref.array_red_green)

        if self.data_idat_ref.array_chip_type != idat_mixed_in.array_chip_type:
            raise Exception("Incompatible chip types: " + str(self.data_idat_ref.array_chip_type) + " and " + str(idat_mixed_in.array_chip_type))
        else:
            mixed_data.set_array_chip_type(self.data_idat_ref.array_chip_type)

        if self.data_idat_ref.array_manifest != idat_mixed_in.array_manifest:
            raise Exception("Differences in old_style_manifest value - odd, since this value seems to be always ''")
        else:
            mixed_data.set_array_manifest(self.data_idat_ref.array_manifest)

        if self.data_idat_ref.array_old_style_manifest != idat_mixed_in.array_old_style_manifest:
            raise Exception("Differences in old_style_manifest value - odd, since this value seems to be always ''")
        else:
            mixed_data.set_array_old_style_manifest(self.data_idat_ref.array_old_style_manifest)

        if self.data_idat_ref.array_unknown_1 != idat_mixed_in.array_unknown_1:
            raise Exception("Differences in unknown_1 value - odd, since this value seems to be always 4x0")
        else:
            mixed_data.set_array_unknown_1(self.data_idat_ref.array_unknown_1)

        if self.data_idat_ref.array_sample_id != idat_mixed_in.array_sample_id:
            raise Exception("Differences in sample_id value - odd, since this value seems to be always ''")
        else:
            mixed_data.set_array_sample_id(self.data_idat_ref.array_sample_id)

        if self.data_idat_ref.array_description != idat_mixed_in.array_description:
            raise Exception("Differences in description value - odd, since this value seems to be always ''")
        else:
            mixed_data.set_array_description(self.data_idat_ref.array_description)

        if self.data_idat_ref.array_plate != idat_mixed_in.array_plate:
            raise Exception("Differences in plate value - odd, since this value seems to be always ''")
        else:
            mixed_data.set_array_plate(self.data_idat_ref.array_plate)

        if self.data_idat_ref.array_well != idat_mixed_in.array_well:
            raise Exception("Differences in well value - odd, since this value seems to be always ''")
        else:
            mixed_data.set_array_well(self.data_idat_ref.array_well)

        if self.data_idat_ref.array_unknown_2 != idat_mixed_in.array_unknown_2:
            raise Exception("Differences in unknown_2 value - odd, since this value seems to be always ''")
        else:
            mixed_data.set_array_unknown_2(self.data_idat_ref.array_unknown_2)

        if re.match("^[0-9]{12}_R[0-9]{2}C[0-9]{2}.+idat$", os.path.basename(output_file)):
            barcode = os.path.basename(output_file).split("_")[0]
            chip_label = os.path.basename(output_file).split("_")[1][0:6]
        else:
            # @todo create some has of all the new data, and convert it to numeric weights (reproducible rather than random identifiers)
            barcode = "20"
            for i in range(10):
                barcode += str(random.randint(0,9))

            chip_label = "R0" + str(random.randint(1,8)) + "C01"

            idattools.log.warning("Output file does not comply with sentrix_id nomenclature ('012345678012_R0x_C0y.idat'), generating random one: "+str(barcode) + "_" + str(chip_label))

        mixed_data.set_array_barcode(barcode)
        mixed_data.set_array_chip_label(chip_label)
        idattools.log.debug("sentrix_id of mixed output file: " + mixed_data.get_sentrix_id())

        if np.any(self.data_idat_ref.per_probe_matrix.columns != idat_mixed_in.per_probe_matrix.columns):
            raise Exception("Different data columns in the arrays")


        data_left = self.data_idat_ref.per_probe_matrix
        data_right = idat_mixed_in.per_probe_matrix
        if self.data_idat_ref.array_n_probes != idat_mixed_in.array_n_probes:
            idattools.log.warning("Different sized arrays are merged ("+str(self.data_idat_ref.array_n_probes)+" ~ "+str(idat_mixed_in.array_n_probes)+") - reduing to intersect:")
            
            shared_probes = set(data_left["probe_ids"]).intersection(set(data_right["probe_ids"]))
            
            data_left = data_left[data_left["probe_ids"].isin(shared_probes)].reset_index(drop=True)
            data_right = data_right[data_right["probe_ids"].isin(shared_probes)].reset_index(drop=True)

            data_left_midblock = data_left[data_left["probe_ids"].isin(shared_probes)].reset_index(drop=True)
            data_right_midblock = data_right[data_right["probe_ids"].isin(shared_probes)].reset_index(drop=True)

            mixed_data.set_array_n_probes(len(shared_probes))
            
            idattools.log.warning("Size intersected array: " + str(len(shared_probes)) + " probes")
            
        else:
            mixed_data.set_array_n_probes(self.data_idat_ref.array_n_probes)

        if np.any(data_left["probe_ids"] != data_right["probe_ids"]):
            raise Exception("Arrays have different probe_ids (or ordering?)")
            
        if np.any(data_left["probe_mid_block"] != data_right["probe_mid_block"]):
            raise Exception("Arrays have different probe_mid_block id's (or ordering?)")


        f = mixed_in_fraction

        if geometric_mean:
            # Geometric mean: I_ref^(1-f) * I_mix^f = exp((1-f)*log(I_ref) + f*log(I_mix))
            # Intensities are clipped at 1 to avoid log(0) for probes with zero signal.
            left_i = np.clip(data_left["probe_mean_intensities"].to_numpy(float), 1, None)
            right_i = np.clip(data_right["probe_mean_intensities"].to_numpy(float), 1, None)
            mixed_intensities = np.round(np.exp((1 - f) * np.log(left_i) + f * np.log(right_i))).astype("<u2")

            left_sd = np.clip(data_left["probe_std_devs"].to_numpy(float), 1, None)
            right_sd = np.clip(data_right["probe_std_devs"].to_numpy(float), 1, None)
            mixed_std_devs = np.round(np.exp((1 - f) * np.log(left_sd) + f * np.log(right_sd))).astype("<u2")
        else:
            mixed_intensities = np.round(
                data_left["probe_mean_intensities"] * (1 - f) + data_right["probe_mean_intensities"] * f
            ).astype("<u2")
            mixed_std_devs = np.round(
                data_left["probe_std_devs"] * (1 - f) + data_right["probe_std_devs"] * f
            ).astype("<u2")

        new_data = pd.DataFrame({
            'probe_ids': data_left["probe_ids"],
            'probe_std_devs': mixed_std_devs,
            'probe_mean_intensities': mixed_intensities,
            'probe_n_beads': np.round(
                data_left["probe_n_beads"] * (1 - f) + data_right["probe_n_beads"] * f
            ).astype("<u1"),
            'probe_mid_block': data_left["probe_mid_block"]
        })

        mixed_data.set_per_probe_matrix(new_data)



        if len(self.data_idat_ref.array_run_info) != len(idat_mixed_in.array_run_info):
            raise Exception("different array_run_info size")
        else:
            ri = []
            for i in range(len(self.data_idat_ref.array_run_info)):
                rir = (
                    self.data_idat_ref.get_sentrix_id() + ":" + self.data_idat_ref.array_run_info[i][0] + "&" + \
                    idat_mixed_in.get_sentrix_id() + ":" + idat_mixed_in.array_run_info[i][0],
                    self.data_idat_ref.get_sentrix_id() + ":" + self.data_idat_ref.array_run_info[i][1] + "&" + \
                    idat_mixed_in.get_sentrix_id() + ":" + idat_mixed_in.array_run_info[i][1],
                    self.data_idat_ref.get_sentrix_id() + ":" + self.data_idat_ref.array_run_info[i][2] + "&" + \
                    idat_mixed_in.get_sentrix_id() + ":" + idat_mixed_in.array_run_info[i][2],
                    self.data_idat_ref.get_sentrix_id() + ":" + self.data_idat_ref.array_run_info[i][3] + "&" + \
                    idat_mixed_in.get_sentrix_id() + ":" + idat_mixed_in.array_run_info[i][3],
                    self.data_idat_ref.get_sentrix_id() + ":" + self.data_idat_ref.array_run_info[i][4] + "&" + \
                    idat_mixed_in.get_sentrix_id() + ":" + idat_mixed_in.array_run_info[i][4]
                )

                ri.append(rir)
            
        mixed_data.set_array_run_info(ri)
        
        w = IDATwriter(mixed_data)
        w.write(output_file)
        
        return mixed_data




