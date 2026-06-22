#!/usr/bin/env python


import numpy as np
from numpy import dtype

from beartype import beartype
from io import BufferedReader, BufferedWriter


@beartype
def bytes_to_int(input_bytes: bytes, signed:bool=False) -> int:
    """Returns the integer represented by the given array of bytes.
    Pre-sets the byteorder to be little-endian.

    Arguments:
        input_bytes -- Holds the array of bytes to convert.  The argument must either
            support the buffer protocol or be an iterable object producing bytes.
            Bytes and bytearray are examples of built-in objects that support the
            buffer protocol.

    Keyword Arguments:
        signed {bool} -- Indicates whether two's complement is used to represent the integer. (default: {False})

    Returns:
        [integer] -- Integer value converted from the supplied bytes.
    """
    return int.from_bytes(input_bytes, byteorder='little', signed=signed)


@beartype
def int_to_bytes(value: int, length: int) -> bytes:
    if length not in [1,2,4,8]:
        raise Exception("Should not happen")
    
    return value.to_bytes(length, 'little')


@beartype
def read_byte(fh_in: BufferedReader) -> int:
    """Converts a single byte to an integer value.

    Arguments:
        fh_in {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [integer] -- Unsigned integer value converted from the supplied bytes.
    """
    return bytes_to_int(fh_in.read(1), signed=False)


@beartype
def read_short(fh_in: BufferedReader) -> int:
    """Converts a two-byte element to an integer value.

    Arguments:
        fh_in {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [integer] -- Unsigned integer value converted from the supplied bytes.
    """
    return bytes_to_int(fh_in.read(2), signed=False)


@beartype
def read_int(fh_in: BufferedReader) -> int:
    """Converts a four-byte element to an integer value.

    Arguments:
        fh_in {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [integer] -- Signed integer value converted from the supplied bytes.
    """
    return bytes_to_int(fh_in.read(4), signed=False) # probably not need to be signed...


@beartype
def read_long(fh_in: BufferedReader) -> int:
    """Converts an eight-byte element to an integer value.

    Arguments:
        fh_in {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [integer] -- Signed integer value converted from the supplied bytes.
    """
    return bytes_to_int(fh_in.read(8), signed=False) # probably not need to be signed...


@beartype
def read_char(fh_in: BufferedReader, num_bytes: int) -> str:
    """Converts an array of bytes to a string.

    Arguments:
        fh_in {file-like} -- The binary file to read the select number of bytes.
        num_bytes {integer} -- The number of bytes to read and parse.

    Returns:
        [string] -- UTF-8 decoded string value.
    """
    return fh_in.read(num_bytes).decode('utf-8')


@beartype
def read_string(fh_in: BufferedReader) -> str:
    """Converts an array of bytes to a string.

    Arguments:
        fh_in {file-like} -- The binary file to read the select number of bytes.

    Returns:
        [string] -- UTF-8 decoded string value.
    """
    num_bytes = read_byte(fh_in)
    num_chars = num_bytes % 128
    shift = 0

    while num_bytes // 128 == 1:
        num_bytes = read_byte(fh_in)
        shift += 7
        offset = (num_bytes % 128) * (2 ** shift)
        num_chars += offset

    return read_char(fh_in, num_chars)



@beartype
def read_numpy_vector(fh_in: BufferedReader, dtype: dtype, n_elements: int):
    # https://stackoverflow.com/questions/72838939/how-to-convert-the-string-between-numpy-array-and-bytes
    """Parses a binary file multiple times, allowing for control if the
    file ends prematurely. This replaces read_results() and runs faster.
    And it provides support for reading gzipped idat files without decompressing.

    Arguments:
        fh_in {file-like} -- The binary file to read the select number of bytes.
        dtype {data type} -- used within idat files, 2-bit, or 4-bit numbers stored in binary at specific addresses
        n {number of snps read} -- see files/idat.py for how this function is applied.

    Raises:
        EOFError: If the end of the file is reached before the number of elements have
            been processed.

    Returns:
        A list of the parsed values.
    """
    dtype = np.dtype(dtype)
    
    alldata = fh_in.read(dtype.itemsize * n_elements)
    
    if len(alldata) != dtype.itemsize * n_elements:
        raise EOFError('End of file reached before number of results parsed')
    readdata=np.frombuffer(alldata, dtype, n_elements)
    if readdata.size != n_elements:
        raise EOFError('End of file reached before number of results parsed')

    return readdata


def write_numpy_vector(fh_out: BufferedReader, np_data):
    out_as_bytes = np.ndarray.tobytes(np_data)
    
    return fh_out.write(np_data)


@beartype
def write_char(fh_out: BufferedWriter, out: str) -> int:
    return fh_out.write(str.encode(out))

@beartype
def write_short(fh_out: BufferedWriter, out: int) -> int:
    return fh_out.write(int_to_bytes(out, 2))

@beartype
def write_int(fh_out: BufferedWriter, out: int) -> int:
    return fh_out.write(int_to_bytes(out, 4))

@beartype
def write_long(fh_out: BufferedWriter, out: int) -> int:
    return fh_out.write(int_to_bytes(out, 8))


@beartype
def binary_string_len(string: str) -> int:
    """
    0 ->     1 + 0
    1 ->     1  + 1
    127 ->   1 + 127
    
    128 ->   2 + 128
    16383 -> 2 + 16383
    
    16384 -> 3 + 16384
    ...
    """
    
    l = len(string)
    l_enc = long_to_7bit_string(l)
    
    return len(l_enc) + l



@beartype
def long_to_7bit_string(var_sized_integer: int) -> bytearray:
    out = bytearray(b'')
    
    offset = var_sized_integer % 128
    if var_sized_integer // 128 > 0:
        out.append(offset + 128)
    else:
        out.append(offset)

    while var_sized_integer // 128 > 0:
        var_sized_integer = var_sized_integer // 128
        offset = var_sized_integer % 128

        if var_sized_integer // 128 > 0:
            out.append(offset + 128)
        else:
            out.append(offset)
    
    return out



@beartype
def write_string(fh_out: BufferedWriter, out: str):
    return fh_out.write(long_to_7bit_string(len(out)) + str.encode(out))



